#include <string.h>
#include "we_op.hpp"
#include "bsplines.hpp"
#include "nlsolver.hpp"
#include "IO.hpp"
#include "mpiWrapper.hpp"
#include "seplib.h"

typedef vecReg<data_t> vec;
typedef cvecReg<data_t> cvec;
typedef axis<data_t> ax;
typedef hypercube<data_t> hyper;

// Executable to run 2D FWI

int main(int argc, char **argv){

    int rank=0, size=1;
    mpiWrapper::init(&argc,&argv);
    mpiWrapper::setSizeRank(&size,&rank);
    fprintf (stderr,"\n====================\nSize of MPI communicator = %d ; current rank = %d\n====================\n",size,rank);
    
    initpar(argc,argv);

// Read parameters for wave propagation and inversion
    param par;
    readParameters(argc, argv, par);
    int verbose=par.verbose;
    if (rank>0) par.verbose=0;

    time_t t1 = time(NULL);
    if (par.verbose>0) fprintf(stderr,"\n====================\n%s\n====================\n",ctime(&t1));

// Set the maximum number of threads
    if (par.nthreads>0) omp_set_num_threads(par.nthreads);

// Read inputs/outputs files
    std::string source_file="none", model_file="none", bsmodel_file="none", data_file="none", output_file="none", ioutput_file="none", obj_func_file="none";
    readParam<std::string>(argc, argv, "source", source_file);
    readParam<std::string>(argc, argv, "model", model_file);
    readParam<std::string>(argc, argv, "bsmodel", bsmodel_file);
    readParam<std::string>(argc, argv, "data", data_file);
    readParam<std::string>(argc, argv, "output", output_file);
    readParam<std::string>(argc, argv, "ioutput", ioutput_file);
    readParam<std::string>(argc, argv, "obj_func", obj_func_file);

    successCheck(source_file!="none","Source wavelet is not provided\n");
    successCheck(model_file!="none","Earth model is not provided\n");
    successCheck(data_file!="none","Data to be inverted is not provided\n");

    std::shared_ptr<vec> src = read<data_t>(source_file, par.format);
    std::shared_ptr<vec> data = read<data_t>(data_file, par.format);
    std::shared_ptr<vec> model = read<data_t>(model_file, par.format);

    successCheck( src->getHyper()->getAxis(1).n==data->getHyper()->getAxis(1).n
                && src->getHyper()->getAxis(1).d==data->getHyper()->getAxis(1).d,
                "source wavelet(s) and data must have the same number of time samples and sampling rate\n" );

    std::shared_ptr<vec> hrz = nullptr;
    std::shared_ptr<vec> gmask = nullptr;
    std::shared_ptr<vec> w = nullptr;
    std::shared_ptr<vec> filter = nullptr;
    std::shared_ptr<vec> invDiagH = nullptr;
    std::shared_ptr<vec> prior = nullptr;
    if (par.mask_file!="none") {gmask = read<data_t>(par.mask_file, par.format); successCheck(gmask->getN123()==model->getN123(),"Gradient mask must have the same number of samples as the model\n");}
    if (par.weights_file!="none") {w = read<data_t>(par.weights_file, par.format); successCheck(w->getN123()==data->getN123(),"Data weights must have the same number of samples as the data\n");}
    if (par.filter_file!="none") {filter = read<data_t>(par.filter_file, par.format); successCheck(filter->getHyper()->getAxis(1).d==data->getHyper()->getAxis(1).d,"Filter and data must have the same sampling rate\n");}
    if (par.inverse_diagonal_hessian_file!="none") {invDiagH = read<data_t>(par.inverse_diagonal_hessian_file, par.format); successCheck(invDiagH->getN123()==model->getN123(),"Inverse diagonal Hessian must have the same number of samples as the model\n");}
    if (par.prior_file!="none") {prior = read<data_t>(par.prior_file, par.format); successCheck(prior->getN123()==model->getN123(),"Prior model must have the same number of samples as the model\n");}

// Analyze the inputs and parameters and modify if necessary
    analyzeGeometry(*model->getHyper(),par, par.verbose>0);
    std::shared_ptr<vec> allsrc = analyzeWavelet(src, par, par.verbose>0);
    analyzeBsplines(*model->getHyper(),par);
    analyzeNLInversion(par);

    ax axR(par.nr*par.nrcomp,0,1);
    hyper range(src->getHyper()->getAxis(1), axR);
    data->setHyper(range);

// Build model precon if 1D inversion is activated
// ----------------------------------------------------------------------------------------//
    std::shared_ptr<vec> model1d = model;
    std::shared_ptr<vec> gmask1d = gmask;
    std::shared_ptr<vec> invDiagH1d = invDiagH;
    std::shared_ptr<vec> prior1d = prior;
    loper * ext1d;

if (par.inversion1d)
{
    if (par.horizon_file!="none") {
        hrz = read<data_t>(par.horizon_file, par.format);
        successCheck(hrz->getHyper()->getNdim()==2,"The provided horizon must have two dimensions (x and y)\n");
        successCheck(hrz->getHyper()->getAxis(1).n==model->getHyper()->getAxis(2).n,"The provided horizon must have the same size as the x-dimension of the model\n");
        successCheck(hrz->getHyper()->getAxis(2).n==model->getHyper()->getAxis(3).n,"The provided horizon must have the same size as the y-dimension of the model\n");
    }
    else{
        hrz = std::make_shared<vec>(hyper(model->getHyper()->getAxis(2),model->getHyper()->getAxis(3)));
        hrz->zero();
    }
    std::vector<ax> axes = model->getHyper()->getAxes();
    axes.erase(axes.begin()+1); // remove x-axis
    axes.erase(axes.begin()+1); // remove y-axis
    model1d = std::make_shared<vec> (hyper(axes));
    model1d->zero();
    ext1d = new extrapolator1d3d(*model1d->getHyper(),hrz);
    ext1d->inverse(false,model1d,model);

    if (par.mask_file != "none") {
        gmask1d = std::make_shared<vec> (hyper(axes));
        gmask1d->zero();
        ext1d->inverse(false,gmask1d,gmask);
    }
    if (par.inverse_diagonal_hessian_file != "none"){
        invDiagH1d = std::make_shared<vec> (hyper(axes));
        invDiagH1d->zero();
        ext1d->inverse(false,invDiagH1d,invDiagH);
    }
    if (par.prior_file != "none"){
        prior1d = std::make_shared<vec> (hyper(axes));
        prior1d->zero();
        ext1d->inverse(false,prior1d,prior);
    }
}
// ----------------------------------------------------------------------------------------//

// Build model parameterization precon
// ----------------------------------------------------------------------------------------//
    std::shared_ptr<vec> model_temp = model1d;
    std::shared_ptr<vec> prior_temp = prior1d;
    int n = model1d->getN123()/par.nmodels;
    data_t vs0=0;
    data_t rho0=0;
    if (par.nmodels>=3){
        vs0 = model1d->sum(n, 2*n);
        rho0 = model1d->sum(2*n, 3*n);
        vs0 /= n;
        rho0 /= n;
    }

    nloper * P = nullptr;
    if (par.model_parameterization==0) P = new lam_mu_rho(*model1d->getHyper());
    else if (par.model_parameterization==2) P = new ip_is_rho(*model1d->getHyper());
    else if (par.model_parameterization==3) {
        P = new vs_vpvs_rho(*model1d->getHyper(), vs0, rho0);
        if (par.verbose>0) fprintf(stderr,"Average Vs and rho used in model parameterization are %.3f and %.3f\n",vs0/n, rho0/n);
    }
    if (P != nullptr){
        model_temp = std::make_shared<vec>(*model1d->getHyper());
        model_temp->zero();
        P->inverse(false,model_temp,model1d);
        if (par.prior_file != "none") {
            prior_temp = std::make_shared<vec>(*model1d->getHyper());
            prior_temp->zero();
            P->inverse(false,prior_temp,prior1d);
        }
    }
    model1d = model_temp;
    prior1d = prior_temp;
// ----------------------------------------------------------------------------------------//

// Build model precon if B-splines are activated
// ----------------------------------------------------------------------------------------//
    std::shared_ptr<vec> bsmodel = model1d;
    std::shared_ptr<vec> bsmask = gmask1d;
    std::shared_ptr<vec> bsinvDiagH = invDiagH1d;
    std::shared_ptr<vec> bsprior = prior1d;
    loper * BD;

if (par.bsplines)
{
    std::vector<data_t> kx;
    std::vector<data_t> ky;
    std::vector<data_t> kz;
    setKnot(kx,par.bs_controlx,par.bs_mx);
    setKnot(ky,par.bs_controly,par.bs_my);
    setKnot(kz,par.bs_controlz,par.bs_mz);

    std::vector<ax > axes = model->getHyper()->getAxes();
    axes[0].n=par.bs_nz; axes[1].n=par.bs_nx; axes[2].n=par.bs_ny;
    if (par.inversion1d) {
        axes.erase(axes.begin()+1);
        axes.erase(axes.begin()+1);
    }
    bsmodel = std::make_shared<vec>(vec(hyper(axes)));
    if (bsmodel_file!="none") {
        bsmodel = read<data_t>(bsmodel_file, par.format);
        successCheck(bsmodel->getHyper()->isCompatible(hyper(axes)),"The B-splines model is not compatible with the input and B-splines parameters\n");
    }
    else if (par.inversion1d) fillin1d(bsmodel,model1d,par.bs_controlz);
    else fillin3d(bsmodel,model1d,par.bs_controlx,par.bs_controly,par.bs_controlz);
    
    if (par.mask_file != "none"){
        bsmask = std::make_shared<vec>(vec(hyper(axes)));
        if (par.inversion1d) fillin1d(bsmask,gmask1d,par.bs_controlz);
        else fillin3d(bsmask,gmask1d,par.bs_controlx,par.bs_controly,par.bs_controlz);
    }
    if (par.inverse_diagonal_hessian_file != "none"){
        bsinvDiagH = std::make_shared<vec>(vec(hyper(axes)));
        if (par.inversion1d) fillin1d(bsinvDiagH,invDiagH1d,par.bs_controlz);
        else fillin3d(bsinvDiagH,invDiagH1d,par.bs_controlx,par.bs_controly,par.bs_controlz);
    }
    if (par.prior_file != "none"){
        bsprior = std::make_shared<vec>(vec(hyper(axes)));
        if (par.inversion1d) fillin1d(bsprior,prior1d,par.bs_controlz);
        else fillin3d(bsprior,prior1d,par.bs_controlx,par.bs_controly,par.bs_controlz);
    }

    if (par.inversion1d){
        duplicate1d D(*bsmodel->getHyper(),par.bs_mz);
        bsplines1d B(*D.getRange(),*model1d->getHyper(),kz);
        BD = new chainLOper(&B,&D);
    }
    else{
        duplicate3d D(*bsmodel->getHyper(),par.bs_mx,par.bs_my,par.bs_mz);
        bsplines3d B(*D.getRange(),*model1d->getHyper(),kx,ky,kz);
        BD = new chainLOper(&B,&D);
    }
}    
// ----------------------------------------------------------------------------------------//

// Build model precon if soft clipping is activated
// ----------------------------------------------------------------------------------------//
    emodelSoftClip * S;
    if (par.soft_clip) 
    {
        if (par.model_parameterization==0) {
            data_t mu_min = par.rhomin*par.vsmin*par.vsmin;
            data_t mu_max = par.rhomax*par.vsmax*par.vsmax;
            data_t lam_min = par.rhomin*(par.vpmin*par.vpmin - 2*par.vsmax*par.vsmax);
            data_t lam_max = par.rhomax*(par.vpmax*par.vpmax - 2*par.vsmin*par.vsmin);
            lam_min = std::max((data_t)0.0 , lam_min);
            S = new emodelSoftClip(*model1d->getHyper(), lam_min, lam_max, mu_min, mu_max, par.rhomin, par.rhomax, 1, 9, 9);
        }
        else if (par.model_parameterization==2){
            data_t ip_min = par.rhomin*par.vpmin;
            data_t ip_max = par.rhomax*par.vpmax;
            data_t is_min = par.rhomin*par.vsmin;
            data_t is_max = par.rhomax*par.vsmax;
            S = new emodelSoftClip(*model1d->getHyper(), ip_min, ip_max, is_min, is_max, par.rhomin, par.rhomax, 1/sqrt(2.00001), 9, 9);
        }
        else if (par.model_parameterization==3) S = new emodelSoftClip(*model1d->getHyper(), log(par.vsmin/vs0), log(par.vsmax/vs0), -10, log(par.vpmax/par.vsmin - sqrt(2)), log(par.rhomin/rho0), log(par.rhomax/rho0), 1, 9, 9);
        else S = new emodelSoftClip(*model1d->getHyper(), par.vpmin, par.vpmax, par.vsmin, par.vsmax, par.rhomin, par.rhomax, 1/sqrt(2.00001), 9, 9);
    }
// ----------------------------------------------------------------------------------------//

// Build the appropriate operators
// ----------------------------------------------------------------------------------------//
    
    nloper * op = nullptr; // full model preconditioner
    nl_we_op * L; // wave modeling operator
    if (par.nmodels==2) L=new nl_we_op_a(*model->getHyper(),allsrc,0,par);
    else if (par.nmodels==3) L=new nl_we_op_e(*model->getHyper(),allsrc,0,par);
    else if (par.nmodels==6) L=new nl_we_op_vti(*model->getHyper(),allsrc,0,par);
    else successCheck(false, "Number of paramaters in the model is not supported\n");

    if (par.bsplines)
    {
        if (P != nullptr)
        {
            if (par.inversion1d)
            {
                if (par.soft_clip)
                {
                    chainNLOper SBD(S,BD);
                    chainNLOper PSBD(P,&SBD);
                    op = new chainNLOper(ext1d,&PSBD);
                }
                else
                {
                    chainNLOper PBD(P,BD);
                    op = new chainNLOper(ext1d,&PBD);
                }
            }
            else
            {
                if (par.soft_clip)
                {
                    chainNLOper SBD(S,BD);
                    op = new chainNLOper(P,&SBD);
                }
                else op = new chainNLOper(P,BD);
            }
        }
        else
        {
            if (par.inversion1d)
            {
                if (par.soft_clip)
                {
                    chainNLOper SBD(S,BD);
                    op = new chainNLOper(ext1d,&SBD);
                }
                else op = new chainNLOper(ext1d,BD);
            }
            else
            {
                if (par.soft_clip) op = new chainNLOper(S,BD);
                else op  = BD->clone();
            }
        }
    }
    else
    {
        if (P != nullptr)
        {
            if (par.inversion1d)
            {
                if (par.soft_clip)
                {
                    chainNLOper PS(P,S);
                    op = new chainNLOper(ext1d,&PS);
                }
                else op = new chainNLOper(ext1d,P);
            }
            else
            {
                if (par.soft_clip) op = new chainNLOper(P,S);
                else op = P->clone();
            }
        }
        else
        {
            if (par.inversion1d)
            {
                if (par.soft_clip) op = new chainNLOper(ext1d,S);
                else op = ext1d->clone();
            }
            else if (par.soft_clip) op = S->clone();
        }
    }

    if (par.bsplines) delete BD;
    if (P != nullptr) delete P;
    if (par.inversion1d) delete ext1d;
    if (par.soft_clip) delete S;

    nloper * D = nullptr;
    loper * R;
    if (par.regularization==0) R = new identity (*model->getHyper());
    else if (par.regularization==1) R = new gradient3d(*model->getHyper(),par.reg_xweight,par.reg_yweight,par.reg_zweight);
    else if (par.regularization==2) R = new laplacian3d(*model->getHyper(),par.reg_xweight,par.reg_yweight,par.reg_zweight);
    else R = nullptr;
    if (R!=nullptr){
        if (op==nullptr) D = R->clone();
        else D = new chainNLOper(R, op);
        delete R;
    }
// ----------------------------------------------------------------------------------------//

// Construct the optimization problem, the solver, then solve the problem
// ----------------------------------------------------------------------------------------//

    par.verbose = verbose;

    nlls_fwi * prob;
    if (D != nullptr) {
        prob = new nlls_fwi_reg(L, par, D, bsmodel, data, par.lambda, bsprior, op, bsmask, w, filter);
    }
    else prob = new nlls_fwi(L, par, bsmodel, data, op, bsmask, w, filter);

    if (rank>0) par.verbose=0;

    lsearch * ls;
    if (par.lsearch=="weak_wolfe") ls = new weak_wolfe(par.ls_c1, par.ls_a0, par.ls_a1, par.ls_version);
    else if(par.lsearch=="strong_wolfe") ls = new strong_wolfe(par.ls_c1, par.ls_c2, par.ls_a0, par.ls_a1, par.ls_max_step, par.ls_version);
    else ls = new regular_wolfe(par.ls_c1, par.ls_c2, par.ls_a0, par.ls_a1, par.ls_max_step, par.ls_version);

    nlsolver * solver;
    if (par.nlsolver=="nlsd") solver = new nlsd(par.niter, par.max_trial, par.threshold, ls); 
    else if (par.nlsolver=="nlcg") solver = new nlcg(par.niter, par.max_trial, par.threshold, ls); 
    else if (par.nlsolver=="bfgs") solver = new bfgs(par.niter, par.max_trial, par.threshold, ls); 
    else solver = new lbfgs(par.niter, par.max_trial, par.threshold, ls, bsinvDiagH, par.lbfgs_m); 
    
    solver->run(prob, par.verbose>0, ioutput_file, par.isave, par.format, par.datapath);
    model=bsmodel;
// ----------------------------------------------------------------------------------------//

// Clean up and save outputs
// ----------------------------------------------------------------------------------------//
    
    if (D != nullptr) delete D;
    if (op != nullptr)
    {
        std::shared_ptr<vec> tmp = std::make_shared<vec>(*op->getRange());
        tmp->zero();
        op->forward(false, model, tmp);
        model = tmp;
        delete op;
    }

    if (rank==0 && output_file!="none") write<data_t>(model, output_file, par.format, par.datapath);
    if (rank==0 && obj_func_file!="none") {
        std::shared_ptr<vec > func = std::make_shared<vec > (hyper(solver->_func.size()));
        memcpy(func->getVals(), solver->_func.data(), solver->_func.size()*sizeof(data_t));
        write<data_t>(func,obj_func_file, par.format, par.datapath);
        if (D != nullptr){
            std::shared_ptr<vec > dfunc = std::make_shared<vec > (hyper(prob->_dfunc.size()));
            std::shared_ptr<vec > mfunc = std::make_shared<vec > (hyper(prob->_mfunc.size()));
            memcpy(dfunc->getVals(), prob->_dfunc.data(), prob->_dfunc.size()*sizeof(data_t));
            memcpy(mfunc->getVals(), prob->_mfunc.data(), prob->_mfunc.size()*sizeof(data_t));
            write<data_t>(dfunc,obj_func_file+".d", par.format, par.datapath);
            write<data_t>(mfunc,obj_func_file+".m", par.format, par.datapath);
        }
    }

    delete L;
    delete prob;
// ----------------------------------------------------------------------------------------//

    mpiWrapper::finalize();

    time_t t2 = time(NULL);
    if (par.verbose>0) fprintf(stderr,"\n====================\n%s\n====================\n",ctime(&t2));

return 0;
}