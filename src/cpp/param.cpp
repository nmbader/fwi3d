#include "param.hpp"
#include "misc.hpp"

void readCoord(std::string srcoord, param &par){

    std::ifstream ifs (srcoord, std::ifstream::in);
    if (ifs.is_open())
    {
        par.sxyz.clear();
        par.rxyz.clear();
        par.srcoord_from_file=true;
        std::string line;
        int counter=0;
        int shot_counter=0;
        int curr=0, prev=0;
        data_t sid = 0;
        while ( getline (ifs, line) )
        {
            if (line[0] != '#' && line[0] != '//')
            {
                std::istringstream ss(line);

                ss >> sid; // read shot id and convert to integer
                curr = (int)sid;
                if (counter==0) // first line
                {
                    // read shot x y z coordinates
                    std::vector<data_t> sxyz(3);
                    ss >> sxyz[0];
                    ss >> sxyz[1];
                    ss >> sxyz[2];

                    // read receiver x y z dip azimuth coordinates
                    std::vector<data_t> rxyz(5);
                    ss >> rxyz[0];
                    ss >> rxyz[1];
                    ss >> rxyz[2];
                    ss >> rxyz[3];
                    ss >> rxyz[4];

                    par.sxyz.push_back(sxyz);
                    std::vector<std::vector<data_t> > vec = {rxyz};
                    par.rxyz.push_back(vec);

                    counter++;
                    prev=curr;
                }
                else // subsequent lines
                {
                    // read shot x z coordinates
                    std::vector<data_t> sxyz(3);
                    ss >> sxyz[0];
                    ss >> sxyz[1];
                    ss >> sxyz[2];

                    // read receiver x y z dip azimuth coordinates
                    std::vector<data_t> rxyz(5);
                    ss >> rxyz[0];
                    ss >> rxyz[1];
                    ss >> rxyz[2];
                    ss >> rxyz[3];
                    ss >> rxyz[4];

                    if (curr==prev)
                    {
                        par.rxyz[shot_counter].push_back(rxyz);
                    }
                    else
                    {   
                        par.sxyz.push_back(sxyz);
                        std::vector<std::vector<data_t> > vec = {rxyz};
                        par.rxyz.push_back(vec);
                        shot_counter++;
                    }
                    counter++;
                    prev=curr;

                }
            }
        }
        ifs.close();
    }
}

void readParam(int argc, char **argv, param &par){
    readParam<std::string>(argc, argv, "resampling", par.resampling);
    readParam<std::string>(argc, argv, "srcoord", par.srcoord);
    readParam<std::string>(argc, argv, "horizon", par.horizon_file);
    readParam<std::string>(argc, argv, "nlsolver", par.nlsolver);
    readParam<std::string>(argc, argv, "lsearch", par.lsearch);
    readParam<std::string>(argc, argv, "mask", par.mask_file);
    readParam<std::string>(argc, argv, "weights", par.weights_file);
    readParam<std::string>(argc, argv, "filter", par.filter_file);
    readParam<std::string>(argc, argv, "inverse_diagonal_hessian", par.inverse_diagonal_hessian_file);
    readParam<std::string>(argc, argv, "prior", par.prior_file);
    readParam<std::string>(argc, argv, "datapath", par.datapath);
    readParam<data_t>(argc, argv, "courant", par.courant);
    readParam<data_t>(argc, argv, "dt", par.dt);
    readParam<data_t>(argc, argv, "fmax", par.fmax);
    readParam<data_t>(argc, argv, "sx0", par.sx0);
    readParam<data_t>(argc, argv, "sy0", par.sy0);
    readParam<data_t>(argc, argv, "sz0", par.sz0);
    readParam<data_t>(argc, argv, "rx0", par.rx0);
    readParam<data_t>(argc, argv, "ry0", par.ry0);
    readParam<data_t>(argc, argv, "rz0", par.rz0);
    readParam<data_t>(argc, argv, "rdip0", par.rdip0);
    readParam<data_t>(argc, argv, "raz0", par.raz0);
    readParam<data_t>(argc, argv, "sxinc", par.sxinc);
    readParam<data_t>(argc, argv, "syinc", par.syinc);
    readParam<data_t>(argc, argv, "szinc", par.szinc);
    readParam<data_t>(argc, argv, "rxinc", par.rxinc);
    readParam<data_t>(argc, argv, "ryinc", par.ryinc);
    readParam<data_t>(argc, argv, "rzinc", par.rzinc);
    readParam<data_t>(argc, argv, "rdipinc", par.rdipinc);
    readParam<data_t>(argc, argv, "razinc", par.razinc);
    readParam<data_t>(argc, argv, "free_surface_stiffness", par.free_surface_stiffness);
    readParam<data_t>(argc, argv, "taper_strength", par.taper_strength);
    readParam<data_t>(argc, argv, "gl", par.gl);
    readParam<data_t>(argc, argv, "mxx", par.mxx);
    readParam<data_t>(argc, argv, "myy", par.myy);
    readParam<data_t>(argc, argv, "mzz", par.mzz);
    readParam<data_t>(argc, argv, "mxy", par.mxy);
    readParam<data_t>(argc, argv, "mxz", par.mxz);
    readParam<data_t>(argc, argv, "myz", par.myz);
    readParam<data_t>(argc, argv, "fangle", par.fangle);
    readParam<data_t>(argc, argv, "fazimuth", par.fazimuth);
    readParam<data_t>(argc, argv, "vpmin", par.vpmin);
    readParam<data_t>(argc, argv, "vpmax", par.vpmax);
    readParam<data_t>(argc, argv, "vsmin", par.vsmin);
    readParam<data_t>(argc, argv, "vsmax", par.vsmax);
    readParam<data_t>(argc, argv, "rhomin", par.rhomin);
    readParam<data_t>(argc, argv, "rhomax", par.rhomax);
    readParam<data_t>(argc, argv, "bs_controlx", par.bs_controlx);
    readParam<data_t>(argc, argv, "bs_controly", par.bs_controly);
    readParam<data_t>(argc, argv, "bs_controlz", par.bs_controlz);
    readParam<data_t>(argc, argv, "threshold", par.threshold);
    readParam<data_t>(argc, argv, "ls_a0", par.ls_a0);
    readParam<data_t>(argc, argv, "ls_a1", par.ls_a1);
    readParam<data_t>(argc, argv, "ls_c1", par.ls_c1);
    readParam<data_t>(argc, argv, "ls_c2", par.ls_c2);
    readParam<data_t>(argc, argv, "ls_max_step", par.ls_max_step);
    readParam<data_t>(argc, argv, "lambda", par.lambda);
    readParam<data_t>(argc, argv, "reg_xweight", par.reg_xweight);
    readParam<data_t>(argc, argv, "reg_yweight", par.reg_yweight);
    readParam<data_t>(argc, argv, "reg_zweight", par.reg_zweight);
    readParam<data_t>(argc, argv, "scale_source_log_clip", par.scale_source_log_clip);
    readParam<int>(argc, argv, "ns", par.ns);
    readParam<int>(argc, argv, "nr", par.nr);
    readParam<int>(argc, argv, "seismotype", par.seismotype);
    readParam<int>(argc, argv, "bc_top", par.bc_top);
    readParam<int>(argc, argv, "bc_bottom", par.bc_bottom);
    readParam<int>(argc, argv, "bc_left", par.bc_left);
    readParam<int>(argc, argv, "bc_right", par.bc_right);
    readParam<int>(argc, argv, "bc_front", par.bc_front);
    readParam<int>(argc, argv, "bc_back", par.bc_back);
    readParam<int>(argc, argv, "taper_top", par.taper_top);
    readParam<int>(argc, argv, "taper_bottom", par.taper_bottom);
    readParam<int>(argc, argv, "taper_left", par.taper_left);
    readParam<int>(argc, argv, "taper_right", par.taper_right);
    readParam<int>(argc, argv, "taper_front", par.taper_front);
    readParam<int>(argc, argv, "taper_back", par.taper_back);
    readParam<int>(argc, argv, "sinc_half_length", par.sinc_half_length);
    readParam<int>(argc, argv, "sub", par.sub);
    readParam<int>(argc, argv, "model_parameterization", par.model_parameterization);
    readParam<int>(argc, argv, "bs_nx", par.bs_nx);
    readParam<int>(argc, argv, "bs_ny", par.bs_ny);
    readParam<int>(argc, argv, "bs_nz", par.bs_nz);
    readParam<int>(argc, argv, "lbfgs_m", par.lbfgs_m);
    readParam<int>(argc, argv, "niter", par.niter);
    readParam<int>(argc, argv, "max_trial", par.max_trial);
    readParam<int>(argc, argv, "isave", par.isave);
    readParam<int>(argc, argv, "envelop", par.envelop);
    readParam<int>(argc, argv, "regularization", par.regularization);
    readParam<int>(argc, argv, "scale_source_times", par.scale_source_times);
    readParam<int>(argc, argv, "verbose", par.verbose);
    readParam<int>(argc, argv, "nthreads", par.nthreads);
    readParam<int>(argc, argv, "compression_rate", par.compression_rate);
    readParam<bool>(argc, argv, "mt", par.mt);
    readParam<bool>(argc, argv, "bsplines", par.bsplines);
    readParam<bool>(argc, argv, "soft_clip", par.soft_clip);
    readParam<bool>(argc, argv, "inversion1d", par.inversion1d);
    readParam<bool>(argc, argv, "normalize", par.normalize);
    readParam<bool>(argc, argv, "normalize_obj_func", par.normalize_obj_func);
    readParam<bool>(argc, argv, "ls_version", par.ls_version);
    readParam<bool>(argc, argv, "format", par.format);
    readParam<int>(argc, argv, "bs_mx", par.bs_mx);
    readParam<int>(argc, argv, "bs_my", par.bs_my);
    readParam<int>(argc, argv, "bs_mz", par.bs_mz);
}

void readParam(std::string parfile, param &par){
    
    std::ifstream ifs (parfile, std::ifstream::in);
    if (ifs.is_open())
    {
        std::string line;
        std::vector<std::string> lines = {"empty_line"};
        while ( getline (ifs, line) )
        {
            line = line + " ";
            if (line [0]!= '#' && line [0]!='//')
            {
                lines.push_back(line);
            }
        }
        int argc=lines.size();
        char * argv[argc];
        for (int i=0; i<argc; i++){
            argv[i] = (char *) (lines[i].c_str());
        }
        readParam(argc, argv, par);
        ifs.close();
    }
}

void readParameters(int argc, char **argv, param &par){
    std::string parfile="none";
    readParam<std::string>(argc, argv, "parfile", parfile);
    readParam(parfile, par);
    readParam(argc, argv, par);
    readCoord(par.srcoord, par);
    if (par.srcoord_from_file==true) {
        par.ns=par.sxyz.size();
    }
    else {
        par.sxyz.clear();
        par.rxyz.clear();
        for (int s=0; s<par.ns; s++){
            // compute shot x z coordinates
            std::vector<data_t> sxyz = {par.sx0+s*par.sxinc, par.sy0+s*par.syinc, par.sz0+s*par.szinc};

            // compute receiver x z dip coordinates
            std::vector<std::vector<data_t> > vec;
            for (int r=0; r<par.nr; r++){
                std::vector<data_t> rxyz = {par.rx0+r*par.rxinc, par.ry0+r*par.ryinc, par.rz0+r*par.rzinc, par.rdip0+r*par.rdipinc, par.raz0+r*par.razinc};
                vec.push_back(rxyz);
            }
            par.sxyz.push_back(sxyz);
            par.rxyz.push_back(vec);
        }
    }
}

void analyzeNLInversion(param &par)
{
    par.isave=std::max(1,par.isave);
    if (par.sub==0) par.sub=-1;
    if (par.verbose>0)
    {
        fprintf(stderr,"\n==========================\n Inversion parameters\n==========================\n");
        fprintf(stderr,"Non-linear solver = %s\n",par.nlsolver.c_str());
        fprintf(stderr,"Line search = %s\n",par.lsearch.c_str());
        fprintf(stderr,"Number of iterations = %d\n",par.niter);
        fprintf(stderr,"Maximum number of trials per iteration = %d\n",par.max_trial);
        fprintf(stderr,"Threshold to stop the inversion = %f\n",par.threshold);
    }
    std::string parameterization[4] = {"lambda, mu, rho, aniso","Vp, Vs, rho, aniso","IP, IS, rho, aniso","log(Vs/Vs0), log(Vp/Vs-sqrt(2)),log(rho/rho0),aniso"};
    if (par.nmodels>=3 && par.verbose>0) fprintf(stderr,"The model parameterization used in the elastic inversion is %s\n",parameterization[par.model_parameterization].c_str());
    if (par.soft_clip && par.verbose>0) fprintf(stderr,"Soft clipping is added to the inversion. It overrides the hard clipping\n");
    if (par.inversion1d && par.verbose>0) fprintf(stderr,"A 1D inversion will be performed in the z-dimension. If \"horizon\" is not provided, the extrapolation will be flat in the xy-dimensions where any possible B-splines will be ignored\n");
    if (par.mask_file != "none" && par.verbose>0) fprintf(stderr,"A gradient mask file is expected and will be applied at each trial\n");
    if (par.weights_file != "none" && par.verbose>0) fprintf(stderr,"A data weights file is expected and will be applied to modeled and observed data\n");
    if (par.filter_file != "none" && par.verbose>0) fprintf(stderr,"A 1D filter file is expected and will be applied to modeled and observed data\n");
    if (par.inverse_diagonal_hessian_file != "none" && par.nlsolver=="lbfgs" && par.verbose>0) fprintf(stderr,"An inverse diagonal Hessian is expected and will be used in the l-BFGS solver\n");
    if (par.prior_file != "none" && par.verbose>0) fprintf(stderr,"A prior model file is expected and will be used in the regularization if any\n");
    if (par.regularization>-1 && par.verbose>0) fprintf(stderr,"A Tikhonov regularization will be used, of order=%d and damping=%f\n",par.regularization,par.lambda);
    else if (par.verbose>0) fprintf(stderr,"No regularization will be used\n");
    if (par.scale_source_times>0 ) {
        par.normalize=0;
        par.envelop=0;
        if (par.verbose>0) fprintf(stderr,"The source time function will be rescaled for the first %d trials as part of Variable Projection method. All non-linear operators on the data will be deactivated\n",par.scale_source_times);
    }
    if (par.normalize && par.verbose>0) fprintf(stderr,"The modeled and observed data will be normalized trace by trace\n");
    if (par.envelop==1 && par.verbose>0) fprintf(stderr,"The envelop of modeled and observed data will be computed trace by trace\n");
    else if (par.envelop==2 && par.verbose>0) fprintf(stderr,"The envelop squared of modeled and observed data will be computed trace by trace\n");
    if (par.verbose>0) fprintf(stderr,"If \"ioutput\" is provided, the model, gradient, and residual will be saved every %d iterations\n",par.isave);
}

void analyzeBsplines(const hypercube<data_t> &domain, param &par)
{
    if (par.bsplines)
    {
        if (par.verbose>0) fprintf(stderr,"\n==========================\n Model parameterization with B-spline functions\n==========================\n");
        if (par.bs_mx.size()==0){
            par.bs_nx = std::max(3,par.bs_nx);
            if (par.verbose>0) fprintf(stderr,"The B-spline nodes are regularly spaced with %d nodes in x\n",par.bs_nx);

            axis<data_t> X = domain.getAxis(2);
            par.bs_controlx.resize(par.bs_nx); par.bs_mx.resize(par.bs_nx);
            par.bs_controlx[0] = X.o;
            data_t dx = X.d*(X.n-1)/(par.bs_nx-1);
            for (int i=1; i<par.bs_nx-1; i++) {
                par.bs_controlx[i] = par.bs_controlx[0] + i*dx;
                par.bs_mx[i] = 1;
            }
            par.bs_controlx[par.bs_nx-1] = X.o + (X.n-1)*X.d;
            par.bs_mx[0] = 2;
            par.bs_mx[par.bs_nx-1] = 2;
        }
        else{
            successCheck(par.bs_mx.size()==par.bs_controlx.size(),"Multiplicity and control vectors in x must be of the same size\n");
            par.bs_nx = par.bs_mx.size();
            successCheck(par.bs_nx>2,"There should be at least 3 B-spline nodes in x dimension\n");
            if (par.verbose>0) fprintf(stderr,"The B-spline nodes and their multiplicity in x are read from parameters list\n");
        }
        if (par.bs_my.size()==0){
            par.bs_ny = std::max(3,par.bs_ny);
            if (par.verbose>0) fprintf(stderr,"The B-spline nodes are regularly spaced with %d nodes in y\n",par.bs_ny);

            axis<data_t> Y = domain.getAxis(3);
            par.bs_controly.resize(par.bs_ny); par.bs_my.resize(par.bs_ny);
            par.bs_controly[0] = Y.o;
            data_t dy = Y.d*(Y.n-1)/(par.bs_ny-1);
            for (int i=1; i<par.bs_ny-1; i++) {
                par.bs_controly[i] = par.bs_controly[0] + i*dy;
                par.bs_my[i] = 1;
            }
            par.bs_controly[par.bs_ny-1] = Y.o + (Y.n-1)*Y.d;
            par.bs_my[0] = 2;
            par.bs_my[par.bs_ny-1] = 2;
        }
        else{
            successCheck(par.bs_my.size()==par.bs_controly.size(),"Multiplicity and control vectors in y must be of the same size\n");
            par.bs_ny = par.bs_my.size();
            successCheck(par.bs_ny>2,"There should be at least 3 B-spline nodes in y dimension\n");
            if (par.verbose>0) fprintf(stderr,"The B-spline nodes and their multiplicity in y are read from parameters list\n");
        }
        if (par.bs_mz.size()==0){
            par.bs_nz = std::max(3,par.bs_nz);
            if (par.verbose>0) fprintf(stderr,"The B-spline nodes are regularly spaced with %d nodes in z\n",par.bs_nz);

            axis<data_t> Z = domain.getAxis(1);
            par.bs_controlz.resize(par.bs_nz); par.bs_mz.resize(par.bs_nz);
            par.bs_controlz[0] = Z.o;
            data_t dz = Z.d*(Z.n-1)/(par.bs_nz-1);
            for (int i=1; i<par.bs_nz-1; i++) {
                par.bs_controlz[i] = par.bs_controlz[0] + i*dz;
                par.bs_mz[i] = 1;
            }
            par.bs_controlz[par.bs_nz-1] = Z.o + (Z.n-1)*Z.d;
            par.bs_mz[0] = 2;
            par.bs_mz[par.bs_nz-1] = 2;
        }
        else{
            successCheck(par.bs_mz.size()==par.bs_controlz.size(),"Multiplicity and control vectors in z must be of the same size\n");
            par.bs_nz = par.bs_mz.size();
            successCheck(par.bs_nz>2,"There should be at least 3 B-spline nodes in z dimension\n");
            if (par.verbose>0) fprintf(stderr,"The B-spline nodes and their multiplicity in z are read from parameters list\n");
        }
    }
}
