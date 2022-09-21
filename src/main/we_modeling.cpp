#include <string.h>
#include "we_op.hpp"
#include "IO.hpp"
#include "MpiWrapper.hpp"
#include "seplib.h"

typedef vecReg<data_t> vec;
typedef cvecReg<data_t> cvec;
typedef axis<data_t> ax;
typedef hypercube<data_t> hyper;

// Executable to model 3D seismic data and optionally save the full wavefield

int main(int argc, char **argv){

    int rank=0, size=1;
    MpiWrapper::init(&argc,&argv);
    MpiWrapper::setSizeRank(&size,&rank);
    fprintf (stderr,"\n====================\nSize of MPI communicator = %d ; current rank = %d\n====================\n",size,rank);

    initpar(argc,argv);

// Read parameters for wave propagation
    param par;
    readParameters(argc, argv, par);
    int verbose=par.verbose;
    if (rank>0) par.verbose=0;

// Set the maximum number of threads
    if (par.nthreads>0) omp_set_num_threads(par.nthreads);

// Read inputs/outputs files
    std::string source_file="none", model_file="none", wavefield_file="none", output_file="none";
    readParam<std::string>(argc, argv, "source", source_file);
    readParam<std::string>(argc, argv, "model", model_file);
    readParam<std::string>(argc, argv, "wavefield", wavefield_file);
    readParam<std::string>(argc, argv, "output", output_file);

    successCheck(source_file!="none","Source wavelet is not provided\n");
    successCheck(model_file!="none","Earth model is not provided\n");

    std::shared_ptr<vec> src = read<data_t>(source_file, par.format);
    std::shared_ptr<vec> model = read<data_t>(model_file, par.format);
    
// Analyze the input source time function and duplicate if necessary, analyze geometry
    analyzeGeometry(*model->getHyper(),par, par.verbose>0);
    std::shared_ptr<vec> allsrc = analyzeWavelet(src, par, par.verbose>0);
    //analyzeModel(*allsrc->getHyper(),model,par);
    //lop.dotProduct();
    
// If more than one shot is modeled, don't save the wavefield
    if (par.ns>1) par.sub=0;

// container for all the data
    std::shared_ptr<vec> allrcv;
    ax T = allsrc->getHyper()->getAxis(1);
    ax R(par.nr*par.nrcomp,0,1);
    if (rank==0) allrcv = std::make_shared<vec> (hyper(T,R));

// cumulative number of samples
    std::vector<int> ns_cumul(par.ns,0);
    ns_cumul[0]=0;
    if (par.ns>1) ns_cumul[1]=par.rxyz[0].size()*par.nrcomp*T.n;
    for (int s=2; s<par.ns; s++) ns_cumul[s] = ns_cumul[s-1] + par.rxyz[s-1].size()*par.nrcomp*T.n;

// loop over shots
    for (int s=rank; s<par.ns; s+=size){

        if (verbose>1) fprintf(stderr,"Start processing shot %d by process %d\n",s, rank);

// Build the appropriate wave equation operator
        nl_we_op_e * op;
        if (par.nmodels==2) op=new nl_we_op_a(*model->getHyper(),allsrc,s,par);
        if (par.nmodels==3) op=new nl_we_op_e(*model->getHyper(),allsrc,s,par);
        else if (par.nmodels==6) op=new nl_we_op_vti(*model->getHyper(),allsrc,s,par);

// Run the forward modeling
        std::shared_ptr<vec> rcv = std::make_shared<vec> (*op->getRange());
        rcv->zero();
        op->forward(false,model,rcv);

// Save one full wavefield if requested
        if ((rank==0) && (wavefield_file!="none") && (par.sub)>0) write<data_t>(op->_full_wfld, wavefield_file, par.format, par.datapath);
        delete op;

// Copy to the full container
        if (rank!=0) {
            successCheck( MpiWrapper::send(rcv->getCVals(), rcv->getN123(), 0) , "MPI error\n");
        }
        else{
            for (int task=1; task<size; task++){
                successCheck( MpiWrapper::recv(allrcv->getVals()+s*ns_cumul[s], rcv->getN123(), task) , "MPI error\n");
            }
            memcpy(allrcv->getVals()+s*ns_cumul[s], rcv->getCVals(),rcv->getN123()*sizeof(data_t));
        }
    }

    if ((rank==0) && (output_file!="none")) write<data_t>(allrcv, output_file, par.format, par.datapath);

    MpiWrapper::finalize();
    
return 0;
}