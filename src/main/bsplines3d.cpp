#include <string.h>
#include <unistd.h>
#include "operator.hpp"
#include "lsolver.hpp"
#include "bsplines.hpp"
#include "param.hpp"
#include "IO.hpp"
#include "seplib.h"

typedef vecReg<data_t> vec;
typedef cvecReg<data_t> cvec;
typedef axis<data_t> ax;
typedef hypercube<data_t> hyper;


void printdoc(){
    std::string doc = "\nDescription:\n"
    "   Apply cubic B-splines smoothing.\n"
    "\nInput/Output:\n"
    "   Provide input as 'input=file.H' or '< file.H' and output as 'output=file.H' or '> file.H'. '<' and '>' are valid for SEPlib format only.\n"
    "   bsmodel - string - ['none'] : optional B-splines model to map back to the input space. Must be compatible with the input and parameters below.\n"
    "\nParameters:\n"
    "   nx,ny,nz - int - [3] :\n\t\tnumber of control points in each direction. If provided, will override the parameters below.\n"
    "   controlx,controly,controlz - [float] :\n\t\tarrays of control points manually entered. Must contain the first and last physical points. To be used in conjunction with mx,my,mz.\n"
    "   mx,my,mz - [int] :\n\t\tarrays of control points multiplicity manually entered.\n"
    "   niter - int - [0]:\n\t\tnumber of iterations to build a least-squares solution for the B-splines model (less smoothing is expected).\n"
    "   threshold - float - [0]:\n\t\tconvergence threshold to stop the solver.\n"
    "   bsoutput - string - ['none']:\n\t\tname of output file to store the B-splines model.\n"
    "   format - bool - [0]:\n\t\tdata format for IO. 0 for SEPlib, 1 for binary with description file.\n"
    "   datapath - string - ['none']:\n\t\tpath for output binaries.\n"
    "\nNote:\n"
    "   The number of control points must be >= 3 in all cases.\n"
    "\nExamples:\n"
    "   BSPLINES3D.x < infile.H nx=11 ny=15 nz=23 > oufile.H.\n"
    "   BSPLINES3D.x < infile.H controlx=0,1,5.5,10 controly=0,3,7,10 controlz=0,5,20 mx=2,1,1,2 my=2,1,1,2 mz=2,1,2 > oufile.H\n"
    "   BSPLINES3D.x < infile.H bsmodel=bsm.H nx=11 ny=15 nz=23 > oufile.H\n"
    "\n";
    fprintf(stderr,doc.c_str());
}

int main(int argc, char **argv){

    if (argc == 1 && isatty(STDIN_FILENO)==true) {printdoc(); return 0;}

	initpar(argc,argv);
    omp_set_num_threads(1);

    std::string input_file="in", bsmodel_file="none", output_file="out", bsoutput_file="none", datapath="none";
    int nx=3, ny=3, nz=3, niter=0;
    data_t threshold = 0;
    bool format=0;
    std::vector<data_t> controlx={0}, controly={0}, controlz={0}; // locations of the control points that define the B-splines
    std::vector<int> mx={0}, my={0}, mz={0}; // multiplicity of the control points that define the B-splines

	readParam<std::string>(argc, argv, "input", input_file);
	readParam<std::string>(argc, argv, "bsmodel", bsmodel_file);
	readParam<std::string>(argc, argv, "output", output_file);
	readParam<std::string>(argc, argv, "bsoutput", bsoutput_file);
	readParam<std::string>(argc, argv, "datapath", datapath);
    readParam<int>(argc, argv, "nx", nx);
    readParam<int>(argc, argv, "ny", ny);
    readParam<int>(argc, argv, "nz", nz);
    readParam<int>(argc, argv, "niter", niter);
    readParam<int>(argc, argv, "mx", mx);
    readParam<int>(argc, argv, "my", my);
    readParam<int>(argc, argv, "mz", mz);
    readParam<bool>(argc, argv, "format", format);
    readParam<data_t>(argc, argv, "controlx", controlx);
    readParam<data_t>(argc, argv, "controly", controly);
    readParam<data_t>(argc, argv, "controlz", controlz);
    readParam<data_t>(argc, argv, "threshold", threshold);
    
    std::shared_ptr<vec> input = read<data_t>(input_file,format);
    std::shared_ptr<vec> output = std::make_shared<vec>(*input->getHyper());
    const data_t* pinput = input->getCVals();
    data_t* poutput = output->getVals();

    // Build the knots, multiplicity, and the B-splines smoothing operator
// ----------------------------------------------------------------------------------------//
    if (nx > 2){
        ax X = input->getHyper()->getAxis(2);
        controlx.resize(nx); mx.resize(nx);
        controlx[0] = X.o;
        data_t dx = X.d*(X.n-1)/(nx-1);
        for (int i=1; i<nx-1; i++) {
            controlx[i] = controlx[0] + i*dx;
            mx[i] = 1;
        }
        controlx[nx-1] = X.o + (X.n-1)*X.d;
        mx[0] = 2;
        mx[nx-1] = 2;
    }
    else successCheck(mx.size() == controlx.size(),"Multiplicity and control vectors must be of the same size\n");
    if (ny > 2){
        ax Y = input->getHyper()->getAxis(3);
        controly.resize(ny); my.resize(ny);
        controly[0] = Y.o;
        data_t dy = Y.d*(Y.n-1)/(ny-1);
        for (int i=1; i<ny-1; i++) {
            controly[i] = controly[0] + i*dy;
            my[i] = 1;
        }
        controly[ny-1] = Y.o + (Y.n-1)*Y.d;
        my[0] = 2;
        my[ny-1] = 2;
    }
    else successCheck(my.size() == controly.size(),"Multiplicity and control vectors must be of the same size\n");
    
    if (nz > 2){
        ax Z = input->getHyper()->getAxis(1);
        controlz.resize(nz); mz.resize(nz);
        controlz[0] = Z.o;
        data_t dz = Z.d*(Z.n-1)/(nz-1);
        for (int i=1; i<nz-1; i++) {
            controlz[i] = controlz[0] + i*dz;
            mz[i] = 1;
        }
        controlz[nz-1] = Z.o + (Z.n-1)*Z.d;
        mz[0] = 2;
        mz[nz-1] = 2;
    }
    else successCheck(mz.size() == controlz.size(),"Multiplicity and control vectors must be of the same size\n");
    
    std::vector<data_t> kx;
    std::vector<data_t> ky;
    std::vector<data_t> kz;
    setKnot(kx,controlx,mx);
    setKnot(ky,controly,my);
    setKnot(kz,controlz,mz);

    std::vector<ax> axes = input->getHyper()->getAxes();
    axes[0].n=controlz.size(); 
    axes[1].n=controlx.size();
    axes[2].n=controly.size();
    
    std::shared_ptr<vec> bsmodel;
    if (bsmodel_file=="none"){
        bsmodel = std::make_shared<vec> (hyper(axes));
        fillin3d(bsmodel,input,controlx,controly,controlz);
    }
    else{
        bsmodel = read<data_t>(bsmodel_file,format);
        successCheck(bsmodel->getHyper()->isCompatible(hyper(axes)),"The B-splines model is not compatible with the input and B-splines parameters\n");
    }
  
    duplicate3d D(*bsmodel->getHyper(),mx,my,mz);
    bsplines3d B(*D.getRange(),*input->getHyper(),kx,ky,kz);
    chainLOper * BD = new chainLOper(&B,&D);

    if (niter>0)
    {
        // set the least-squares problem
        llsq prob(BD, bsmodel, input);

        // Set the linear solver
        lsolver * sol = new cgls(niter,threshold);

        // Invert the B-splines model
        sol->run(&prob,true);

        delete sol;
    }

    BD->forward(false,bsmodel,output);
    delete BD;

// ----------------------------------------------------------------------------------------//

    if (output_file!="none") write<data_t>(output,output_file,format,datapath);
    if (bsoutput_file!="none") write<data_t>(bsmodel,bsoutput_file,format,datapath);

    return 0;
}