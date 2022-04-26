#pragma once

#include "vecReg.hpp"
#include "misc.hpp"

// operator to inject or extract sources/receivers from wavefield
class injector{
public:
    std::vector<data_t> _xw; // x ,y, z weights based on the type of injector
    std::vector<data_t> _yw;
    std::vector<data_t> _zw;
    std::vector<int> _xind; // array indices where injection/extraction occur
    std::vector<int> _yind;
    std::vector<int> _zind;
    int _npts;
    
    injector(){}
    virtual ~injector(){}

    // forward operator: inject one sample of the source into the wavefield (single component)
    virtual void inject(bool add, const data_t ** in, data_t * out, int nx, int ny, int nz, int nt, int npts, int it, int itr_min, int itr_max, const int * xind, const int * yind, const int * zind, const data_t * xw, const data_t * yw, const data_t * zw) const = 0;

    // adjoint operator: extract one sample from the wavefield into the receiver
    virtual void extract(bool add, const data_t * in, data_t ** out, int nx, int ny, int nz, int nt, int npts, int it, int itr_min, int itr_max, const int * xind, const int * yind, const int * zind, const data_t * xw, const data_t * yw, const data_t * zw) const = 0;

    //virtual void inject_gpu(bool add, const data_t ** in, data_t * out, int nx, int ny, int nz, int nt, int npts, int it, int itr_min, int itr_max, const int * xind, const int * yind, const int * zind, const data_t * xw, const data_t * yw, const data_t * zw) const = 0;
    //virtual void extract_gpu(bool add, const data_t * in, data_t ** out, int nx, int ny, int nz, int nt, int npts, int it, int itr_min, int itr_max, const int * xind, const int * yind, const int * zind, const data_t * xw, const data_t * yw, const data_t * zw) const = 0;

};

// injector/extractor by approximating the Dirac distribution using H norm operator and 3 moment conditions
// in 1D, 3 points are required to approximate delta: 
// di-1 = (alpha^2 - alpha)/(2*h*wi-1)
// di = (1 - alpha^2)/(h*wi)
// di+1 =  (alpha^2 + alpha)/(2*h*wi+1)
// alpha=(x-xi)/h where the singularity is @ x ; h is the sampling ; wi-1, wi, wi+1 are the coefficients from H (not including the sampling)
// in 3D, it is the product of the 1D approximations
//
// Near the boundaries, the pivot xi may be shifted one sample to make room for 3 pts ; alpha is computed accordingly
// Note that the adjoint is relative to the H norm, it is not the transpose of the forward
class delta_m3 : public injector{
public:
    delta_m3(){}
    virtual ~delta_m3(){};
    delta_m3(const hypercube<data_t> &range, const std::vector<std::vector<data_t> > &loc){
        int npts=loc.size();
        _npts = npts;
        _xind.resize(3*npts, 0);
        _yind.resize(3*npts, 0);
        _zind.resize(3*npts, 0);
        _xw.resize(6*npts, 0.0);
        _yw.resize(6*npts, 0.0);
        _zw.resize(6*npts, 0.0);
        findIndicesWeights(range, loc);
    }

    // find needed indices and weights on given grid
    void findIndicesWeights(const hypercube<data_t> &range, const std::vector<std::vector<data_t> > &loc);
    void inject(bool add, const data_t ** in, data_t * out, int nx, int ny, int nz, int nt, int npts, int it, int itr_min, int itr_max, const int * xind, const int * yind, const int * zind, const data_t * xw, const data_t * yw, const data_t * zw) const = 0;
    void extract(bool add, const data_t * in, data_t ** out, int nx, int ny, int nz, int nt, int npts, int it, int itr_min, int itr_max, const int * xind, const int * yind, const int * zind, const data_t * xw, const data_t * yw, const data_t * zw) const = 0;

    //void inject_gpu(bool add, const data_t ** in, data_t * out, int nx, int ny, int nz, int nt, int npts, int it, int itr_min, int itr_max, const int * xind, const int * yind, const int * zind, const data_t * xw, const data_t * yw, const data_t * zw) const = 0;
    //void extract_gpu(bool add, const data_t * in, data_t ** out, int nx, int ny, int nz, int nt, int npts, int it, int itr_min, int itr_max, const int * xind, const int * yind, const int * zind, const data_t * xw, const data_t * yw, const data_t * zw) const = 0;

};

// injector/extractor by approximating the Dirac distribution and its spatial derivative using H norm operator and 3 moment conditions
// this class is used to inject/extract moment tensor
// in 1D, 3 points are required to approximate the derivative of delta: 
// di-1= (1-2*alpha)/(2*h*wi-1) ; di=2*alpha/(h*wi) ; di+1=(-1-2*alpha)/(2*h*wi+1)
// alpha=(x-xi)/h where the singularity is @ x ; h is the sampling ; wi-1, wi, wi+1 & wi+2 are the coefficients from H (including the sampling)
// the first component of the source will have a delta derivative in x direction, the second one in y direction, the third in z direction
class ddelta_m3 : public injector{
public:
    ddelta_m3(){}
    virtual ~ddelta_m3(){};
    ddelta_m3(const hypercube<data_t> &range, const std::vector<std::vector<data_t> > &loc){
        int npts=loc.size();
        _npts = npts;
        _xind.resize(3*npts, 0);
        _zind.resize(3*npts, 0);
        _xw.resize(12*npts, 0.0);
        _zw.resize(12*npts, 0.0);
        findIndicesWeights(range, loc);
    }

    void findIndicesWeights(const hypercube<data_t> &range, const std::vector<std::vector<data_t> > &loc);
    void inject(bool add, const data_t ** in, data_t * out, int nx, int ny, int nz, int nt, int npts, int it, int itr_min, int itr_max, const int * xind, const int * yind, const int * zind, const data_t * xw, const data_t * yw, const data_t * zw) const = 0;
    void extract(bool add, const data_t * in, data_t ** out, int nx, int ny, int nz, int nt, int npts, int it, int itr_min, int itr_max, const int * xind, const int * yind, const int * zind, const data_t * xw, const data_t * yw, const data_t * zw) const = 0;

    //void inject_gpu(bool add, const data_t ** in, data_t * out, int nx, int ny, int nz, int nt, int npts, int it, int itr_min, int itr_max, const int * xind, const int * yind, const int * zind, const data_t * xw, const data_t * yw, const data_t * zw) const = 0;
    //void extract_gpu(bool add, const data_t * in, data_t ** out, int nx, int ny, int nz, int nt, int npts, int it, int itr_min, int itr_max, const int * xind, const int * yind, const int * zind, const data_t * xw, const data_t * yw, const data_t * zw) const = 0;

};

// injector/extractor for a linear DAS fiber along arbitrary direction defined by a dip and azimuth angles
// Approximating of double Dirac_x, Dirac_y, and Dirac_z using H norm operator and 3 moment conditions
// 3 points are required to approximate delta (see above)
// DAS measures int_{-gl/2}^{gl/2} d/dt(du(l,t)/dl)dl = d/dt(u(gl/2,t) - u(-gl/2,t)) (curvilinear)
// Note that the adjoint is relative to the H norm, it is not the transpose of the forward
// The class below performs only dipole injection/extraction from a scalar field, ex: ux(gl/2,t) - ux(-gl/2,t)
// Another function is needed to convert dipole data to/from strain data
class dipole_m3 : public injector{
public:
    dipole_m3(){}
    virtual ~dipole_m3(){};
    dipole_m3(const hypercube<data_t> &range, const std::vector<std::vector<data_t> > &loc, data_t gl){
        int npts=loc.size();
        _npts = npts;
        _xind.resize(6*npts, 0);
        _zind.resize(6*npts, 0);
        _xw.resize(12*npts, 0.0);
        _zw.resize(12*npts, 0.0);
        findIndicesWeights(range, loc, gl);
    }

    void findIndicesWeights(const hypercube<data_t> &range, const std::vector<std::vector<data_t> > &loc, data_t gl);
    void inject(bool add, const data_t ** in, data_t * out, int nx, int ny, int nz, int nt, int npts, int it, int itr_min, int itr_max, const int * xind, const int * yind, const int * zind, const data_t * xw, const data_t * yw, const data_t * zw) const = 0;
    void extract(bool add, const data_t * in, data_t ** out, int nx, int ny, int nz, int nt, int npts, int it, int itr_min, int itr_max, const int * xind, const int * yind, const int * zind, const data_t * xw, const data_t * yw, const data_t * zw) const = 0;

    //void inject_gpu(bool add, const data_t ** in, data_t * out, int nx, int ny, int nz, int nt, int npts, int it, int itr_min, int itr_max, const int * xind, const int * yind, const int * zind, const data_t * xw, const data_t * yw, const data_t * zw) const = 0;
    //void extract_gpu(bool add, const data_t * in, data_t ** out, int nx, int ny, int nz, int nt, int npts, int it, int itr_min, int itr_max, const int * xind, const int * yind, const int * zind, const data_t * xw, const data_t * yw, const data_t * zw) const = 0;

};