#pragma once

#include "vecReg.hpp"
#include "misc.hpp"
#ifdef ENABLE_ZFP
    #include "zfp.hpp"
    #include "array3.hpp"
#endif

struct zfpWrapper
{
private:
    int nx, ny, nz, nc, nt, rate;

public:

    std::shared_ptr<vecReg<data_t> > _full_wfld;
#ifdef ENABLE_ZFP
    std::vector<std::vector<zfp::array3<data_t> > > _full_wfld_compressed;
#endif
    zfpWrapper(){}
    ~zfpWrapper(){}

    void initialize(const hypercube<data_t> & hyper, int compression_rate){

        successCheck(hyper.getNdim()==5,"The initialization of the full wavefield must have 5 dimensions, Z, X, Y, Component, Time\n");
        std::vector<axis<data_t> > axes = hyper.getAxes();
        nx = axes[1].n;
        nz = axes[0].n;
        ny = axes[2].n;
        nc = axes[3].n;
        nt = axes[4].n;
        rate = compression_rate;

#ifdef ENABLE_ZFP
        if (rate>0){
            _full_wfld_compressed = std::vector<std::vector<zfp::array3<data_t> > >(nt, std::vector<zfp::array3<data_t> >(nc, zfp::array3<data_t>(nz,nx,ny,rate)));
        }
        else{
            if (axes[3].n > 1) _full_wfld=std::make_shared<vecReg<data_t> >(hyper);
            else _full_wfld=std::make_shared<vecReg<data_t> >(hypercube<data_t>(axes[0], axes[1], axes[2], axes[4]));
            _full_wfld->zero();
        }
#else
        if (axes[3].n > 1) _full_wfld=std::make_shared<vecReg<data_t> >(hyper);
        else _full_wfld=std::make_shared<vecReg<data_t> >(hypercube<data_t>(axes[0], axes[1], axes[2], axes[4]));
        _full_wfld->zero();
#endif

    }

    void allocate(data_t * &v) const {

        int nxyz = nx*ny*nz;

#ifdef ENABLE_ZFP
        if (rate>0){
            v = new data_t[nc*nxyz];
        }
#endif
    }

    void deallocate(data_t * v) const {

#ifdef ENABLE_ZFP
        if (rate>0){
            delete [] v;
        }
#endif
    }

    void set(const data_t * v, int it) {

        int nxyz = nx*ny*nz;

#ifdef ENABLE_ZFP
        if (rate>0){
            for (int ic=0; ic<nc; ic++){
                _full_wfld_compressed[it][ic].set(v + ic*nxyz);
            }
        }
        else{
            data_t (* pfull) [nc][nxyz] = (data_t (*) [nc][nxyz]) _full_wfld->getVals();
            memcpy(pfull[it], v, nc*nxyz*sizeof(data_t));
        }
#else
        data_t (* pfull) [nc][nxyz] = (data_t (*) [nc][nxyz]) _full_wfld->getVals();
        memcpy(pfull[it], v, nc*nxyz*sizeof(data_t));
#endif
    }

    void get(data_t * &v, int it) const {

        int nxyz = nx*ny*nz;

#ifdef ENABLE_ZFP
        if (rate>0){
            for (int ic=0; ic<nc; ic++){
                _full_wfld_compressed[it][ic].get(v + ic*nxyz);
            }
        }
        else{
            v = _full_wfld->getVals() + it*nc*nxyz;
        }
#else
        v = _full_wfld->getVals() + it*nc*nxyz;
#endif
    }

};