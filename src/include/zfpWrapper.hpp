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
    int nx=1, ny=1, nz=1, nc=1, nt=1, rate=0, ot=0, dt=1;

public:

    std::vector<std::shared_ptr<vecReg<data_t> > > _full_wfld;
#ifdef ENABLE_ZFP
    std::vector<std::vector<zfp::array3<data_t> > > _full_wfld_compressed;
#endif
    zfpWrapper(){}
    ~zfpWrapper(){}

    void initialize(const hypercube<data_t> & hyper, int compression_rate){

        successCheck(hyper.getNdim()==5,"The hypercube for initializing the zfpWrapper must have 5 dimensions, Z, X, Y, Component, Time\n");
        std::vector<axis<data_t> > axes = hyper.getAxes();
        nx = axes[1].n;
        nz = axes[0].n;
        ny = axes[2].n;
        nc = axes[3].n;
        nt = axes[4].n;
        ot = axes[4].o;
        dt = axes[4].d;
        rate = compression_rate;

#ifdef ENABLE_ZFP
        if (rate>0){
            _full_wfld_compressed = std::vector<std::vector<zfp::array3<data_t> > >(nt, std::vector<zfp::array3<data_t> >(nc, zfp::array3<data_t>(nz,nx,ny,rate)));
        }
        else{
            for (int it=0; it<nt; it++){
                std::shared_ptr<vecReg<data_t> > v = std::make_shared<vecReg<data_t> >(hypercube<data_t>(axes[0], axes[1], axes[2], axes[3]));
                v->zero();
                _full_wfld.push_back(v);
            }
        }
#else
        for (int it=0; it<nt; it++){
            std::shared_ptr<vecReg<data_t> > v = std::make_shared<vecReg<data_t> >(hypercube<data_t>(axes[0], axes[1], axes[2], axes[3]));
            v->zero();
            _full_wfld.push_back(v);
        }
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

        if (it>nt-1) return;

        int nxyz = nx*ny*nz;

#ifdef ENABLE_ZFP
        if (rate>0){
            for (int ic=0; ic<nc; ic++){
                _full_wfld_compressed[it][ic].set(v + ic*nxyz);
            }
        }
        else{
            memcpy(_full_wfld[it]->getVals(), v, nc*nxyz*sizeof(data_t));
        }
#else
        memcpy(_full_wfld[it]->getVals(), v, nc*nxyz*sizeof(data_t));
#endif
    }

    void get(data_t * &v, int it) const {

        if (it>nt-1) return;

        int nxyz = nx*ny*nz;

#ifdef ENABLE_ZFP
        if (rate>0){
            for (int ic=0; ic<nc; ic++){
                _full_wfld_compressed[it][ic].get(v + ic*nxyz);
            }
        }
        else{
            v = _full_wfld[it]->getVals();
        }
#else
        v = _full_wfld[it]->getVals();
#endif
    }

    void getFullWfld(std::shared_ptr<vecReg<data_t> > &full) const{

        std::vector<axis<data_t> > axes = _full_wfld[0]->getHyper()->getAxes();
        axis<data_t> T(nt,ot,dt);
        full = std::make_shared<vecReg<data_t> >(hypercube<data_t>(axes[0],axes[1],axes[2],axes[3],T));

        int ncxyz = nc*nx*ny*nz;
        data_t (* pfull) [ncxyz] = (data_t (*) [ncxyz]) full->getVals();
        for (int it=0; it<nt; it++) memcpy(pfull[it], _full_wfld[it]->getCVals(), ncxyz*sizeof(data_t));

    }

};