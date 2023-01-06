#pragma once

#include "operator.hpp"
#include "injector.hpp"
#include "param.hpp"
#include "zfpWrapper.hpp"

// Analyze model, sources,, and receiver geomtery, as well as boundary conditions
void analyzeGeometry(const hypercube<data_t> &model, param &par, bool verbose=true);

// check how many wavelets are provided and how many are needed ; duplicate when necessary
std::shared_ptr<vecReg<data_t> > analyzeWavelet(std::shared_ptr<vecReg<data_t> > src, const param &par, bool verbose=true);

// Analyze model and modify as necessary
void analyzeModel(const hypercube<data_t> &domain, std::shared_ptr<vecReg<data_t> > model, param &par);

// general non-linear wave equation operator: taking a model and computing data
class nl_we_op : virtual public nloper {
public:
    param _par;
    std::shared_ptr<vecReg<data_t> > _src;
    std::shared_ptr<vecReg<data_t> > _allsrc;
    // std::shared_ptr<vecReg<data_t> > _full_wfld;
    mutable zfpWrapper _zfp_wfld;

    nl_we_op(){}
    virtual ~nl_we_op(){}
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        
        successCheck( domain == _domain,"New domain must be the same as the original\n");
        successCheck( range.getAxis(1) == _range.getAxis(1),"New range must have the same time axis as the original\n");
        successCheck( range.getN123()/range.getAxis(1).n == _range.getN123()/_range.getAxis(1).n,  "New range must have the same number of traces as the original\n");

        _domain=domain;
        _range=range;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        if ((*mod->getHyper() != _domain) || ( dat->getHyper()->getAxis(1) != _range.getAxis(1)) || ( dat->getN123()/dat->getHyper()->getAxis(1).n != _range.getN123()/_range.getAxis(1).n) ) return false;
        else return true;
    }
    
    virtual void forward(bool add, const std::shared_ptr<vecReg<data_t> > mod, std::shared_ptr<vecReg<data_t> > dat) {
        successCheck(checkDomainRange(mod,dat),"Vectors hypercube do not match the operator domain and range\n");
        apply_forward(add, mod->getCVals(), dat->getVals());
    }
    virtual void jacobianT(bool add, std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > mod0, const std::shared_ptr<vecReg<data_t> > dat) {
        successCheck(checkDomainRange(mod,dat),"Vectors hypercube do not match the operator domain and range\n");
        successCheck(mod->getHyper()->isCompatible(*mod0->getHyper()),"Input and background models hypercubes are not compatible\n");
        apply_jacobianT(add, mod->getVals(), mod0->getCVals(), dat->getCVals());
    }
    virtual void apply_forward(bool add, const data_t * pmod, data_t * pdat) = 0;
    virtual void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat) = 0;
};

// non-linear isotropic elastic wave equation operator: taking an elastic model and computing data
class nl_we_op_e : virtual public nl_we_op {
public:
    nl_we_op_e(){}
    virtual ~nl_we_op_e(){}
    nl_we_op_e(const hypercube<data_t> &domain, const std::shared_ptr<vecReg<data_t> > allsrc, int s, param &par){
        _allsrc = allsrc;
        _src = std::make_shared<vecReg<data_t> >(hypercube<data_t>(allsrc->getHyper()->getAxis(1),allsrc->getHyper()->getAxis(2)));
        int n=_src->getN123();
        memcpy(_src->getVals(),allsrc->getVals()+s*n,n*sizeof(data_t));
        _domain = domain;
        _par = par;
        _par.sxyz={par.sxyz[s]};
        _par.rxyz={par.rxyz[s]};
        if (par.gl>0){
            _par.rdip.clear();
            _par.raz.clear();
            for (int r=0; r<par.rxyz[s].size(); r++){
                _par.rdip.push_back(par.rxyz[s][r][3]);
                _par.raz.push_back(par.rxyz[s][r][4]);
            }
        }
        axis<data_t> X(_par.rxyz[0].size(),0,1);
        axis<data_t> C(_par.nrcomp,0,1);
        _range = hypercube<data_t>(_src->getHyper()->getAxis(1),X,C);
    }
    nl_we_op_e * clone() const {
        param par = _par;
        nl_we_op_e * op = new nl_we_op_e(_domain,_src,0,par);
        return op;
    }
    
    // convert Vp, Vs, rho to lambda, mu, rho and vice versa
    // mu = rho.vs2
    // lambda = rho.(vp2 - 2.vs2)
    virtual void convert_model(data_t * m, int n, bool forward) const;

    // refer to the SBP notes for gradients expression
    virtual void compute_gradients(const data_t * model, const data_t ** u_full, const data_t * curr, const data_t * u_x,  const data_t * u_y, const data_t * u_z, data_t * tmp, data_t * grad, const param &par, int nx, int ny, int nz, int it, data_t dx, data_t dy, data_t dz, data_t dt) const;
    virtual void propagate(bool adj, const data_t * model, const data_t * src, data_t * rcv, const injector * inj, const injector * ext, data_t * grad, const param &par, int nx, int ny, int nz, data_t dx, data_t dy, data_t dz, data_t ox, data_t oy, data_t oz) const;

    //virtual void compute_gradients_gpu(const data_t * model, const data_t * u_full, const data_t * curr, const data_t * u_x,  const data_t * u_y, const data_t * u_z, data_t * tmp, data_t * grad, const param &par, int nx, int ny, int nz, int it, data_t dx, data_t dy, data_t dz, data_t dt) const;
    //virtual void propagate_gpu(bool adj, const data_t * model, const data_t * src, data_t * rcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int ny, int nz, data_t dx, data_t dy, data_t dz, data_t ox, data_t oy, data_t oz) const;

    virtual void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    virtual void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat);
};

// non-linear VTI elastic wave equation operator: taking an elastic model and computing data
class nl_we_op_vti : virtual public nl_we_op_e {
    
public:
    using nl_we_op_e::nl_we_op_e;
    ~nl_we_op_vti(){}
    nl_we_op_vti * clone() const {
        param par = _par;
        nl_we_op_vti * op = new nl_we_op_vti(_domain,_src,0,par);
        return op;
    }
    
    // convert Vp, Vs, rho, delta, epsilon, gamma to generalized lambda, generalized mu, rho, c13, eps, c66 and vice versa
    // c33 = lambda+2.mu
    // c11 = (1+2.eps).c33
    // c55 = mu
    // c13 = sqrt[2.c33.(c33-c55).del + (c33-c55)^2] - c55
    // c66 = (1+2.gam).mu
    // c12 = c11 - 2.c66
    // mu = rho.vs2
    // lambda = rho.(vp2 - 2.vs2)
    void convert_model(data_t * m, int n, bool forward) const;

    // refer to the SBP notes for gradients expression
    void compute_gradients(const data_t * model, const data_t ** u_full, const data_t * curr, const data_t * u_x,  const data_t * u_y, const data_t * u_z, data_t * tmp, data_t * grad, const param &par, int nx, int ny, int nz, int it, data_t dx, data_t dy, data_t dz, data_t dt) const;
    void propagate(bool adj, const data_t * model, const data_t * src, data_t * rcv, const injector * inj, const injector * ext, data_t * grad, const param &par, int nx, int ny, int nz, data_t dx, data_t dy, data_t dz, data_t ox, data_t oy, data_t oz) const;

    //void compute_gradients_gpu(const data_t * model, const data_t * u_full, const data_t * curr, const data_t * u_x,  const data_t * u_y, const data_t * u_z, data_t * tmp, data_t * grad, const param &par, int nx, int ny, int nz, int it, data_t dx, data_t dy, data_t dz, data_t dt) const;
    //void propagate_gpu(bool adj, const data_t * model, const data_t * src, data_t * rcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int ny, int nz, data_t dx, data_t dy, data_t dz, data_t ox, data_t oy, data_t oz) const;
};

// non-linear acoustic wave equation operator: taking an acoustic model and computing data
class nl_we_op_a : virtual public nl_we_op_e {
public:
    std::shared_ptr<vecReg<data_t> > _transmission; // vector containing transmission coefficients (in [0,1]) of the top boundary
    nl_we_op_a(){}
    virtual ~nl_we_op_a(){}
    nl_we_op_a(const hypercube<data_t> &domain, const std::shared_ptr<vecReg<data_t> > allsrc, int s, param &par){
        _allsrc = allsrc;
        _src = std::make_shared<vecReg<data_t> >(hypercube<data_t>(allsrc->getHyper()->getAxis(1)));
        int n=_src->getN123();
        memcpy(_src->getVals(),allsrc->getVals()+s*n,n*sizeof(data_t));
        _domain = domain;
        _par = par;
        _par.sxyz={par.sxyz[s]};
        _par.rxyz={par.rxyz[s]};
        axis<data_t> X(_par.rxyz[0].size(),0,1);
        axis<data_t> C(1,0,1);
        _range = hypercube<data_t>(_src->getHyper()->getAxis(1),X,C);
    }
    nl_we_op_a * clone() const {
        param par = _par;
        nl_we_op_a * op = new nl_we_op_a(_domain,_src,0,par);
        return op;
    }
    
    // convert Vp, rho to bulk modulus K, 1/rho and vice versa
    void convert_model(data_t * m, int n, bool forward) const;

    // refer to the SBP notes for gradients expression
    void compute_gradients(const data_t * model, const data_t ** u_full, const data_t * curr, data_t * tmp, data_t * grad, const param &par, int nx, int ny, int nz, int it, data_t dx, data_t dy, data_t dz, data_t dt) const;
    void propagate(bool adj, const data_t * model, const data_t * src, data_t * rcv, const injector * inj, const injector * ext, data_t * grad, const param &par, int nx, int ny, int nz, data_t dx, data_t dy, data_t dz, data_t ox, data_t oy, data_t oz) const;

    //void compute_gradients_gpu(const data_t * model, const data_t * u_full, const data_t * curr, const data_t * u_x, const data_t * u_z, data_t * tmp, data_t * grad, const param &par, int nx, int nz, int it, data_t dx, data_t dz, data_t dt) const;
    //void propagate_gpu(bool adj, const data_t * model, const data_t * allsrc, data_t * allrcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int nz, data_t dx, data_t dz, data_t ox, data_t oz, int s) const;
};