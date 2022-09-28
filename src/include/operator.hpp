#pragma once

#include "vecReg.hpp"
#include "misc.hpp"

// generic class for non-linear operators
class nloper {
protected:
    hypercube<data_t> _domain;
    hypercube<data_t> _range;
public:
    nloper(){}
    virtual ~nloper(){}
    nloper(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        _domain = domain;
        _range = range;
    }
    const hypercube<data_t> * getDomain() const {return &_domain;}
    const hypercube<data_t> * getRange() const {return &_range;}
    virtual void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        _domain = domain;
        _range = range;
    }
    virtual bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {return checkN123(mod, dat);}
    bool checkSame(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        if ((*mod->getHyper() != _domain) || (*dat->getHyper() != _range)) return false;
        else return true;
    }
    bool checkCompatible(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        if ((mod->getHyper()->isCompatible(_domain)) && (dat->getHyper()->isCompatible(_range))) return true;
        else return false;
    }
    bool checkN123(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        if ((mod->getN123() != _domain.getN123()) || (dat->getN123() != _range.getN123())) return false;
        else return true;
    }
    // clone the object
    virtual nloper * clone() const = 0;

    // d = F(m)
    virtual void forward(bool add, const std::shared_ptr<vecReg<data_t> > mod, std::shared_ptr<vecReg<data_t> > dat) {
        successCheck(checkDomainRange(mod,dat),"Vectors hypercube do not match the operator domain and range\n");
        apply_forward(add, mod->getCVals(), dat->getVals());
    }
    // m = invF(d) if relevant
    virtual void inverse(bool add, std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) {
        successCheck(checkDomainRange(mod,dat),"Vectors hypercube do not match the operator domain and range\n");
        apply_inverse(add, mod->getVals(), dat->getCVals());
    }
    // d = (dF/dm0)*m
    virtual void jacobian(bool add, const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > mod0, const std::shared_ptr<vecReg<data_t> > dat) {
        successCheck(checkDomainRange(mod,dat),"Vectors hypercube do not match the operator domain and range\n");
        successCheck(mod->getHyper()->isCompatible(*mod0->getHyper()),"Input and background models hypercubes are not compatible\n");
        apply_jacobian(add, mod->getCVals(), mod0->getCVals(), dat->getVals());
    }
    //  m = (dF/dm0)' * d
    virtual void jacobianT(bool add, std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > mod0, const std::shared_ptr<vecReg<data_t> > dat) {
        successCheck(checkDomainRange(mod,dat),"Vectors hypercube do not match the operator domain and range\n");
        successCheck(mod->getHyper()->isCompatible(*mod0->getHyper()),"Input and background models hypercubes are not compatible\n");
        apply_jacobianT(add, mod->getVals(), mod0->getCVals(), dat->getCVals());
    }
    virtual void apply_forward(bool add, const data_t * pmod, data_t * pdat) {
        successCheck(false,"This function from the abstract non-linear operator must not be called\n");
    }
    virtual void apply_inverse(bool add, data_t * mod, const data_t * dat) {
        successCheck(false,"This function from the abstract non-linear operator must not be called\n");
    }
    virtual void apply_jacobian(bool add, const data_t * pmod, const data_t * pmod0, data_t * pdat) {
        successCheck(false,"This function from the abstract non-linear operator must not be called\n");
    }
    virtual void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat) {
        successCheck(false,"This function from the abstract non-linear operator must not be called\n");
    }

    // Approximate dot product test
    // Forward linear operator df/dm
    // Adjoint linear operator (df/dm)' implicit in the jacobian (df/dm)'.d
    // (F(m + eps.dm) - F(m))/eps ~= (dF/dm).dm     for eps << 1
    // Perform dot product < (dF/dm).dm , d > ~= < dm , (dF/dm)'.d >
    void dotProduct();
};

// generic class for linear operators
class loper : virtual public nloper {
protected:

public:
    loper(){}
    virtual ~loper(){}
    virtual void adjoint(bool add, std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) {
        successCheck(checkDomainRange(mod,dat),"Vectors hypercube do not match the operator domain and range\n");
        apply_adjoint(add, mod->getVals(), dat->getCVals());
    }
    virtual loper * clone() const = 0;
    virtual void apply_adjoint(bool add, data_t * pmod, const data_t * pdat) {
        successCheck(false,"This function from the abstract linear operator must not be called.\n");
    }
    void apply_jacobian(bool add, const data_t * pmod, const data_t * pmod0, data_t * pdat) {
        apply_forward(add, pmod, pdat);
    }
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat) {
        apply_adjoint(add, pmod, pdat);
    }
    // dot product test
    void dotProduct();
};

// identity operator
class identity : public loper {
protected:

public:
    identity(){}
    ~identity(){}
    identity(const hypercube<data_t> &domain){
        _domain = domain;
        _range = domain;
    }
    identity * clone() const {
        identity * op = new identity(_domain);
        return op;
    }
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        successCheck(domain.getN123()==_domain.getN123(),"The new domain must have the same number of samples\n");
        successCheck(range.getN123()==_range.getN123(),"The new range must have the same number of samples\n");
        successCheck(domain.getN123()==range.getN123(),"The new domain and range must have the same number of samples\n");
        _domain = domain;
        _range = range;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        bool ans1 = (mod->getN123()==_domain.getN123());
        bool ans2 = (dat->getN123()==_range.getN123());
        bool ans3 = (dat->getN123()==mod->getN123());
        return (ans1 && ans2 && ans3);
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat) {
        for (long i=0; i<_domain.getN123(); i++) pdat[i] = add*pdat[i] + pmod[i];
    }
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat) {
        for (long i=0; i<_domain.getN123(); i++) pmod[i] = add*pmod[i] + pdat[i];
    }
};

// Chain of two non-linear (or linear) operators L(R(m))
class chainNLOper : public nloper {
protected:
    nloper * _L; // left operator
    nloper * _R; // right operator
    std::shared_ptr<vecReg<data_t> > _m0; // intermediate vectors used for chained forward and jacobian
    std::shared_ptr<vecReg<data_t> > _v;

public:
    chainNLOper(){}
    ~chainNLOper(){delete _L; delete _R;}
    chainNLOper(nloper * L, nloper * R){
        successCheck(R->getRange()->getN123()==L->getDomain()->getN123(),"Number of samples in range of R must be the same as the one in domain of L\n");
        _domain = *R->getDomain();
        _range = *L->getRange();
        _L=L->clone();
        _R=R->clone();
        _m0 = std::make_shared<vecReg<data_t> >(*_R->getRange());
        _v = std::make_shared<vecReg<data_t> >(*_R->getRange());
        _m0->zero();
        _v->zero();
    }
    chainNLOper * clone() const {
        chainNLOper * op = new chainNLOper(_L, _R);
        return op;
    }
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        successCheck(domain.getN123()==_domain.getN123(),"The new domain must have the same number of samples\n");
        successCheck(range.getN123()==_range.getN123(),"The new range must have the same number of samples\n");
        _domain = domain;
        _range = range;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        bool ans1 = (mod->getN123()==_domain.getN123());
        bool ans2 = (dat->getN123()==_range.getN123());
        return (ans1 && ans2);
    }    
    void apply_forward(bool add, const data_t * pmod, data_t * pdat) {
        _R->apply_forward(false, pmod, _v->getVals());
        _L->apply_forward(add, _v->getCVals(), pdat);
    }
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat) {
        _R->apply_forward(false, pmod0, _m0->getVals());
        _L->apply_jacobianT(false, _v->getVals(), _m0->getCVals(), pdat);
        _R->apply_jacobianT(add, pmod, pmod0, _v->getCVals());
    }
};

// Chain of two linear operators L.R.m
class chainLOper : public loper {
protected:
    loper * _L; // left operator
    loper * _R; // right operator
    std::shared_ptr<vecReg<data_t> > _v;
public:
    chainLOper(){}
    ~chainLOper(){delete _L; delete _R;}
    chainLOper(loper * L, loper * R){
        successCheck(R->getRange()->getN123()==L->getDomain()->getN123(),"Number of samples in range of R must be the same as the one in domain of L\n");
        _domain = *R->getDomain();
        _range = *L->getRange();
        _L=L->clone();
        _R=R->clone();
        _v = std::make_shared<vecReg<data_t> >(*_R->getRange());
        _v->zero();
    }
    chainLOper * clone() const {
        chainLOper * op = new chainLOper(_L, _R);
        return op;
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat) {
        _R->apply_forward(false, pmod, _v->getVals());
        _L->apply_forward(add, _v->getCVals(), pdat);
    }
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat) {
        _L->apply_adjoint(false, _v->getVals(), pdat);
        _R->apply_adjoint(add, pmod, _v->getCVals());
    }
};





// ----------------------------------------------------------------------------------------//
// non-linear operators
// ----------------------------------------------------------------------------------------//

// non-linear soft clip operator h(x) = ig(f(g(x)))
// g(x) = (2.(x-xmean)/Dx)^p    Dx = xmax-xmin  xmean = (xmax+xmin)/2
// f(g) = g-g^q/q if |g|<1 and = sign(g).(q-1)/q otherwise
// ig(f) = Dx/2.f^(1/p) + xmean
// h(x) = Dx/2.(g-g^q/q)^(1/p) + xmean if |g|<1  and = Dx/2.(sign(g).(q-1)/q)^(1/p) + xmean otherwise 
// h'(x) = (2.(x-xmean)/Dx)^(p-1).(1-g^(q-1)).(g-g^q/q)^(1/p-1) if |g|<1 and = 0 otherwise
// p and q are positive odd integers
// xmax and xmin are actually re-scaled so that the final clipping is within the provided range
class softClip : public nloper {
protected:
    data_t _xmin;
    data_t _xmax;
    int _p;
    int _q;
public:
    softClip(){}
    ~softClip(){}
    softClip(const hypercube<data_t> &domain, data_t xmin, data_t xmax, int p=1, int q=9){
        successCheck(xmax>xmin,"The upper bound must be strictly greater than the lower bound\n");
        successCheck((p%2==1) && (q%2==1),"The powers p and q must be odd numbers\n");
        _domain = domain;
        _range = domain;
        _xmin = xmin;
        _xmax = xmax;
        _p = p;
        _q = q;
    }
    softClip * clone() const {
        softClip * op = new softClip(_domain, _xmin, _xmax, _p, _q);
        return op;
    }   
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat);
};

// soft clip of the elastic model and the Vs/Vp ratio
class emodelSoftClip : public nloper {
    data_t _vpmin;
    data_t _vpmax;
    data_t _vsmin;
    data_t _vsmax;
    data_t _rhomin;
    data_t _rhomax;
    data_t _spratio;
    int _p;
    int _q;
public:
    emodelSoftClip(){}
    ~emodelSoftClip(){}
    emodelSoftClip(const hypercube<data_t> &domain, data_t vpmin, data_t vpmax, data_t vsmin, data_t vsmax, data_t rhomin, data_t rhomax, data_t spratio=1/sqrt(2.00001), int p=9, int q=9){
        successCheck((domain.getNdim()==4) && (domain.getAxis(4).n>=3),"The domain must be 4D with 4th dimension containing at least 3 fields\n");
        _domain = domain;
        _range = domain;
        _vpmin=vpmin; _vpmax=vpmax; _vsmin=vsmin; _vsmax=vsmax; _rhomin=rhomin; _rhomax=rhomax; _spratio=spratio;
        _p = p;
        _q = q;
    }
    emodelSoftClip * clone() const {
        emodelSoftClip * op = new emodelSoftClip(_domain, _vpmin, _vpmax, _vsmin, _vsmax, _rhomin, _rhomax, _spratio, _p, _q);
        return op;
    }   
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat);
};

// Elastic model parameterization using lambda, mu, density, then anisotropy if applicable
class lam_mu_rho : public nloper {
public:
    lam_mu_rho(){}
    ~lam_mu_rho(){}
    lam_mu_rho(const hypercube<data_t> &domain){
        successCheck((domain.getNdim()>=2) && (domain.getAxis(domain.getNdim()).n>=3),"The domain must be at least 2D with the last dimension containing at least 3 fields\n");
        _domain = domain;
        _range = domain;
    }
    lam_mu_rho * clone() const {
        lam_mu_rho * op = new lam_mu_rho(_domain);
        return op;
    }   
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat);
    void apply_inverse(bool add, data_t * pmod, const data_t * pdat);
};

// Elastic model parameterization using P-impedance, S-impedance, density, then anisotropy if applicable
class ip_is_rho : public nloper {
public:
    ip_is_rho(){}
    ~ip_is_rho(){}
    ip_is_rho(const hypercube<data_t> &domain){
        successCheck((domain.getNdim()>=2) && (domain.getAxis(domain.getNdim()).n>=3),"The domain must be at least 2D with the last dimension containing at least 3 fields\n");
        _domain = domain;
        _range = domain;
    }
    ip_is_rho * clone() const {
        ip_is_rho * op = new ip_is_rho(_domain);
        return op;
    }   
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat);
    void apply_inverse(bool add, data_t * pmod, const data_t * pdat);
};

// Elastic model parameterization using log(Vs/Vs0), log(Vp/Vs - sqrt(2)), log(rho/rho0), then anisotropy if applicable
// Vs0 and rho0 are some fixed reference values
class vs_vpvs_rho : public nloper {
protected:
    data_t _vs0;
    data_t _rho0;
public:
    vs_vpvs_rho(){}
    ~vs_vpvs_rho(){}
    vs_vpvs_rho(const hypercube<data_t> &domain, data_t vs0=1.0, data_t rho0=1.0){
        successCheck((domain.getNdim()>=2) && (domain.getAxis(domain.getNdim()).n>=3),"The domain must be at least 2D with the last dimension containing at least 3 fields\n");
        _domain = domain;
        _range = domain;
        _vs0=vs0;
        _rho0=rho0;
    }
    vs_vpvs_rho * clone() const {
        vs_vpvs_rho * op = new vs_vpvs_rho(_domain,_vs0,_rho0);
        return op;
    }   
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat);
    void apply_inverse(bool add, data_t * pmod, const data_t * pdat);
};


// ----------------------------------------------------------------------------------------//
// linear operators
// ----------------------------------------------------------------------------------------//

// operator to resampler data in time (along the fast axis)
class resampler : public loper {
public:
    resampler(){}
    virtual ~resampler(){}
    resampler(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        setDomainRange(domain,range);
    }
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        
        successCheck( domain.getAxis(1).o == range.getAxis(1).o,"Domain and range must have the same time origin\n");
        successCheck( domain.getN123()/domain.getAxis(1).n == range.getN123()/range.getAxis(1).n,"Domain and range must have the same number of traces\n");

        _domain=domain;
        _range=range;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        if ((mod->getHyper()->getAxis(1).o != dat->getHyper()->getAxis(1).o) || (mod->getN123()/mod->getHyper()->getAxis(1).n != dat->getN123()/dat->getHyper()->getAxis(1).n)) return false;
        else return true;
    }
}; 

// linear resampler
class linear_resampler : public resampler{
public:
    // inherit constructor from the base class
    using resampler::resampler;
    ~linear_resampler(){}
    linear_resampler * clone() const{
        linear_resampler * resampler = new linear_resampler(_domain, _range);
        return resampler;
    }

    // forward operator (often interpolation)
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    // adjoint operator (often decimation)
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);

};

// sinc resampler
class sinc_resampler : public resampler{
protected:
    int _hl;
    data_t _alpha;

public:
    sinc_resampler(){}
    ~sinc_resampler(){}
    sinc_resampler(const hypercube<data_t> &domain, const hypercube<data_t> &range, int half_length, data_t alpha=0.5):resampler(domain, range){
   
        successCheck(half_length>0,"The half length of the sinc filter must be larger than 0\n");
        _hl = half_length;
        _alpha = alpha;
    }
    sinc_resampler * clone() const{
        sinc_resampler * resampler = new sinc_resampler(_domain, _range, _hl);
        return resampler;
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
};

// Matrix multiplication by a vector; the matrix is stored in row major
class matrix : public loper {
protected:
    std::shared_ptr<vecReg<data_t> > _mat;
public:
    matrix(){}
    ~matrix(){}
    matrix(const std::shared_ptr<vecReg<data_t> > mat, const hypercube<data_t> &domain, const hypercube<data_t> &range){
        successCheck(mat->getHyper()->getNdim()==2,"The matrix must be 2D\n");
        successCheck(mat->getHyper()->getAxis(1).n==domain.getN123(),"The domain dimension must match the matrix row size\n");
        successCheck(mat->getHyper()->getAxis(2).n==range.getN123(),"The range dimension must match the matrix column size\n");
        _domain = domain;
        _range = range;
        _mat=mat->clone();
    }
    matrix(const std::shared_ptr<vecReg<data_t> > mat){
        successCheck(mat->getHyper()->getNdim()==2,"The matrix must be 2D\n");
        _mat=mat->clone();
        _domain = hypercube<data_t>(mat->getHyper()->getAxis(1));
        _range = hypercube<data_t>(mat->getHyper()->getAxis(2));
    }
    matrix(int n, data_t val){
        _domain = hypercube<data_t>(n);
        _range = hypercube<data_t>(n);
        _mat = std::make_shared<vecReg<data_t> >(hypercube<data_t>(n,n));
        for (int i=0; i<n; i++) _mat->getVals()[i*n+i] = val;
    }
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        successCheck(_domain.getN123()==domain.getN123(),"The domain dimension must match the matrix row size\n");
        successCheck(_range.getN123()==range.getN123(),"The range dimension must match the matrix column size\n");
        _domain=domain;
        _range=range;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        if ((_domain.getN123() != mod->getN123()) || (_range.getN123() != dat->getN123())) return false;
        else return true;
    }
    matrix * clone() const {
        matrix * mat = new matrix(_mat,_domain,_range);
        return mat;
    }
    std::shared_ptr<vecReg<data_t> > getMat() const {return _mat;}
    
    void apply_forward(bool add, const data_t * pmod, data_t * pdat){
        const data_t * pmat = _mat->getCVals();
        int n=_domain.getN123();
        int m=_range.getN123();
        data_t val;
        #pragma omp parallel for
        for (int i=0; i<m; i++){
            val=0;
            for (int j=0; j<n; j++){
                val += pmat[i*m+j]*pmod[j]; 
            }
            pdat[i] = add*pdat[i] + val;
        } 
    }
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){
        const data_t * pmat = _mat->getCVals();
        int n=_domain.getN123();
        int m=_range.getN123();
        data_t val;
        #pragma omp parallel for
        for (int j=0; j<n; j++){
            val=0;
            for (int i=0; i<m; i++){
                val += pmat[i*m+j]*pdat[i]; 
            }
            pmod[j] = add*pmod[j] + val;
        }
    }
};

// Integral operator along the fast axis using the Trapezoidal quadrature
// The add option is inactive
class integral : public loper {
public:
    integral(){}
    ~integral(){}
    integral(const hypercube<data_t> &domain){
        _domain = domain;
        _range = domain;
    }
    integral * clone() const {
        integral * op = new integral(_domain);
        return op;
    }    
    void apply_forward(bool add, const data_t * pmod, data_t * pdat){
        int nt = _domain.getAxis(1).n;
        int nx = _domain.getN123() / nt;
        data_t dt = _domain.getAxis(1).d;
        #pragma omp parallel for
        for (int ix=0; ix<nx; ix++){
            int i1=ix*nt;
            pdat[i1] = 0.5*pmod[i1]*dt;
            for (int it=1; it<nt-1; it++){
                pdat[i1+it] = pdat[i1+it-1] + pmod[i1+it]*dt;
            }
            pdat[i1+nt-1] = pdat[i1+nt-2] + 0.5*pmod[i1+nt-1]*dt;
        }
    }
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){
        int nt = _domain.getAxis(1).n;
        int nx = _domain.getN123() / nt;
        data_t dt = _domain.getAxis(1).d;
        for (int ix=0; ix<nx; ix++){
            int i1=ix*nt;
            pmod[i1+nt-1] = pdat[i1+nt-1]*dt;
            for (int it=nt-2; it>=0; it--){
                pmod[i1+it] = pmod[i1+it+1] + pdat[i1+it]*dt;
            }
            pmod[i1+nt-1] *= 0.5;
            pmod[i1] *= 0.5; 
        }
    }
};


// class to tranform from t-x (or z-x) to f-x space and vice versa
class fxTransform{
protected:
    hypercube<data_t> _domain;
    hypercube<data_t> _range;
    hypercube<data_t> _crange;

public:
	fxTransform(){}
    ~fxTransform(){}
    fxTransform(const hypercube<data_t> &domain){
        _domain = domain;
        std::vector<axis<data_t> > axes = domain.getAxes();
        axis<data_t> Z = axes[0];
        Z.n = Z.n/2 + 1;
        Z.d = 1.0/((domain.getAxis(1).n-1)*Z.d);
        Z.o = 0.0;
        axes[0] = Z;
        _range = hypercube<data_t>(axes);
        Z.n = domain.getAxis(1).n;
        Z.o = -(Z.n-1)/2 * Z.d;
        axes[0] = Z;
        _crange = hypercube<data_t> (axes);
    }
    fxTransform * clone() const{
        fxTransform * op = new fxTransform(_domain);
        return op;
    }
    const hypercube<data_t> * getDomain() const {return &_domain;}
    const hypercube<data_t> * getRange() const {return &_range;}
    const hypercube<data_t> * getCRange() const {return &_crange;}
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        successCheck(domain.getN123()==_domain.getN123(),"The new domain must have the same number of samples\n");
        successCheck(range.getN123()==_range.getN123(),"The new range must have the same number of samples\n");
        successCheck(domain.getAxis(1)==_domain.getAxis(1),"The first axis in the domain must be the same\n");
        successCheck(range.getAxis(1)==_range.getAxis(1),"The first axis in the range must be the same\n");
        _domain = domain;
        _range = range;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat) const {
        bool ans1 = (mod->getN123()==_domain.getN123());
        bool ans2 = (dat->getN123()==_range.getN123());
        bool ans3 = (mod->getHyper()->getAxis(1)==_domain.getAxis(1));
        bool ans4 = (dat->getHyper()->getAxis(1)==_range.getAxis(1));
        return (ans1 && ans2 && ans3 && ans4);
    }
    bool checkCDomainRange(const std::shared_ptr<cvecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat) const {
        bool ans1 = (mod->getN123()==_domain.getN123());
        bool ans2 = (dat->getN123()==_crange.getN123());
        bool ans3 = (mod->getHyper()->getAxis(1)==_domain.getAxis(1));
        bool ans4 = (dat->getHyper()->getAxis(1)==_crange.getAxis(1));
        return (ans1 && ans2 && ans3 && ans4);
    }
    void apply_forward(bool add, const data_t * pmod, std::complex<data_t> * pdat);
    void apply_adjoint(bool add, data_t * pmod, const std::complex<data_t> * pdat){
        apply_inverse(add, pmod, pdat);
    }
    void apply_inverse(bool add, data_t * pmod, const std::complex<data_t> * pdat);
    void forward(bool add, const std::shared_ptr<vecReg<data_t> > mod, std::shared_ptr<cvecReg<data_t> > dat){
        successCheck(checkDomainRange(mod,dat),"Vectors hypercube do not match the operator domain and range\n");
        apply_forward(add, mod->getCVals(), dat->getVals());
    }
    void inverse(bool add, std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat){
        successCheck(checkDomainRange(mod,dat),"Vectors hypercube do not match the operator domain and range\n");
        apply_inverse(add, mod->getVals(), dat->getCVals());
    }
    void adjoint(bool add, std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat){
        inverse(add, mod, dat);
    }
    void cforward(bool add, const std::shared_ptr<cvecReg<data_t> > mod, std::shared_ptr<cvecReg<data_t> > dat);
    void cinverse(bool add, std::shared_ptr<cvecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat);
    void cadjoint(bool add, std::shared_ptr<cvecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat){
        cinverse(add, mod, dat);
    }
};


// low order 3D gradient operator to be used in first order Tikhonov regularization
// the boundary gradient is discarded
class gradient3d : public loper {
protected:
    data_t _xw, _yw, _zw; // directional weights
public:
    gradient3d(){}
    ~gradient3d(){}
    gradient3d(const hypercube<data_t> &domain, data_t xw=1, data_t yw=1, data_t zw=1){
        successCheck(domain.getNdim()>=3,"The domain must contain at least 3 dimensions\n");
        _domain = domain;
        std::vector<axis<data_t> > axes = domain.getAxes();
        axes.push_back(axis<data_t>(3,0,1));
        _range = hypercube<data_t>(axes);
        _xw=xw;
        _yw=yw;
        _zw=zw;
    }
    gradient3d * clone() const {
        gradient3d * op = new gradient3d(_domain);
        return op;
    }    
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
};

// low order 3D laplacian operator to be used in second order Tikhonov regularization
// the boundary laplacian is discarded
class laplacian3d : public loper {
protected:
    data_t _xw, _yw, _zw; // directional weights
public:
    laplacian3d(){}
    ~laplacian3d(){}
    laplacian3d(const hypercube<data_t> &domain, data_t xw=1, data_t yw=1, data_t zw=1){
        successCheck(domain.getNdim()>=3,"The domain must contain at least 3 dimensions\n");
        _domain = domain;
        _range = domain;
        _xw=xw;
        _yw=yw;
        _zw=zw;
    }
    laplacian3d * clone() const {
        laplacian3d * op = new laplacian3d(_domain);
        return op;
    }    
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
};

// Model extrapolation from 1D to 3D along a given horizon in the second and third dimensions
class extrapolator1d3d : public loper {
protected:
    std::shared_ptr<vecReg<data_t> > _hrz; // 2D horizon (surface) for extrapolation
public:
    extrapolator1d3d(){}
    ~extrapolator1d3d(){}
    extrapolator1d3d(const hypercube<data_t> &domain, std::shared_ptr<vecReg<data_t> > hrz){
        successCheck(domain.getNdim()==2,"Domain must contain 2 axes (depth and components)\n");
        successCheck(hrz->getHyper()->getNdim()==2,"The horizon must be 2D\n");
        _domain = domain;
        std::vector<axis<data_t> > axes = domain.getAxes();
        axis<data_t> X = hrz->getHyper()->getAxis(1);
        axis<data_t> Y = hrz->getHyper()->getAxis(2);
        axes.insert(axes.begin()+1, X);
        axes.insert(axes.begin()+2, Y);
        _range = hypercube<data_t>(axes);
        _hrz = hrz->clone();
    }
    extrapolator1d3d * clone() const {
        extrapolator1d3d * op = new extrapolator1d3d(_domain, _hrz);
        return op;
    }    
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
    void apply_inverse(bool add, data_t * pmod, const data_t * pdat);
};

// Time domain Convolution of 1D filter with N-dimensional vector
class conv1dnd : public loper {
protected:
    std::shared_ptr<vecReg<data_t> > _f; // filter always treated as 1D
    bool _centered; // center the convolution

public:
    conv1dnd(){}
    ~conv1dnd(){}
    conv1dnd(const hypercube<data_t> domain, const std::shared_ptr<vecReg<data_t> > f, bool centered = true){
        successCheck(f->getHyper()->getAxis(1).d==domain.getAxis(1).d,"Filter must have the same sampling as the domain fast axis\n");
        _domain=domain;
        _range=domain;
        _f=f->clone();
        _centered = centered;
    }
    conv1dnd * clone() const {
        conv1dnd * op = new conv1dnd(_domain,_f,_centered);
        return op;
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
};

// transform a filter to zero phase
std::shared_ptr<vecReg<data_t> > zero_phase(const std::shared_ptr<vecReg<data_t> > dat);
// transform a filter to minimum phase
std::shared_ptr<vecReg<data_t> > minimum_phase(const std::shared_ptr<vecReg<data_t> > dat, const data_t eps);