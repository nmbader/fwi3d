#pragma once

#define ZERO 1e-16
#define M_INF -1e+16

#include <time.h>
#include "operator.hpp"
#include "we_op.hpp"
#include "mpiWrapper.hpp"

#ifdef CUDA
    #include "cudaMisc.h"
#endif

// General class for optimization problems: min(f(m))
class optimization {
protected:
    std::shared_ptr<vecReg<data_t> > _d; // data vector
    std::shared_ptr<vecReg<data_t> > _m; // model vector
    std::shared_ptr<vecReg<data_t> > _r; // residual vector
    std::shared_ptr<vecReg<data_t> > _g; // functional gradient
public:
    optimization(){}
    virtual ~optimization(){}

    // default functional = 1/2 ||r||^2
    virtual data_t getFunc() {
        return 0.5*_r->norm2();
    }
    void initGrad() {
        _g = std::make_shared<vecReg<data_t> > (*_m->getHyper());
        _g->zero();
    }
    virtual void initRes() {
        _r = std::make_shared<vecReg<data_t> > (*_d->getHyper());
        _r->zero();
    }
    std::shared_ptr<vecReg<data_t> > getMod(){return _m;}
    std::shared_ptr<vecReg<data_t> > getGrad(){return _g;}
    std::shared_ptr<vecReg<data_t> > getRes(){return _r;}
    virtual void res() {}
    virtual void grad() {}
    virtual void hessian(const std::shared_ptr<vecReg<data_t> > H){};
    virtual data_t getZero(){return ZERO;}
};


// #################################################### Non-linear part ################################################################## //
// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// //
// #################################################### ############### ################################################################## //

// Rosenbrock function: f(x,y)=(a-x)^2 + b(y-x^2)^2     m = (x,y)
// minimizer at (a,a^2) and fmin=0
// grad(f)=( 2(x-a) + 4bx(x^2-y) , 2b(y-x^2) )
class rosenbrock : public optimization{
protected:
    data_t _a, _b;
public:
    rosenbrock(){}
    virtual ~rosenbrock(){}
    rosenbrock(data_t a, data_t b, std::shared_ptr<vecReg<data_t> > m){
        successCheck(m->getN123()==2,"Rosenbrock model must contain two elements\n");
        _a = a;
        _b = b;
        _m = m;
        _d = std::make_shared<vecReg<data_t> >(hypercube<data_t>(1));
        _d->set(0);
        initGrad();
        initRes();
    }
    // f(m)
    data_t getFunc(){
        data_t x = _m->getVals()[0];
        data_t y = _m->getVals()[1];
        return (_a-x)*(_a-x) + _b*(y-x*x)*(y-x*x);
    }
    // grad(f)
    void grad() {
        data_t x = _m->getVals()[0];
        data_t y = _m->getVals()[1];
        _g->getVals()[0] = 2*(x-_a) + 4*_b*x*(x*x-y);
        _g->getVals()[1] = 2*_b*(y-x*x); 
    }
    // Hessian(f)
    void hessian(std::shared_ptr<vecReg<data_t> > H){
        data_t x = _m->getVals()[0];
        data_t y = _m->getVals()[1];
        H->getVals()[0] = 2+4*_b*(x*x-y) + 8*_b*x*x;
        H->getVals()[1] = -4*_b*x;
        H->getVals()[2] = H->getVals()[1];
        H->getVals()[3] = 2*_b;
    }
};

// Paraboloid: f(x,y)=1/2(a.x)^2 + 1/2(b.y)^2     m = (x,y)
// minimizer at (0,0) and fmin=0
// grad(f)=(a^2.x,b^2.y)
class paraboloid : public optimization{
protected:
    data_t _a, _b;
public:
    paraboloid(){}
    virtual ~paraboloid(){}
    paraboloid(data_t a, data_t b, std::shared_ptr<vecReg<data_t> > m){
        successCheck(m->getN123()==2,"Paraboloid model must contain two elements\n");
        _a = a;
        _b = b;
        _m = m;
        _d = std::make_shared<vecReg<data_t> >(hypercube<data_t>(1));
        _d->set(0);
        initGrad();
        initRes();
    }
    // f(m)
    data_t getFunc(){
        data_t x = _m->getVals()[0];
        data_t y = _m->getVals()[1];
        return 0.5*((_a*_a*x*x) + (_b*_b*y*y));
    }
    // grad(f)
    void grad() {
        data_t x = _m->getVals()[0];
        data_t y = _m->getVals()[1];
        _g->getVals()[0] = _a*_a*x;
        _g->getVals()[1] = _b*_b*y; 
    }
    // Hessian(f)
    void hessian(std::shared_ptr<vecReg<data_t> > H){
        H->getVals()[0] = _a*_a;
        H->getVals()[1] = 0;
        H->getVals()[2] = 0;
        H->getVals()[3] = _b*_b;
    }
};

// Non-linear least-squares problem: f(m)=1/2.|L(m)-d|^2
class nlls : public optimization{
protected:
    nloper * _L; // non-linear operator

public:
    nlls(){}
    virtual ~nlls(){}
    nlls(nloper * L, std::shared_ptr<vecReg<data_t> > m, std::shared_ptr<vecReg<data_t> > d){
       successCheck(L->checkDomainRange(m,d),"Vectors hypercube do not match the operator domain and range\n");
        _L = L;
        _m = m;
        _d = d;
        initGrad();
        initRes();
    }
    // L(m) - d
    virtual void res(){
        _L->forward(false,_m,_r);
        _r->scaleAdd(_d,1,-1);
    }
    // (dL(m)/dm)'.(Lm-d)
    virtual void grad() {
        _L->jacobianT(false,_g,_m,_r);
    }
};

// Non-linear least-squares problem with time quadrature, model preconditioner, data weighting and gradient mask for FWI
// the functional is normalized: f(m)=1/2.[(L(m)-d)'.Ht.(L(m)-d)]/(d'.Ht.d) = 1/2.r'.Ht.r / (d'.Ht.d)
// optionally other functional may be considered:
// 1) with weighting: r <- W.r
// 2) with trace normalization: r <- (u/|u| - d/|d|) where u = L(m)
// 3) envelop or envelop squared: r <- E(u) - E(d)
// The WE operator and model preconditioner are provided separately
class nlls_fwi : public optimization{
protected:
    nl_we_op * _L; // non-linear WE operator
    nloper * _P; // model preconditioner
    std::shared_ptr<vecReg<data_t> > _p; // preconditioned model
    std::shared_ptr<vecReg<data_t> > _pg; // gradient before applying preconditioner Jacobian
    std::shared_ptr<vecReg<data_t> > _gmask; // gradient mask vector
    std::shared_ptr<vecReg<data_t> > _w; // residual weighting vector
    std::shared_ptr<vecReg<data_t> > _filter; // 1D filter
    data_t _dnorm; // data normalization factor
    data_t _f; // objective function
    bool _flag; // flag to re-evaluate or not the objective function
    int _scale_source_times; // for source inversion by reduced Variable Projection
    param _par; // parameters object for the non-linear WE operator

public:
    std::vector <data_t> _dfunc; // vector that stores data objective function values for all trials (used in the regularization class)
    std::vector <data_t> _mfunc; // vector that stores model objective function values for all trials (used in the regularization class)
    std::vector<data_t> _scalers; // scalers applied to each source at a given trial
    nlls_fwi(){}
    virtual ~nlls_fwi(){}
    nlls_fwi(nl_we_op * L, param &par, std::shared_ptr<vecReg<data_t> > m, std::shared_ptr<vecReg<data_t> > d, nloper * P = nullptr, std::shared_ptr<vecReg<data_t> > gmask = nullptr, std::shared_ptr<vecReg<data_t> > w = nullptr, std::shared_ptr<vecReg<data_t> > filter = nullptr){
        _L = L;
        _par = par;
        _P = P;
        _m = m;
        _d = d;
        _gmask = gmask;
        _w = w;
        _filter = filter;
        _f = 0;
        _flag = true;
        initGrad();
        initRes();
       
        // Model precconditioner (or parameterizer)
        if (P!= nullptr){
            _p = std::make_shared<vecReg<data_t> >(*P->getRange());
            _pg = std::make_shared<vecReg<data_t> >(*P->getRange());
            _p->zero();
            _pg->zero();
            successCheck(P->checkDomainRange(m,_p),"Vectors hypercube do not match the operator domain and range\n");
        }
        else{
            _p = _m;
            _pg = _g;
        }
        if (gmask != nullptr) {
            successCheck(m->getN123()==gmask->getN123(),"The gradient mask and model vectors must have the same size\n");
            successCheck(gmask->min()>=0,"Gradient mask must be non-negative\n");
        }
        if (w != nullptr) {
            successCheck(d->getN123()==w->getN123(),"The weighting and data vectors must have the same size\n");
            successCheck(w->min()>=0,"Data weights  must be non-negative\n");
        }

        if (_w != nullptr) _d->mult(_w); // apply weighting

        if (_filter != nullptr){ // apply filtering
            data_t eps=1e-07;
            if (_par.filter_phase=="zero") _filter = zero_phase(filter);
            else if (_par.filter_phase=="minimum") _filter = minimum_phase(filter,eps);
            axis<data_t> Tf = _filter->getHyper()->getAxis(1);
            axis<data_t> T = _d->getHyper()->getAxis(1);
            Tf.d=T.d;
            _filter->setHyper(hypercube<data_t>(Tf)); 
            conv1dnd op(*_d->getHyper(), _filter, _par.filter_phase!="minimum");
            std::shared_ptr<vecReg<data_t> > output = std::make_shared<vecReg<data_t> >(*op.getRange());
            output->zero();
            op.forward(false,_d,output);
            _d=output;
        }

        if (_par.normalize) { // apply trace normalization
            int nt=_d->getHyper()->getAxis(1).n;
            int ntr=_d->getN123()/nt;
            data_t * norms = new data_t[ntr];
            ttnormalize(_d->getVals(), norms, nt, ntr);
            delete [] norms;
        }
        if (_par.envelop==1) envelop1(_d); // compute the envelop (or envelop squared) of the data
        else if (_par.envelop==2) envelop2(_d);

        if (_par.normalize_obj_func){
            _dnorm = _d->norm2();
            _dnorm *= _d->getHyper()->getAxis(1).d;
            if (_dnorm<ZERO) _dnorm=1;
        }
        else _dnorm=1.0;

        _scale_source_times=0;
        if (_par.scale_source_times>0) {
            for (int s=0; s<_par.ns; s++) _scalers.push_back(1.0);
        }
    }

    void compute_res_and_grad(data_t * r);

    virtual void res(){
        _r->zero();
        compute_res_and_grad(_r->getVals()); _flag = true;
        // gather all residual
        mpiWrapper::allReduceSum(_r->getCVals(), _r->getVals(), _r->getN123());
    }

    virtual data_t getFunc(){
        if (_flag)
        {
            int n123 = _r->getN123();
            int nt = _r->getHyper()->getAxis(1).n;
            int nx = n123 / nt;
            data_t dt = _r->getHyper()->getAxis(1).d;
            data_t * pr = _r->getVals();
            data_t f = _r->norm2();
            for (int i=0; i<nx; i++) f = f - 0.5*(pr[i*nt]*pr[i*nt] + pr[i*nt+nt-1]*pr[i*nt+nt-1]);
            _f = 0.5*dt*f/_dnorm;
            _flag = false;
        }
        return _f;
    }
    virtual void grad() {
        // already computed in the compute_res_and_grad() method above
        if (_P != nullptr) _P->jacobianT(false,_g,_m,_pg);
        if (_gmask != nullptr) _g->mult(_gmask);
    }
};

// Same class as nlls_fwi but with model regularizations
// f(m)=1/2.[(L(m)-d)'.Ht.(L(m)-d)]/(d'.Ht.d)  + 1/2.lambda^2.|D(m)-D(m_prior)|^2 / |D(m_prior)|^2
class nlls_fwi_reg : public nlls_fwi{
protected:
    nloper * _D; // regularization operator
    data_t _lambda; // regularization damping parameter
    std::shared_ptr<vecReg<data_t> > _Dmp; // prior model pre-multiplied by D
    std::shared_ptr<vecReg<data_t> > _dg; // gradient component corresponding to the regularization D
    data_t _mnorm; // model normalization factor

public:
    nlls_fwi_reg(){}
    virtual ~nlls_fwi_reg(){}
    nlls_fwi_reg(nl_we_op * L, param &par, nloper * D, std::shared_ptr<vecReg<data_t> > m, std::shared_ptr<vecReg<data_t> > d, data_t lambda, std::shared_ptr<vecReg<data_t> > mprior = nullptr, nloper * P = nullptr, std::shared_ptr<vecReg<data_t> > gmask = nullptr, std::shared_ptr<vecReg<data_t> > w = nullptr, std::shared_ptr<vecReg<data_t> > filter = nullptr){
        _L = L;
        _par = par;
        _D = D;
        _P = P;
        _m = m;
        _d = d;
        _lambda = lambda;    
        _gmask = gmask;
        _w = w;
        _filter = filter;
        _f = 0;
        _flag = true;
        initGrad();
        initRes();

        _Dmp = std::make_shared<vecReg<data_t> > (*_D->getRange());
        _Dmp->zero();
        if (mprior!=nullptr) _D->forward(false, mprior, _Dmp); // mprior assumed = 0 if not provided

        _dg = std::make_shared<vecReg<data_t> >(*m->getHyper());
        _dg->zero();
       
        if (P!= nullptr){
            _p = std::make_shared<vecReg<data_t> >(*P->getRange());
            _pg = std::make_shared<vecReg<data_t> >(*P->getRange());
            _p->zero();
            _pg->zero();
            successCheck(P->checkDomainRange(m,_p),"Vectors hypercube do not match the operator domain and range\n");
        }
        else{
            _p = _m;
            _pg = _g;
        }
        if (gmask != nullptr) {
            successCheck(m->getN123()==gmask->getN123(),"The gradient mask and model vectors must have the same size\n");
            successCheck(gmask->min()>=0,"Gradient mask must be non-negative\n");
        }
        if (w != nullptr) {
            successCheck(d->getN123()==w->getN123(),"The weighting and data vectors must have the same size\n");
            successCheck(w->min()>=0,"Data weights  must be non-negative\n");
        }

        if (_w != nullptr) _d->mult(_w);

        if (_filter != nullptr){ // apply filtering
            data_t eps=1e-07;
            if (_par.filter_phase=="zero") _filter = zero_phase(filter);
            else if (_par.filter_phase=="minimum") _filter = minimum_phase(filter,eps);
            axis<data_t> Tf = _filter->getHyper()->getAxis(1);
            axis<data_t> T = _d->getHyper()->getAxis(1);
            Tf.d=T.d;
            _filter->setHyper(hypercube<data_t>(Tf)); 
            conv1dnd op(*_d->getHyper(), _filter, _par.filter_phase!="minimum");
            std::shared_ptr<vecReg<data_t> > output = std::make_shared<vecReg<data_t> >(*op.getRange());
            output->zero();
            op.forward(false,_d,output);
            _d=output;
        }

        if (_par.normalize) {
            int nt=_d->getHyper()->getAxis(1).n;
            int ntr=_d->getN123()/nt;
            data_t * norms = new data_t[ntr];
            ttnormalize(_d->getVals(), norms, nt, ntr);
            delete [] norms;
        }
        if (_par.envelop==1) envelop1(_d);
        else if (_par.envelop==2) envelop2(_d);

        if (L->_par.normalize_obj_func){
            _dnorm = _d->norm2();
            _dnorm *= _d->getHyper()->getAxis(1).d;
            if (_dnorm<ZERO) _dnorm=1;

            _mnorm = _Dmp->norm2();
            if (_mnorm<ZERO) _mnorm=1;
            //_mnorm = std::max((data_t)1.0, _mnorm);
        }
        else{
            _dnorm=1.0;
            _mnorm=1.0;
        }

        _scale_source_times=0;
        if (_par.scale_source_times>0) {
            for (int s=0; s<_par.ns; s++) _scalers.push_back(1.0);
        }
    }

    void initRes() {
        _r = std::make_shared<vecReg<data_t> >(hypercube<data_t>(_d->getN123()+_D->getRange()->getN123()));
        _r->zero();
        _dfunc.clear();
        _mfunc.clear();
    }

    // res = ( L(m) - d ; lambda*(D(m) - D(m_prior)) )
    // grad = (dL(m)/dm)'.Ht.(L(m)-d)/(d'.Ht.d) + lambda.(dD(m)/dm)'.(D(m)-D(m_prior))/|D(m_prior)|^2
    void res(){
        _r->zero();
        int nd = _d->getN123();
        int nm = _D->getRange()->getN123();
        data_t * pr = _r->getVals();
        const data_t * pd = _d->getCVals();
        const data_t * pdmp = _Dmp->getCVals();
        compute_res_and_grad(_r->getVals());

        // gather all residual
        mpiWrapper::allReduceSum(_r->getCVals(), _r->getVals(), _d->getN123());

        _D->apply_forward(false,_m->getCVals(),pr+nd);
        #pragma omp parallel for
        for (int i=0; i<nm; i++) pr[nd+i] = _lambda*(pr[nd+i] - pdmp[i]);

        _flag = true;
    }

    data_t getFunc(){
        if (_flag)
        {
            int n123 = _d->getN123();
            int nt = _d->getHyper()->getAxis(1).n;
            int nx = n123 / nt;
            data_t dt = _d->getHyper()->getAxis(1).d;
            data_t * pr = _r->getVals();
            data_t fd = _r->norm2(0,n123);
            data_t fm = _r->norm2(n123,-1);
            for (int i=0; i<nx; i++) fd = fd - 0.5*(pr[i*nt]*pr[i*nt] + pr[i*nt+nt-1]*pr[i*nt+nt-1]);
            fd = 0.5*dt*fd/_dnorm;
            fm = 0.5*fm/_mnorm;
            _f = fd+fm;
            _flag = false;
            
            if (_par.verbose>0) fprintf(stderr,"Data functional = %f; Model functional = %f\n",fd,fm);
            _dfunc.push_back(fd);
            _mfunc.push_back(fm);
        }
        return _f;
    }
    void grad() {
        // the first component is already computed in the compute_res_and_grad() method above
        _D->apply_jacobianT(false,_dg->getVals(),_m->getCVals(),_r->getVals()+_d->getN123());

        if (_P != nullptr) _P->jacobianT(false,_g,_m,_pg);
        _g->scaleAdd(_dg,1,_lambda/_mnorm);

        if (_gmask != nullptr) _g->mult(_gmask);
    }
};

#undef ZERO
#undef M_INF