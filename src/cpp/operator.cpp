#include "operator.hpp"
#include "fftw3.h"

void nloper::dotProduct(){

    std::shared_ptr<vecReg<data_t> > m (new vecReg<data_t>(_domain));
    std::shared_ptr<vecReg<data_t> > dm (new vecReg<data_t>(_domain));
    std::shared_ptr<vecReg<data_t> > d (new vecReg<data_t>(_range));
    std::shared_ptr<vecReg<data_t> > m1;

    m->random(-1,1,123);
    d->random(-1,1,123);
    dm->random(-1,1,123);

    std::shared_ptr<vecReg<data_t> > mtilde = m->clone();
    std::shared_ptr<vecReg<data_t> > d1 = d->clone();
    std::shared_ptr<vecReg<data_t> > dtilde = d->clone();

    data_t eps = 0.1;

    m1 = dm->clone();
    m1->scaleAdd(m,eps,1); // m + eps.dm
    forward(false, m, d1); // f(m)
    forward(false,m1,dtilde); // f(m + eps.dm)
    dtilde->scaleAdd(d1,1,-1); // f(m + eps.dm) - f(m)
    dtilde->scale(1.0/eps); // (df/dm).dm
    jacobianT(false,mtilde,m,d); // (df/dm)'.d

    data_t sum1 = d->dot(dtilde);
    data_t sum2 = dm->dot(mtilde);

    fprintf(stderr,"Dot product with eps = %f:\n",eps);
    fprintf(stderr,"sum1 = %f\nsum2 = %f\ndiff (x1000) = %f\n",sum1,sum2,1000*(sum1 - sum2));

    eps = 0.01;

    m1 = dm->clone();
    m1->scaleAdd(m,eps,1); // m + eps.dm
    forward(false,m1,dtilde); // f(m + eps.dm)
    dtilde->scaleAdd(d1,1,-1); // f(m + eps.dm) - f(m)
    dtilde->scale(1.0/eps); // (df/dm).dm
    jacobianT(false,mtilde,m,d); // (df/dm)'.d

    sum1 = d->dot(dtilde);
    sum2 = dm->dot(mtilde);

    fprintf(stderr,"Dot product with eps = %f:\n",eps);
    fprintf(stderr,"sum1 = %f\nsum2 = %f\ndiff (x1000) = %f\n",sum1,sum2,1000*(sum1 - sum2));

    eps = 0.0001;
    
    m1 = dm->clone();
    m1->scaleAdd(m,eps,1); // m + eps.dm
    forward(false,m1,dtilde); // f(m + eps.dm)
    dtilde->scaleAdd(d1,1,-1); // f(m + eps.dm) - f(m)
    dtilde->scale(1.0/eps); // (df/dm).dm
    jacobianT(false,mtilde,m,d); // (df/dm)'.d

    sum1 = d->dot(dtilde);
    sum2 = dm->dot(mtilde);

    fprintf(stderr,"Dot product with eps = %f:\n",eps);
    fprintf(stderr,"sum1 = %f\nsum2 = %f\ndiff (x1000) = %f\n",sum1,sum2,1000*(sum1 - sum2));   
}

void loper::dotProduct(){

    std::shared_ptr<vecReg<data_t> > m (new vecReg<data_t>(_domain));
    std::shared_ptr<vecReg<data_t> > d (new vecReg<data_t>(_range));

    //time_t timer; unsigned int seed=time(&timer);
    m->random(-1,1,123);
    d->random(-1,1,123);

    std::shared_ptr<vecReg<data_t> > mtild = m->clone();
    std::shared_ptr<vecReg<data_t> > dtild = d->clone();

    this->forward(false, m, dtild);
    this->adjoint(false, mtild, d);

    data_t sum1 = d->dot(dtild);
    data_t sum2 = m->dot(mtild);

    fprintf(stderr,"Dot product with add=false:\n");
    fprintf(stderr,"sum1 = %f\nsum2 = %f\ndiff (x1000) = %f\n",sum1,sum2,1000*(sum1 - sum2));

    this->forward(true, m, dtild);
    this->adjoint(true, mtild, d);

    sum1 = d->dot(dtild);
    sum2 = m->dot(mtild);

    fprintf(stderr,"Dot product with add=true:\n");
    fprintf(stderr,"sum1 = %f\nsum2 = %f\ndiff (x1000) = %f\n",sum1,sum2,1000*(sum1 - sum2));
}

void softClip::apply_forward(bool add, const data_t * pmod, data_t * pdat) {
    
    int n = _domain.getN123();
    data_t xbar = 0.5*(_xmax+_xmin);
    data_t dx = 2*(_xmax-xbar)*std::pow(1.0*_q/(_q-1),1.0/_p);
    data_t xmax = xbar + dx;
    data_t xmin = xbar - dx;
    data_t g = 0;

    #pragma omp parallel for private(g)
    for (int i=0; i<n; i++)
    {
        g = std::pow(2*(pmod[i]-xbar)/dx,_p);
        if (std::abs(g) < 1) {
            if (g>0) pdat[i] = add*pdat[i] + xbar + 0.5*dx*std::pow(g-std::pow(g,_q)/_q,1.0/_p);
            else if (g<0) pdat[i] = add*pdat[i] + xbar - 0.5*dx*std::pow(-g+std::pow(g,_q)/_q,1.0/_p);
            else pdat[i] = add*pdat[i] + xbar;
        }
        else if (g>=1) pdat[i] = add*pdat[i] + xbar + 0.5*dx*std::pow(1.0*(_q-1)/_q,1.0/_p);
        else pdat[i] = add*pdat[i] + xbar - 0.5*dx*std::pow(1.0*(_q-1)/_q,1.0/_p);
    }
}
void softClip::apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat) {

    int n = _domain.getN123();
    data_t xbar = 0.5*(_xmax+_xmin);
    data_t dx = 2*(_xmax-xbar)*std::pow(1.0*_q/(_q-1),1.0/_p);
    data_t g = 0;
    #pragma omp parallel for private(g)
    for (int i=0; i<n; i++)
    {
        g = std::pow(2*(pmod0[i]-xbar)/dx,_p);
        if (std::abs(g) < 1 && g!=0) {
            if  (g>0) pmod[i] = add*pmod[i] + std::pow(2*(pmod0[i]-xbar)/dx,_p-1)*(1-std::pow(g,_q-1))*std::pow(g-std::pow(g,_q)/_q,1.0*(1-_p)/_p)*pdat[i];
            else pmod[i] = add*pmod[i] + std::pow(2*(pmod0[i]-xbar)/dx,_p-1)*(1-std::pow(g,_q-1))*std::pow(-g+std::pow(g,_q)/_q,1.0*(1-_p)/_p)*pdat[i];
        }
        else pmod[i] = add*pmod[i];
    }
}

void emodelSoftClip::apply_forward(bool add, const data_t * pmod, data_t * pdat){
    
    hypercube<data_t> hyper(_domain.getAxis(1),_domain.getAxis(2));
    softClip vpClip(hyper, _vpmin, _vpmax, _p, _q);
    softClip vsClip(hyper, _vsmin, _vsmax, _p, _q);
    softClip rhoClip(hyper, _rhomin, _rhomax, _p, _q);
    softClip spratioClip(hyper, -_spratio, _spratio, _p, _q);

    int n = hyper.getN123();

    data_t (* pm) [n] = (data_t (*) [n]) pmod;
    data_t (* pd) [n] = (data_t (*) [n]) pdat;

    data_t * v0;
    if (add) {
        v0 = new data_t [n];
        memcpy(v0, pd[1], sizeof(data_t)*n);
    }

    // vp, vs, rho clip
    vpClip.apply_forward(add, pm[0], pd[0]);
    vsClip.apply_forward(false, pm[1], pd[1]);
    rhoClip.apply_forward(add, pm[2], pd[2]);

    // vs/vp ratio clip
    #pragma omp parallel for
    for (int i=0; i<n; i++) pd[1][i] = pd[1][i] / pd[0][i];

    spratioClip.apply_forward(false, pd[1], pd[1]);
    
    #pragma omp parallel for
    for (int i=0; i<n; i++) {
        pd[1][i] = pd[1][i] * pd[0][i];
    }

    if (add){
        #pragma omp parallel for
        for (int i=0; i<n; i++) {
            pd[1][i] += v0[i];
        }
        delete [] v0;
    }

    // copy other parameters if any
    #pragma omp parallel for
    for (int i=3*n; i<_domain.getN123(); i++){
        pdat[i] = add*pdat[i] + pmod[i];
    }
}

void emodelSoftClip::apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat){

    hypercube<data_t> hyper(_domain.getAxis(1),_domain.getAxis(2));
    softClip vpClip(hyper, _vpmin, _vpmax, _p, _q);
    softClip vsClip(hyper, _vsmin, _vsmax, _p, _q);
    softClip rhoClip(hyper, _rhomin, _rhomax, _p, _q);
    softClip spratioClip(hyper, -_spratio, _spratio, _p, _q);

    int n = hyper.getN123();

    data_t (* pm) [n] = (data_t (*) [n]) pmod;
    const data_t (* pm0) [n] = (const data_t (*) [n]) pmod0;
    const data_t (* pd) [n] = (const data_t (*) [n]) pdat;

    // rho jacobian
    rhoClip.apply_jacobianT(add, pm[2], pm0[2], pd[2]);

    // adjoint for other parameters if any
    #pragma omp parallel for
    for (int i=3*n; i<_domain.getN123(); i++){
        pmod[i] = add*pmod[i] + pdat[i];
    }

    // forward vp,vs clip -> vp,spratio -> spratio clip
    data_t * v0 = new data_t [2*n];
    data_t * v1 = new data_t [n];
    data_t * m0 = new data_t [2*n];
    memset(v0, 0, sizeof(data_t)*2*n);
    memset(m0, 0, sizeof(data_t)*2*n);
    vpClip.apply_forward(false, pmod0, v0);
    vsClip.apply_forward(false, pmod0+n, v0+n);

    #pragma omp parallel for
    for (int i=0; i<n; i++) v0[n+i] = v0[n+i] / v0[i];

    memcpy(v1, v0+n, sizeof(data_t)*n);

    spratioClip.apply_forward(false, v1, v1);

    // jacobian of vp,spratio -> vp,vs
    #pragma omp parallel for
    for (int i=0; i<n; i++) {
        m0[i] = pd[0][i] + v1[i]*pd[1][i];
        m0[n+i] = v0[i]*pd[1][i];
    }

    // jacobian of spratio clip
    spratioClip.apply_jacobianT(false,m0+n, v0+n, m0+n);

    // jacobian of vp,vs -> vp,spratio
    #pragma omp parallel for
    for (int i=0; i<n; i++) {
        m0[i] = m0[i] - v0[n+i]/v0[i]*m0[n+i];
        m0[n+i] /= v0[i];
    }

    // jacobian of vp,vs clip
    vpClip.apply_jacobianT(add, pm[0], pm0[0], m0);
    vsClip.apply_jacobianT(add, pm[1], pm0[1], m0+n);

    delete [] v0;
    delete [] v1;
    delete [] m0;
}

inline data_t sinc(data_t x){
    if (x==0.0) return 1.0;
    else return sin(x)/x;
}

void linear_resampler::apply_forward(bool add, const data_t * pmod, data_t * pdat) {

    axis<data_t> Tm = _domain.getAxis(1);
    axis<data_t> Td = _range.getAxis(1);
    int nx = _domain.getN123()/Tm.n;

    int i, j;
    data_t wi, wj;
    #pragma omp parallel for private(i,j,wi,wj)
    for (int ix=0; ix<nx; ix++){
        for (int it=0; it<Td.n; it++){

            i = std::min(Tm.n - 1, (int) floor(it*Td.d / Tm.d));
            j = std::min(Tm.n - 1, (int) ceil(it*Td.d / Tm.d));
            wj = (it*Td.d - i*Tm.d)/Tm.d;
            wi = 1 - wj;
            pdat[ix*Td.n+it] = add*pdat[ix*Td.n+it] + wi*pmod[ix*Tm.n+i] + wj*pmod[ix*Tm.n+j];
        }
    }

}

void linear_resampler::apply_adjoint(bool add, data_t * pmod, const data_t * pdat) {

    axis<data_t> Tm = _domain.getAxis(1);
    axis<data_t> Td = _range.getAxis(1);
    int nx = _domain.getN123()/Tm.n;

    if (add == false) {
        for (int i=0; i<Tm.n*nx; i++) pmod[i] = 0;
    }

    int i, j;
    data_t wi, wj;
    #pragma omp parallel for private(i,j,wi,wj)
    for (int ix=0; ix<nx; ix++){
        for (int it=0; it<Td.n; it++){

            i = std::min(Tm.n - 1, (int) floor(it*Td.d / Tm.d));
            j = std::min(Tm.n - 1, (int) ceil(it*Td.d / Tm.d));
            wj = (it*Td.d - i*Tm.d)/Tm.d;
            wi = 1 - wj;
            pmod[ix*Tm.n+i] += wi*pdat[ix*Td.n+it];
            pmod[ix*Tm.n+j] += wj*pdat[ix*Td.n+it];
        }
    }
}

void sinc_resampler::apply_forward(bool add, const data_t * pmod, data_t * pdat) {

    axis<data_t> Tm = _domain.getAxis(1);
    axis<data_t> Td = _range.getAxis(1);
    int nx = _domain.getN123()/Tm.n;

    int itm0;
    data_t t;
    data_t a = _alpha; // 0 < a < 1 for the windowing of the cosine taper

    #pragma omp parallel for private(itm0,t)
    for (int ix=0; ix<nx; ix++){
        for (int itd=0; itd<Td.n; itd++){
            t = itd * Td.d;
            itm0 = floor(t / Tm.d);
            pdat[ix*Td.n+itd] = add*pdat[ix*Td.n+itd];
            for (int itm = std::max(0, itm0 -_hl + 1) ; itm <= std::min(Tm.n-1, itm0 + _hl); itm++){
                pdat[ix*Td.n+itd] += pmod[ix*Tm.n+itm]*sinc(M_PI * (t-itm*Tm.d) / Tm.d)*(a +(1-a)*cos(M_PI*(itm-t/Tm.d)/_hl));
            }
        }
    }

}

void sinc_resampler::apply_adjoint(bool add, data_t * pmod, const data_t * pdat) {

    axis<data_t> Tm = _domain.getAxis(1);
    axis<data_t> Td = _range.getAxis(1);
    int nx = _domain.getN123()/Tm.n;

    if (add == false) {
       for (int i=0; i<Tm.n*nx; i++) pmod[i] = 0;
    }

    int itm0;
    data_t t;
    data_t a = _alpha; // 0 < a < 1 for the windowing of the cosine taper
    
    #pragma omp parallel for private(itm0,t)
    for (int ix=0; ix<nx; ix++){
        for (int itd=0; itd<Td.n; itd++){
            t = itd * Td.d;
            itm0 = floor(t / Tm.d);
            for (int itm = std::max(0, itm0 -_hl + 1) ; itm <= std::min(Tm.n-1, itm0 + _hl); itm++){
                pmod[ix*Tm.n+itm] += pdat[ix*Td.n+itd]*sinc(M_PI * (t-itm*Tm.d) / Tm.d)*(a +(1-a)*cos(M_PI*(itm-t/Tm.d)/_hl));
            }
        }
    }
}

void fxTransform::apply_forward(bool add, const data_t * pmod, std::complex<data_t> * pdat)
{
    int nt = _domain.getAxis(1).n;
    int ntr = _domain.getN123()/nt;
    int nf = _range.getAxis(1).n;
    const data_t (*m) [nt] = (const data_t (*) [nt]) pmod;
    std::complex<data_t> (*d) [nf] = (std::complex<data_t> (*) [nf]) pdat;

    // double precision
    if ((int)sizeof(data_t)==8){

        double * in = new double[nt];
        fftw_complex * out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf);
        fftw_plan p = fftw_plan_dft_r2c_1d(nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, m[itr], sizeof(data_t)*nt);

            // forward transform
            fftw_execute(p);

            // copy to output
            for (int f=0; f<nf; f++) {
                d[itr][f].real(add*d[itr][f].real() + out[f][0]);
                d[itr][f].imag(add*d[itr][f].imag() + out[f][1]);
            }
        }

        delete [] in;
        fftw_destroy_plan(p);
        fftw_free(out);
        fftw_cleanup();
    }
    else
    {
        float * in = new float[nt];
        fftwf_complex * out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nf);
        fftwf_plan p = fftwf_plan_dft_r2c_1d(nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, m[itr], sizeof(data_t)*nt);

            // forward transform
            fftwf_execute(p);

            // copy to output
            for (int f=0; f<nf; f++) {
                d[itr][f].real(add*d[itr][f].real() + out[f][0]);
                d[itr][f].imag(add*d[itr][f].imag() + out[f][1]);
            }
        }

        delete [] in;
        fftwf_destroy_plan(p);
        fftwf_free(out);
        fftwf_cleanup();
    }
}

void fxTransform::apply_inverse(bool add, data_t * pmod, const std::complex<data_t> * pdat)
{

    int nt = _domain.getAxis(1).n;
    int ntr = _domain.getN123()/nt;
    int nf = _range.getAxis(1).n;
    data_t (*m) [nt] = (data_t (*) [nt]) pmod;
    const std::complex<data_t> (*d) [nf] = (const std::complex<data_t> (*) [nf]) pdat;

    // double precision
    if ((int)sizeof(data_t)==8){

        double * out = new double[nt];
        fftw_complex * in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf);
        fftw_plan p = fftw_plan_dft_c2r_1d(nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, d[itr], sizeof(fftw_complex)*nf);

            // inverse transform
            fftw_execute(p);

            // copy to output
            for (int it=0; it<nt; it++) m[itr][it]= add*m[itr][it] + out[it]/nt;
        }

        delete [] out;
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_cleanup();
    }
    else
    {
        float * out = new float[nt];
        fftwf_complex * in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nf);
        fftwf_plan p = fftwf_plan_dft_c2r_1d(nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, d[itr], sizeof(fftwf_complex)*nf);

            // inverse transform
            fftwf_execute(p);

            // copy to output
            for (int it=0; it<nt; it++) m[itr][it]= add*m[itr][it] + out[it]/nt;
        }

        delete [] out;
        fftwf_destroy_plan(p);
        fftwf_free(in);
        fftwf_cleanup();
    }
}

void fxTransform::cforward(bool add, const std::shared_ptr<cvecReg<data_t> > mod, std::shared_ptr<cvecReg<data_t> > dat)
{
    successCheck(checkCDomainRange(mod,dat),"Vectors hypercube do not match the operator domain and range\n");

    int nt = _domain.getAxis(1).n;
    int ntr = _domain.getN123()/nt;
    int nf = _crange.getAxis(1).n;
    const std::complex<data_t> (*m) [nt] = (const std::complex<data_t> (*) [nt]) mod->getCVals();
    std::complex<data_t> (*d) [nf] = (std::complex<data_t> (*) [nf]) dat->getVals();

    // double precision
    if ((int)sizeof(data_t)==8){

        fftw_complex * in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nt);
        fftw_complex * out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nt);
        fftw_plan p = fftw_plan_dft_1d(nt, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, m[itr], sizeof(fftw_complex)*nt);

            // forward transform
            fftw_execute(p);

            // copy to output
            for (int f=0; f<(nf-1)/2; f++) {
                d[itr][f].real(add*d[itr][f].real() + out[nf-(nf-1)/2+f][0]);
                d[itr][f].imag(add*d[itr][f].imag() + out[nf-(nf-1)/2+f][1]);
            }
            for (int f=(nf-1)/2; f<nf; f++) {
                d[itr][f].real(add*d[itr][f].real() + out[f-(nf-1)/2][0]);
                d[itr][f].imag(add*d[itr][f].imag() + out[f-(nf-1)/2][1]);
            }
        }

        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
        fftw_cleanup();
    }
    else
    {
        fftwf_complex * in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt);
        fftwf_complex * out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt);
        fftwf_plan p = fftwf_plan_dft_1d(nt, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, m[itr], sizeof(fftwf_complex)*nt);

            // forward transform
            fftwf_execute(p);

            // copy to output
            for (int f=0; f<(nf-1)/2; f++) {
                d[itr][f].real(add*d[itr][f].real() + out[nf-(nf-1)/2+f][0]);
                d[itr][f].imag(add*d[itr][f].imag() + out[nf-(nf-1)/2+f][1]);
            }
            for (int f=(nf-1)/2; f<nf; f++) {
                d[itr][f].real(add*d[itr][f].real() + out[f-(nf-1)/2][0]);
                d[itr][f].imag(add*d[itr][f].imag() + out[f-(nf-1)/2][1]);
            }
        }

        fftwf_destroy_plan(p);
        fftwf_free(in);
        fftwf_free(out);
        fftwf_cleanup();
    }
}

void fxTransform::cinverse(bool add, std::shared_ptr<cvecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat)
{
    successCheck(checkCDomainRange(mod,dat),"Vectors hypercube do not match the operator domain and range\n");

    int nt = _domain.getAxis(1).n;
    int ntr = _domain.getN123()/nt;
    int nf = _crange.getAxis(1).n;
    std::complex<data_t> (*m) [nt] = (std::complex<data_t> (*) [nt]) mod->getCVals();
    const std::complex<data_t> (*d) [nf] = (const std::complex<data_t> (*) [nf]) dat->getVals();

    // double precision
    if ((int)sizeof(data_t)==8){

        fftw_complex * in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nt);
        fftw_complex * out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nt);
        fftw_plan p = fftw_plan_dft_1d(nt, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in+nt-(nt-1)/2, d[itr], sizeof(fftw_complex)*(nt-1)/2);
            memcpy(in, d[itr]+(nt-1)/2, sizeof(fftw_complex)*(nt-(nt-1)/2));

            // forward transform
            fftw_execute(p);

            // copy to output
            for (int f=0; f<nf; f++) {
                m[itr][f].real(add*m[itr][f].real() + out[f][0]/nt);
                m[itr][f].imag(add*m[itr][f].imag() + out[f][1]/nt);
            }
        }

        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
        fftw_cleanup();
    }
    else
    {
        fftwf_complex * in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt);
        fftwf_complex * out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt);
        fftwf_plan p = fftwf_plan_dft_1d(nt, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in+nt-(nt-1)/2, d[itr], sizeof(fftwf_complex)*(nt-1)/2);
            memcpy(in, d[itr]+(nt-1)/2, sizeof(fftwf_complex)*(nt-(nt-1)/2));

            // forward transform
            fftwf_execute(p);

            // copy to output
            for (int f=0; f<nf; f++) {
                m[itr][f].real(add*m[itr][f].real() + out[f][0]/nt);
                m[itr][f].imag(add*m[itr][f].imag() + out[f][1]/nt);
            }
        }

        fftwf_destroy_plan(p);
        fftwf_free(in);
        fftwf_free(out);
        fftwf_cleanup();
    }
}

void fxTransform::testInverse() {

    std::shared_ptr<vecReg<data_t> > mod1 (new vecReg<data_t> (_domain));
    std::shared_ptr<vecReg<data_t> > mod1b (new vecReg<data_t> (_domain));

    mod1->random();

    std::shared_ptr<cvecReg<data_t> > dat1 (new cvecReg<data_t> (_range));
    dat1->zero();
    mod1b->zero();

    this->forward(false,mod1,dat1);
    this->inverse(false,mod1b,dat1);

    data_t trueNorm1, invNorm1, error1;
    trueNorm1 = mod1->norm2();
    invNorm1 = mod1b->norm2();
    mod1b->scaleAdd(mod1,1,-1);
    error1 = mod1b->norm2();

    fprintf(stderr,"input = %f output = %f inverse = %f relative error: %f\n", trueNorm1, dat1->norm2(), invNorm1, error1/trueNorm1);
}

void fxTransform::testCInverse() {

    std::shared_ptr<cvecReg<data_t> > mod1 (new cvecReg<data_t> (_domain));
    std::shared_ptr<cvecReg<data_t> > mod1b (new cvecReg<data_t> (_domain));

    mod1->random();
    
    std::shared_ptr<cvecReg<data_t> > dat1 (new cvecReg<data_t> (_crange));
    dat1->zero();
    mod1b->zero();

    this->cforward(false,mod1,dat1);
    this->cinverse(false,mod1b,dat1);
    
    data_t trueNorm1, invNorm1, error1;
    trueNorm1 = mod1->norm2();
    invNorm1 = mod1b->norm2();
    mod1b->scaleAdd(mod1,1,-1);
    error1 = mod1b->norm2();

    fprintf(stderr,"input = %f output = %f inverse = %f relative error: %f\n", trueNorm1, dat1->norm2(), invNorm1, error1/trueNorm1);
}

void fxTransform::testAdjoint(){

    std::shared_ptr<vecReg<data_t> > m (new vecReg<data_t>(_domain));
    std::shared_ptr<cvecReg<data_t> > d (new cvecReg<data_t>(_range));

    //time_t timer; unsigned int seed=time(&timer);
    m->random();
    d->random();
    
    int nf = _range.getAxis(1).n;
    int nx = _range.getN123()/nf;
    for (int i=0; i<nx; i++) d->getVals()[i*nf].imag(0);

    std::shared_ptr<vecReg<data_t> > mtild = m->clone();
    std::shared_ptr<cvecReg<data_t> > dtild = d->clone();

    this->forward(false, m, dtild);
    this->adjoint(false, mtild, d);

    data_t sum1 = dtild->real()->dot(d->real()) + dtild->imag()->dot(d->imag());
    data_t sum2 = m->dot(mtild);

    fprintf(stderr,"Dot product with add=false:\n");
    fprintf(stderr,"sum1 = %f\nsum2 = %f\n",sum1, sum2);

    this->forward(true, m, dtild);
    this->adjoint(true, mtild, d);

    sum1 = dtild->real()->dot(d->real()) + dtild->imag()->dot(d->imag());
    sum2 = m->dot(mtild);

    fprintf(stderr,"Dot product with add=true:\n");
    fprintf(stderr,"sum1 = %f\nsum2 = %f\n",sum1, sum2);
}

void fxTransform::testCAdjoint(){

    std::shared_ptr<cvecReg<data_t> > m (new cvecReg<data_t>(_domain));
    std::shared_ptr<cvecReg<data_t> > d (new cvecReg<data_t>(_crange));

    //time_t timer; unsigned int seed=time(&timer);
    m->random();
    d->random();

    std::shared_ptr<cvecReg<data_t> > mtild = m->clone();
    std::shared_ptr<cvecReg<data_t> > dtild = d->clone();

    this->cforward(false, m, dtild);
    dtild->scale(1.0/sqrt(_domain.getAxis(1).n));
    this->cadjoint(false, mtild, d);
    mtild->scale(sqrt(_domain.getAxis(1).n));

    std::complex<data_t> sum1 = dtild->dot(d);
    std::complex<data_t> sum2 = m->dot(mtild);

    fprintf(stderr,"Dot product with add=false:\n");
    fprintf(stderr,"sum1 = %f + i%f\nsum2 = %f + i%f\n",sum1.real(), sum1.imag(), sum2.real(), sum2.imag());

    m->scale(1.0/sqrt(_domain.getAxis(1).n));
    this->cforward(true, m, dtild);
    m->scale(sqrt(_domain.getAxis(1).n));
    d->scale(sqrt(_domain.getAxis(1).n));
    this->cadjoint(true, mtild, d);
    d->scale(1.0/sqrt(_domain.getAxis(1).n));

    sum1 = dtild->dot(d);
    sum2 = m->dot(mtild);

    fprintf(stderr,"Dot product with add=true:\n");
    fprintf(stderr,"sum1 = %f + i%f\nsum2 = %f + i%f\n",sum1.real(), sum1.imag(), sum2.real(), sum2.imag());
}

void fkTransform::forward(bool add, const std::shared_ptr<vecReg<data_t> > mod, std::shared_ptr<cvecReg<data_t> > dat)
{
    successCheck(checkDomainRange(mod,dat),"Vectors hypercube do not match the operator domain and range\n");

    int nt = _domain.getAxis(1).n;
    int nx = _domain.getAxis(2).n;
    int ntr = _domain.getN123()/(nt*nx);
    int nf = _range.getAxis(1).n;
    const data_t (*m) [nx][nt] = (const data_t (*) [nx][nt]) mod->getCVals();
    std::complex<data_t> (*d) [nx][nf] = (std::complex<data_t> (*) [nx][nf]) dat->getVals();

    // double precision
    if ((int)sizeof(data_t)==8){

        double * in = new double[nt*nx];
        fftw_complex * out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf*nx);
        fftw_plan p = fftw_plan_dft_r2c_2d(nx, nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, m[itr], sizeof(data_t)*nt*nx);

            // forward transform
            fftw_execute(p);

            // copy to output
            int x=0, f=0;
            for (int i=0; i<nx*nf; i++){
                x = i/nf;
                f = i - x*nf;
                x += (nx-1)/2;
                x = x - nx*(x/nx);
                d[itr][nx-x-1][f].real(add*d[itr][nx-x-1][f].real() + out[i][0]);
                d[itr][nx-x-1][f].imag(add*d[itr][nx-x-1][f].imag() + out[i][1]);
            }
        }

        delete [] in;
        fftw_destroy_plan(p);
        fftw_free(out);
        fftw_cleanup();
    }
    else
    {
        float * in = new float[nt*nx];
        fftwf_complex * out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nf*nx);
        fftwf_plan p = fftwf_plan_dft_r2c_2d(nx, nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, m[itr], sizeof(data_t)*nt*nx);

            // forward transform
            fftwf_execute(p);

            // copy to output
            int x=0, f=0;
            for (int i=0; i<nx*nf; i++){
                x = i/nf;
                f = i - x*nf;
                x += (nx-1)/2;
                x = x - nx*(x/nx);
                d[itr][nx-x-1][f].real(add*d[itr][nx-x-1][f].real() + out[i][0]);
                d[itr][nx-x-1][f].imag(add*d[itr][nx-x-1][f].imag() + out[i][1]);
            }
        }

        delete [] in;
        fftwf_destroy_plan(p);
        fftwf_free(out);
        fftwf_cleanup();
    }
}

void fkTransform::inverse(bool add, std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat)
{
    successCheck(checkDomainRange(mod,dat),"Vectors hypercube do not match the operator domain and range\n");

    int nt = _domain.getAxis(1).n;
    int nx = _domain.getAxis(2).n;
    int ntr = _domain.getN123()/(nt*nx);
    int nf = _range.getAxis(1).n;
    data_t (*m) [nx][nt] = (data_t (*) [nx][nt]) mod->getCVals();
    const std::complex<data_t> (*d) [nx][nf] = (const std::complex<data_t> (*) [nx][nf]) dat->getVals();

    // double precision
    if ((int)sizeof(data_t)==8){

        double * out = new double[nt*nx];
        fftw_complex * in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf*nx);
        fftw_plan p = fftw_plan_dft_c2r_2d(nx, nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            int x=0, f=0;
            for (int i=0; i<nx*nf; i++){
                x = i/nf;
                f = i - x*nf;
                x += (nx-1)/2;
                x = x - nx*(x/nx);
                in[i][0] = d[itr][nx-1-x][f].real();
                in[i][1] = d[itr][nx-1-x][f].imag();
            }

            // inverse transform
            fftw_execute(p);

            // copy to output
            int it=0;
            for (int i=0; i<nx*nt; i++){
                x = i/nt;
                it = i - x*nt;
                m[itr][x][it]= add*m[itr][x][it] + out[i]/(nt*nx);
            }
        }

        delete [] out;
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_cleanup();
    }
    else
    {
        float * out = new float[nt*nx];
        fftwf_complex * in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nf*nx);
        fftwf_plan p = fftwf_plan_dft_c2r_2d(nx, nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            int x=0, f=0;
            for (int i=0; i<nx*nf; i++){
                x = i/nf;
                f = i - x*nf;
                x += (nx-1)/2;
                x = x - nx*(x/nx);
                in[i][0] = d[itr][nx-1-x][f].real();
                in[i][1] = d[itr][nx-1-x][f].imag();
            }

            // inverse transform
            fftwf_execute(p);

            // copy to output
            int it=0;
            for (int i=0; i<nx*nt; i++){
                x = i/nt;
                it = i - x*nt;
                m[itr][x][it]= add*m[itr][x][it] + out[i]/(nt*nx);
            }
        }

        delete [] out;
        fftwf_destroy_plan(p);
        fftwf_free(in);
        fftwf_cleanup();
    }
}

void fkTransform::cforward(bool add, const std::shared_ptr<cvecReg<data_t> > mod, std::shared_ptr<cvecReg<data_t> > dat)
{
    successCheck(checkCDomainRange(mod,dat),"Vectors hypercube do not match the operator domain and range\n");

    int nt = _domain.getAxis(1).n;
    int nx = _domain.getAxis(2).n;
    int ntr = _domain.getN123()/(nt*nx);
    int nf = _crange.getAxis(1).n;
    const data_t (*m) [nx][nt] = (const data_t (*) [nx][nt]) mod->getCVals();
    std::complex<data_t> (*d) [nx][nf] = (std::complex<data_t> (*) [nx][nf]) dat->getVals();

    // double precision
    if ((int)sizeof(data_t)==8){

        fftw_complex * in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf*nx);
        fftw_complex * out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf*nx);
        fftw_plan p = fftw_plan_dft_2d(nx, nt, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, m[itr], sizeof(data_t)*nt*nx);

            // forward transform
            fftw_execute(p);

            // copy to output
            int x=0, f=0;
            for (int i=0; i<nx*nf; i++){
                x = i/nt;
                f = i - x*nt;
                x += (nx-1)/2;
                x = x - nx*(x/nx);
                f += (nt-1)/2;
                f = f - nt*(f/nt);
                d[itr][x][f].real(add*d[itr][x][f].real() + out[i][0]);
                d[itr][x][f].imag(add*d[itr][x][f].imag() + out[i][1]);
            }
        }

        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
        fftw_cleanup();
    }
    else
    {
        fftwf_complex * in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nf*nx);
        fftwf_complex * out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nf*nx);
        fftwf_plan p = fftwf_plan_dft_2d(nx, nt, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, m[itr], sizeof(data_t)*nt*nx);

            // forward transform
            fftwf_execute(p);

            // copy to output
            int x=0, f=0;
            for (int i=0; i<nx*nf; i++){
                x = i/nt;
                f = i - x*nt;
                x += (nx-1)/2;
                x = x - nx*(x/nx);
                f += (nt-1)/2;
                f = f - nt*(f/nt);
                d[itr][x][f].real(add*d[itr][x][f].real() + out[i][0]);
                d[itr][x][f].imag(add*d[itr][x][f].imag() + out[i][1]);
            }
        }

        fftwf_destroy_plan(p);
        fftwf_free(in);
        fftwf_free(out);
        fftwf_cleanup();
    }
}

void fkTransform::cinverse(bool add, std::shared_ptr<cvecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat)
{
    successCheck(checkCDomainRange(mod,dat),"Vectors hypercube do not match the operator domain and range\n");

    int nt = _domain.getAxis(1).n;
    int nx = _domain.getAxis(2).n;
    int ntr = _domain.getN123()/(nt*nx);
    int nf = _crange.getAxis(1).n;
    std::complex<data_t> (*m) [nx][nt] = (std::complex<data_t> (*) [nx][nt]) mod->getCVals();
    const std::complex<data_t> (*d) [nx][nf] = (const std::complex<data_t> (*) [nx][nf]) dat->getVals();

    // double precision
    if ((int)sizeof(data_t)==8){

        fftw_complex * in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nt*nx);
        fftw_complex * out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nt*nx);
        fftw_plan p = fftw_plan_dft_2d(nx, nt, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            int x=0, f=0;
            for (int i=0; i<nx*nt; i++){
                x = i/nt;
                f = i - x*nt;
                x += (nx-1)/2;
                x = x - nx*(x/nx);
                f += (nt-1)/2;
                f = f - nt*(f/nt);
                in[i][0] = d[itr][x][f].real();
                in[i][1] = d[itr][x][f].imag();
            }

            // forward transform
            fftw_execute(p);

            // copy to output
            int it=0;
            for (int i=0; i<nx*nt; i++){
                x = i/nt;
                it = i - x*nt;
                m[itr][x][f].real(add*m[itr][x][f].real() + out[i][0]/(nt*nx));
                m[itr][x][f].imag(add*m[itr][x][f].imag() + out[i][1]/(nt*nx));
            }
        }

        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
        fftw_cleanup();
    }
    else
    {
        fftwf_complex * in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt*nx);
        fftwf_complex * out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt*nx);
        fftwf_plan p = fftwf_plan_dft_2d(nx, nt, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            int x=0, f=0;
            for (int i=0; i<nx*nt; i++){
                x = i/nt;
                f = i - x*nt;
                x += (nx-1)/2;
                x = x - nx*(x/nx);
                f += (nt-1)/2;
                f = f - nt*(f/nt);
                in[i][0] = d[itr][x][f].real();
                in[i][1] = d[itr][x][f].imag();
            }

            // forward transform
            fftwf_execute(p);

            // copy to output
            int it=0;
            for (int i=0; i<nx*nt; i++){
                x = i/nt;
                it = i - x*nt;
                m[itr][x][f].real(add*m[itr][x][f].real() + out[i][0]/(nt*nx));
                m[itr][x][f].imag(add*m[itr][x][f].imag() + out[i][1]/(nt*nx));
            }
        }

        fftwf_destroy_plan(p);
        fftwf_free(in);
        fftwf_free(out);
        fftwf_cleanup();
    }
}

void gradient3d::apply_forward(bool add, const data_t * pmod, data_t * pdat){

    if (!add) memset(pdat,0,_range.getN123()*sizeof(data_t));

    int nz = _domain.getAxis(1).n;
    int nx = _domain.getAxis(2).n;
    int ny = _domain.getAxis(3).n;
    int nc = _domain.getN123() / (nz*nx*ny);
    data_t dz = 2*_domain.getAxis(1).d;
    data_t dx = 2*_domain.getAxis(2).d;
    data_t dy = 2*_domain.getAxis(3).d;

    data_t (* pm)[ny][nx][nz] = (data_t (*) [ny][nx][nz]) pmod;
    data_t (* __restrict pd)[nc][ny][nx][nz] = (data_t (*) [nc][ny][nx][nz]) pdat;

    for (int ic=0; ic<nc; ic++)
    {
        #pragma omp parallel for
        for (int iy=1; iy<ny-1; iy++){
            for (int ix=1; ix<nx-1; ix++){
                for (int iz=1; iz<nz-1; iz++){
                    pd[0][ic][iy][ix][iz] += _xw*(pm[ic][iy][ix+1][iz] - pm[ic][iy][ix-1][iz])/dx;
                    pd[1][ic][iy][ix][iz] += _yw*(pm[ic][iy+1][ix][iz] - pm[ic][iy-1][ix][iz])/dy;
                    pd[2][ic][iy][ix][iz] += _zw*(pm[ic][iy][ix][iz+1] - pm[ic][iy][ix][iz-1])/dz;
                }
            }
        }
    }
}

void gradient3d::apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

    if (!add) memset(pmod,0,_domain.getN123()*sizeof(data_t));

    int nz = _domain.getAxis(1).n;
    int nx = _domain.getAxis(2).n;
    int ny = _domain.getAxis(3).n;
    int nc = _domain.getN123() / (nz*nx*ny);
    data_t dz = 2*_domain.getAxis(1).d;
    data_t dx = 2*_domain.getAxis(2).d;
    data_t dy = 2*_domain.getAxis(3).d;

    data_t (* pm)[ny][nx][nz] = (data_t (*) [ny][nx][nz]) pmod;
    data_t (* pd)[nc][ny][nx][nz] = (data_t (*) [nc][ny][nx][nz]) pdat;

    for (int ic=0; ic<nc; ic++)
    {
        for (int iy=1; iy<ny-1; iy++){
            for (int ix=1; ix<nx-1; ix++){
                for (int iz=1; iz<nz-1; iz++){
                    pm[ic][iy][ix+1][iz] += _xw*pd[0][ic][iy][ix][iz] / dx;
                    pm[ic][iy][ix-1][iz] -= _xw*pd[0][ic][iy][ix][iz] / dx;
                    pm[ic][iy+1][ix][iz] += _yw*pd[1][ic][iy][ix][iz] / dy;
                    pm[ic][iy-1][ix][iz] -= _yw*pd[1][ic][iy][ix][iz] / dy;
                    pm[ic][iy][ix][iz+1] += _zw*pd[2][ic][iy][ix][iz] / dz;
                    pm[ic][iy][ix][iz-1] -= _zw*pd[2][ic][iy][ix][iz] / dz;
                }
            }
        }
    }
}

void laplacian3d::apply_forward(bool add, const data_t * pmod, data_t * pdat){

    if (!add) memset(pdat,0,_range.getN123()*sizeof(data_t));

    int nz = _domain.getAxis(1).n;
    int nx = _domain.getAxis(2).n;
    int ny = _domain.getAxis(3).n;
    int nc = _domain.getN123() / (nz*nx*ny);
    data_t dz = _domain.getAxis(1).d;
    data_t dx = _domain.getAxis(2).d;
    data_t dy = _domain.getAxis(3).d;
    dz *=dz;
    dx *=dx;
    dy *=dy;

    data_t (* pm)[ny][nx][nz] = (data_t (*) [ny][nx][nz]) pmod;
    data_t (* __restrict pd)[ny][nx][nz] = (data_t (*) [ny][nx][nz]) pdat;

    for (int ic=0; ic<nc; ic++)
    {
        #pragma omp parallel for
        for (int iy=1; iy<ny-1; iy++){
            for (int ix=1; ix<nx-1; ix++){
                for (int iz=1; iz<nz-1; iz++){
                    pd[ic][iy][ix][iz] += _xw*(pm[ic][iy][ix+1][iz] -2*pm[ic][iy][ix][iz] + pm[ic][iy][ix-1][iz])/dx
                                            + _yw*(pm[ic][iy+1][ix][iz] -2*pm[ic][iy][ix][iz] + pm[ic][iy-1][ix][iz])/dy
                                            + _zw*(pm[ic][iy][ix][iz+1] -2*pm[ic][iy][ix][iz] + pm[ic][iy][ix][iz-1])/dz;
                }
            }
        }
    }
}

void laplacian3d::apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

    if (!add) memset(pmod,0,_domain.getN123()*sizeof(data_t));

    int nz = _domain.getAxis(1).n;
    int nx = _domain.getAxis(2).n;
    int ny = _domain.getAxis(3).n;
    int nc = _domain.getN123() / (nz*nx*ny);
    data_t dz = _domain.getAxis(1).d;
    data_t dx = _domain.getAxis(2).d;
    data_t dy = _domain.getAxis(3).d;
    dz *=dz;
    dx *=dx;
    dy *=dy;

    data_t (* pm)[ny][nx][nz] = (data_t (*) [ny][nx][nz]) pmod;
    data_t (* __restrict pd)[ny][nx][nz] = (data_t (*) [ny][nx][nz]) pdat;

    for (int ic=0; ic<nc; ic++)
    {
        for (int iy=1; iy<ny-1; iy++){
            for (int ix=1; ix<nx-1; ix++){
                for (int iz=1; iz<nz-1; iz++){
                    pm[ic][iy][ix][iz] -=     2*pd[ic][iy][ix][iz]*(_xw/dx + _yw/dy + _zw/dz);
                    pm[ic][iy][ix+1][iz] += _xw*pd[ic][iy][ix][iz] / dx;
                    pm[ic][iy][ix-1][iz] += _xw*pd[ic][iy][ix][iz] / dx;
                    pm[ic][iy+1][ix][iz] += _yw*pd[ic][iy][ix][iz] / dy;
                    pm[ic][iy-1][ix][iz] += _yw*pd[ic][iy][ix][iz] / dy;
                    pm[ic][iy][ix][iz+1] += _zw*pd[ic][iy][ix][iz] / dz;
                    pm[ic][iy][ix][iz-1] += _zw*pd[ic][iy][ix][iz] / dz;
                }
            }
        }
    }
}

void conv1dnd::apply_forward(bool add, const data_t * pmod, data_t * pdat){

    if (!add) memset(pdat,0,_range.getN123()*sizeof(data_t));

    int nt = _domain.getAxis(1).n;
    int nx = _domain.getN123()/nt;
    int nf = _f->getN123();

    data_t * pf = _f->getVals();
    data_t scale = 1;
    //data_t scale = nt;

    int sh = 0;
    if (_centered) sh = std::floor((nf+1)/2)-1;

    #pragma omp parallel for
    for (int ix=0; ix<nx; ix++){
        for (int it=0; it<nt; it++){
            for (int j=std::max(0,it-nt+1+sh); j<std::min(nf,it+1+sh); j++) pdat[ix*nt+it] += pf[j]*pmod[ix*nt+it-j+sh]/scale;
        }
    }
}

void conv1dnd::apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

    if (!add) memset(pmod,0,_domain.getN123()*sizeof(data_t));

    int nt = _domain.getAxis(1).n;
    int nx = _domain.getN123()/nt;
    int nf = _f->getN123();

    data_t * pf = _f->getVals();
    data_t scale = 1;
    //data_t scale = nt;

    int sh = 0;
    if (_centered) sh = std::floor((nf+1)/2)-1;

    #pragma omp parallel for
    for (int ix=0; ix<nx; ix++){
        for (int it=0; it<nt; it++){
            for (int j=std::max(0,it-nt+1+sh); j<std::min(nf,it+1+sh); j++) pmod[ix*nt+it-j+sh] += pf[j]*pdat[ix*nt+it]/scale;
        }
    }
}

void convnd1d::apply_forward(bool add, const data_t * pmod, data_t * pdat){

    if (!add) memset(pdat,0,_range.getN123()*sizeof(data_t));

    int nt = _range.getAxis(1).n;
    int nx = _range.getN123()/nt;
    int nf = _domain.getN123();

    data_t * pf = _f->getVals();
    data_t scale = 1;
    //data_t scale = nt;

    int sh = 0;
    if (_centered) sh = std::floor((nf+1)/2)-1;
    
    #pragma omp parallel for
    for (int ix=0; ix<nx; ix++){
        for (int it=0; it<nt; it++){
            for (int j=std::max(0,it-nt+1+sh); j<std::min(nf,it+1+sh); j++) pdat[ix*nt+it] += pmod[j]*pf[ix*nt+it-j+sh]/scale;
        }
    }
}

void convnd1d::apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

    if (!add) memset(pmod,0,_domain.getN123()*sizeof(data_t));

    int nt = _range.getAxis(1).n;
    int nx = _range.getN123()/nt;
    int nf = _domain.getN123();

    data_t * pf = _f->getVals();
    data_t scale = 1; 
    //data_t scale = nt; 

    int sh = 0;
    if (_centered) sh = std::floor((nf+1)/2)-1;

    if (!_centered){
        #pragma omp parallel for
        for (int j=0; j<nf; j++){
            for (int ix=0; ix<nx; ix++){
                for (int it=j; it<nt; it++) pmod[j] += pdat[ix*nt+it]*pf[ix*nt+it-j]/scale;
            }
        }
    }
    else{
        for (int ix=0; ix<nx; ix++){
            for (int it=0; it<nt; it++){
                for (int j=std::max(0,it-nt+1+sh); j<std::min(nf,it+1+sh); j++) pmod[j] += pdat[ix*nt+it]*pf[ix*nt+it-j+sh]/scale;
            }
        }
    }
}





std::shared_ptr<vecReg<data_t> > zero_phase(const std::shared_ptr<vecReg<data_t> > dat){

    axis<data_t>T = dat->getHyper()->getAxis(1);
    axis<data_t>X(1,0,1);
    
    // Make a new vector with odd length
    std::shared_ptr<vecReg<data_t> > vec0;
	if (2*(T.n/2)==T.n) {
		T.n+=1;
        vec0 = std::make_shared<vecReg<data_t> >(hypercube<data_t>(T));

        for (int iz=0; iz<T.n-1; iz++){
        vec0->getVals()[iz] = dat->getVals()[iz];
        }
        vec0->getVals()[T.n-1] = 0;
    }
	
    else {
        vec0 = dat;
        vec0->setHyper((hypercube<data_t>(T,X)));
    }

    // Fourier transform the vector
    fxTransform fx(hypercube<data_t>(T,X));
    axis<data_t>F;
    F.d=fx.getRange()->getAxis(1).d;
    F.o=0;
    F.n=T.n/2 + 1;
    std::shared_ptr<cvecReg<data_t> > fxvec = std::make_shared<cvecReg<data_t> >(hypercube<data_t>(F,X));
    fx.forward(false,vec0,fxvec);

    int N = F.n;
    data_t dw = 2*M_PI/T.n;
    data_t r;
    std::complex<data_t> z (0,0);

    // Modify the FT by setting a linear phase while keeping the same amplitude spectrum
	for (int i=0; i<N; i++){

        r = abs(fxvec->getVals()[i]);

        // add linear phase shift (* exp (-j.omega.dt.(Nt-1)/2)))
        z.real(r*cos(-i*dw*(T.n - 1)/2));
        z.imag(r*sin(-i*dw*(T.n - 1)/2));

        fxvec->getVals()[i] = z;
	}

    // inverse FT and modify the original vector
    fx.inverse(false,vec0,fxvec);
    T.o = (int)(-T.n/2) * T.d; 
    vec0->setHyper(hypercube<data_t>(T));
    
    return vec0;
}
std::shared_ptr<vecReg<data_t> > minimum_phase(const std::shared_ptr<vecReg<data_t> > dat, const data_t eps){
    
    axis<data_t>T = dat->getHyper()->getAxis(1);
    axis<data_t>X(1,0,1);
    
    // Make a new vector with odd length
    std::shared_ptr<vecReg<data_t> > vec0;
	if (2*(T.n/2)==T.n) {
		T.n+=1;
        vec0 = std::make_shared<vecReg<data_t> >(hypercube<data_t>(T));

        for (int iz=0; iz<T.n-1; iz++){
        vec0->getVals()[iz] = dat->getVals()[iz];
        }
        vec0->getVals()[T.n-1] = 0;
    }
	
    else {
        vec0 = dat;
        vec0->setHyper((hypercube<data_t>(T,X)));
    }

    // Fourier transform the vector
    fxTransform fx(hypercube<data_t>(T,X));
    axis<data_t>F;
    F.d=fx.getRange()->getAxis(1).d;
    F.o=0;
    F.n=T.n/2 + 1;
    std::shared_ptr<cvecReg<data_t> > fxvec = std::make_shared<cvecReg<data_t> >(hypercube<data_t>(F,X));
    fx.forward(false,vec0,fxvec);

    int N = F.n;
    data_t r;
    std::complex<data_t> z (0,0);

    // Modify the FT by making it real and equal to the logarithm of the power spectrum
	for (int i=0; i<N; i++){

        r = sqrt(norm(fxvec->getVals()[i]));

        // set FT = log(r+eps)
        z.real(log(r+eps));
        z.imag(0);

        fxvec->getVals()[i] = z;
	}

    // inverse FT and keep the causal part
    fx.inverse(false,vec0,fxvec);
    for (int iz=1; iz<T.n/2+1; iz++){
        vec0->getVals()[iz] *= 2;
    }
    for (int iz=T.n/2+1; iz<T.n; iz++){
        vec0->getVals()[iz] = 0;
    }

    // FT again and take the exponent
    fx.forward(false,vec0,fxvec);
    for (int iz=0; iz<N; iz++){
        fxvec->getVals()[iz] = exp(fxvec->getVals()[iz]);
    }

    // inverse FT 
    fx.inverse(false,vec0,fxvec);
    T.o = 0; 
    vec0->setHyper(hypercube<data_t>(T));
    
    return vec0;
}