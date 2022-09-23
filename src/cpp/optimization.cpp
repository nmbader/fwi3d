#include "optimization.hpp"

#define ZERO 1e-16

void nlls_fwi::compute_res_and_grad(data_t * r){       

    _g->zero();

    if (_P != nullptr)  {
        _P->forward(false, _m, _p);
        _pg->zero();
    }
    else _pg=_g;

    hypercube<data_t> hyp(_p->getHyper()->getAxis(1),_p->getHyper()->getAxis(2),_p->getHyper()->getAxis(3),_p->getHyper()->getAxis(4));
    int ncxyz = hyp.getN123();

    int ns = _L->_par.ns;
    int nt = _d->getHyper()->getAxis(1).n;
    data_t dt = _d->getHyper()->getAxis(1).d;

    int size=1, rank=0;
    mpiWrapper::setSizeRank(&size,&rank);

    // cumulative number of traces
    std::vector<int> ntr_cumul(ns,0);
    ntr_cumul[0]=0;
    if (ns>1) ntr_cumul[1]=_L->_par.rxyz[0].size()*_L->_par.nrcomp;
    for (int s=2; s<ns; s++) ntr_cumul[s] = ntr_cumul[s-1] + _L->_par.rxyz[s-1].size()*_L->_par.nrcomp;

    time_t t = time(NULL);
    if (_L->_par.verbose>0 && rank==0) fprintf(stderr,"\n====================\n%s\n====================\n",ctime(&t));

    // loop over shots
    for (int s=rank; s<ns; s+=size)
    {
        if (_L->_par.verbose>1) fprintf(stderr,"Modeling shot %d by process %d\n",s, rank);

        // Build the appropriate wave equation operator
        nl_we_op_e * L;
        if (_L->_par.nmodels==2) L=new nl_we_op_a(hyp,_L->_allsrc,s,_L->_par);
        if (_L->_par.nmodels==3) L=new nl_we_op_e(hyp,_L->_allsrc,s,_L->_par);
        else if (_L->_par.nmodels==6) L=new nl_we_op_vti(hyp,_L->_allsrc,s,_L->_par);

        std::shared_ptr<vecReg<data_t> > rcv = std::make_shared<vecReg<data_t> >(*L->getRange());
        rcv->zero();
        int ntr = rcv->getN123()/nt;
        L->apply_forward(false,_p->getVals(),rcv->getVals());
        data_t * pr = rcv->getVals();

        if (_L->_par.verbose>1) fprintf(stderr,"Computing residual and adjoint sources for shot %d by process %d\n",s, rank);

        if (_w != nullptr) {
            data_t * pw = _w->getVals()+ntr_cumul[s]*nt;
            #pragma omp parallel for
            for (int i=0; i<ntr*nt; i++) pr[i] *= pw[i];
        }

        if (_filter != nullptr){
            data_t eps=1e-07;
            conv1dnd op(*rcv->getHyper(), _filter, _L->_par.filter_phase!="minimum");
            std::shared_ptr<vecReg<data_t> > output = std::make_shared<vecReg<data_t> >(*op.getRange());
            output->zero();
            op.forward(false,rcv,output);
            rcv=output;
            pr = rcv->getVals();
        }

        // rescale the synthetics by u'd/u'u (equivalent to Variable Projection with the variable being a single scaler multiplying the source time function)
        data_t * pd = _d->getVals()+ntr_cumul[s]*nt;
        if (_L->_par.scale_source_times>0){
            if (_L->_par.scale_source_times>_scale_source_times){
                data_t scaler1=0;
                data_t scaler2=0;
                // first component
                #pragma omp parallel for reduction(+: scaler1,scaler2)
                for (int i=0; i<ntr*nt; i++) {
                    scaler1 += pr[i]*pd[i];
                    scaler2 += pr[i]*pr[i];
                }
                _scalers[s] = scaler1/scaler2;
                if (std::abs(std::log(_scalers[s]))>_L->_par.scale_source_log_clip) _scalers[s]=std::exp(((_scalers[s]>=1) - (_scalers[s]<1))*_L->_par.scale_source_log_clip);
            }
            rcv->scale(_scalers[s]);
            if (_L->_par.verbose>0) fprintf(stderr,"Shot %d rescaled by a factor of %f\n",s,_scalers[s]);
        }

        std::shared_ptr<vecReg<data_t> > norms;
        std::shared_ptr<vecReg<data_t> > syn1;
        if (_L->_par.normalize){
            norms = std::make_shared<vecReg<data_t> >(hypercube<data_t>(ntr));
            syn1 = rcv->clone();
            ttnormalize(rcv->getVals(), norms->getVals(), nt, ntr);
        }

        std::shared_ptr<vecReg<data_t> > syn2;
        if (_L->_par.envelop != 0){
            syn2 = rcv->clone();
            if (_L->_par.envelop==1) envelop1(rcv);
            else envelop2(rcv);
        }

        pr = rcv->getVals();

        // compute the residual u-d
        // first component
        #pragma omp parallel for
        for (int i=0; i<ntr*nt; i++) pr[i] -= pd[i];

        memcpy(r+ntr_cumul[s]*nt, rcv->getVals(), ntr*nt*sizeof(data_t));

        // compute the gradient per shot
        applyHt(false, false, pr, pr, ntr, nt, dt, 0, ntr);

        if (_L->_par.envelop != 0){
            std::shared_ptr<vecReg<data_t> > temp = syn2->clone();
            hilbert(temp);
            std::shared_ptr<vecReg<data_t> > temp2;
            data_t * ptemp = temp->getVals();
            data_t * psyn = syn2->getVals();
            if (_L->_par.envelop==1) {
                temp2 = temp->clone();
                data_t * ptemp2 = temp2->getVals();
                for (int i=0; i<temp->getN123(); i++) {
                    ptemp2[i] = std::max((data_t)ZERO,(data_t)sqrt(ptemp[i]*ptemp[i]+psyn[i]*psyn[i]));
                    ptemp[i] *= pr[i]/ptemp2[i];
                }
            }
            else{
                for (int i=0; i<temp->getN123(); i++) ptemp[i] *= pr[i];
            }
            hilbert(temp);
            ptemp = temp->getVals();
            if (_L->_par.envelop==1) {
                data_t * ptemp2 = temp2->getVals();
                for (int i=0; i<temp->getN123(); i++) pr[i] = psyn[i]/ptemp2[i] * pr[i] - ptemp[i];
            }
            else {
                for (int i=0; i<temp->getN123(); i++) pr[i] = 2*(psyn[i]*pr[i]-ptemp[i]);
            }
        }

        if (_L->_par.normalize){
            data_t * psyn = syn1->getVals();
            data_t * pnorm = norms->getVals();
            for (int ix=0; ix<ntr; ix++){
                data_t val = 0;
                int i;
                for (int it=0; it<nt; it++){
                    i = ix*nt+it;
                    val += pr[i] * psyn[i]; // u'.r
                }
                for (int it=0; it<nt; it++){
                    i = ix*nt+it;
                    pr[i] = pr[i]/pnorm[ix] - val/(pnorm[ix]*pnorm[ix]*pnorm[ix]) * psyn[i];
                }
            }
        }

        if (_filter != nullptr){
            data_t eps=1e-07;
            conv1dnd op(*rcv->getHyper(), _filter, _L->_par.filter_phase!="minimum");
            std::shared_ptr<vecReg<data_t> > output = std::make_shared<vecReg<data_t> >(*op.getDomain());
            output->zero();
            op.adjoint(false,output,rcv);
            rcv=output;
            pr = rcv->getVals();
        }

        if (_w != nullptr) {
            data_t * pw = _w->getVals()+ntr_cumul[s]*nt;
            // first component
            #pragma omp parallel for
            for (int i=0; i<ntr*nt; i++) pr[i] *= pw[i];
        }

        rcv->scale(1.0/_dnorm);

        if (_L->_par.verbose>1) fprintf(stderr,"Computing gradient for shot %d by process %d\n",s, rank);

        L->apply_jacobianT(true,_pg->getVals(),_p->getVals(),rcv->getVals());

        delete L;

        if (_L->_par.verbose>1) fprintf(stderr,"Finish processing shot %d by process %d\n",s, rank);

    } // end of loop over shots


    // Sum all gradients
    mpiWrapper::allReduceSum(_pg->getVals(), _pg->getVals(), _pg->getN123());

    if (_P != nullptr) _P->jacobianT(false,_g,_m,_pg);

    if (_gmask != nullptr) _g->mult(_gmask);

    if (_L->_par.scale_source_times>0) _scale_source_times++;
}

#undef ZERO