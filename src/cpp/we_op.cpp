#include "we_op.hpp"
#include "misc.hpp"
#include "spatial_operators.hpp"

void analyzeGeometry(const hypercube<data_t> &model, param &par, bool verbose)
{
    successCheck(model.getNdim()==4,"The model must have four axes\n");
    par.nmodels = model.getAxis(4).n;
    successCheck(par.nmodels==2 || par.nmodels==3 || par.nmodels==6,"The model must contain 2, 3, or 6 parameters\n");

    if (par.nmodels==2){
        par.mt=false;
        par.gl=0;
        par.free_surface_stiffness=std::max(par.free_surface_stiffness,(data_t)1.0);
        par.soft_clip=0;
        par.nscomp=1;
        par.nrcomp=1;
    }
    else{
        if (par.mt) par.nscomp=6;
        else par.nscomp=3;
        if (par.gl>0) par.nrcomp=1;
        else par.nrcomp=3;
    }

    if (verbose) fprintf(stderr,"\n==========================\n Subsurface model geometry\n==========================\n");
    axis<data_t> Z = model.getAxis(1);
    axis<data_t> X = model.getAxis(2);
    axis<data_t> Y = model.getAxis(3);
    if (verbose){
        fprintf(stderr,"xmin=%.5f km, xmax=%.5f km, dx=%0.5f km, nx=%d\n",X.o,(X.n-1)*X.d+X.o,X.d,X.n);
        fprintf(stderr,"ymin=%.5f km, ymax=%.5f km, dy=%0.5f km, ny=%d\n",Y.o,(Y.n-1)*Y.d+Y.o,Y.d,Y.n);
        fprintf(stderr,"zmin=%.5f km, zmax=%.5f km, dz=%0.5f km, nz=%d\n",Z.o,(Z.n-1)*Z.d+Z.o,Z.d,Z.n);
    }

    if (verbose) fprintf(stderr,"\n==========================\n Boundary conditions\n==========================\n");
    par.bc_top = std::max(0, par.bc_top);
    par.bc_bottom = std::max(0, par.bc_bottom);
    par.bc_left = std::max(0, par.bc_left);
    par.bc_right = std::max(0, par.bc_right);
    par.bc_front = std::max(0, par.bc_front);
    par.bc_back = std::max(0, par.bc_back);
    std::string bc[5] = {"none","free surface","locally absorbing","rigid-wall","mixed rigid-absorbing"};
    if (verbose){
        fprintf(stderr,"Top boundary condition = %s\t taper size = %d\n",bc[par.bc_top].c_str(), par.taper_top); 
        fprintf(stderr,"Bottom boundary condition = %s\t taper size = %d\n",bc[par.bc_bottom].c_str(), par.taper_bottom);
        fprintf(stderr,"Left boundary condition = %s\t taper size = %d\n",bc[par.bc_left].c_str(), par.taper_left); 
        fprintf(stderr,"Right boundary condition = %s\t taper size = %d\n",bc[par.bc_right].c_str(), par.taper_right);
        fprintf(stderr,"Front boundary condition = %s\t taper size = %d\n",bc[par.bc_front].c_str(), par.taper_front); 
        fprintf(stderr,"Back boundary condition = %s\t taper size = %d\n",bc[par.bc_back].c_str(), par.taper_back);
    }

    if (verbose) {
        fprintf(stderr,"\n==========================\n Sources' and receivers' geometry\n==========================\n");
        if (par.srcoord_from_file==true) {
            fprintf(stderr,"Sources' and receivers' coordinates are read from file %s\n",par.srcoord.c_str());
        }
        else {
            fprintf(stderr,"Sources' and receivers' coordinates are read from parameters list\n");
        }
        std::string seismotype1[2] = {"displacement","velocity"};
        std::string seismotype2[2] = {"strain","strain rate"};
        std::string seismotype3[2] = {" "," derivative of "};
        if (par.nrcomp==3) fprintf(stderr,"Receivers are point measurement of type particle %s\n",seismotype1[par.seismotype].c_str());
        else{
            if (par.nmodels==2) fprintf(stderr,"Receivers are point measurement of type%shydraulic pressure\n",seismotype3[par.seismotype].c_str());
            else fprintf(stderr,"Receivers are DAS measurement of type %s with gauge length = %f km\n",seismotype2[par.seismotype].c_str(), par.gl);
        } 
    }

    if (verbose) fprintf(stderr,"Number of sources = %d\n",par.ns);
    bool check=true;
    int s=0;
    par.nr=0;
    data_t precision = std::min(X.d,Z.d)*1e-3;
    while (s<par.ns)
    {
        int nr = par.rxyz[s].size();
        par.nr+=nr;
        // check source inside model boundaries
        if (verbose) fprintf(stderr,"Shot %d located at x=%.5f km, y=%.5f km, z=%.5f km, has %d receivers\n",s,par.sxyz[s][0],par.sxyz[s][1],par.sxyz[s][2],nr);
        check = ((par.sxyz[s][0] - X.o >= -precision) && (-precision <= -par.sxyz[s][0] + X.o+(X.n-1)*X.d) && 
                (par.sxyz[s][1] - Y.o >= -precision) && (-precision <= -par.sxyz[s][1] + Y.o+(Y.n-1)*Y.d) &&
                (par.sxyz[s][2] -Z.o >= -precision) && (-precision <= -par.sxyz[s][2] + Z.o+(Z.n-1)*Z.d));
        successCheck(check,"Shot falls outside subsurface model boundaries\n");

        // check receivers inside model boundaries
        int r=0;
        while (r<nr)
        {
            data_t rx=par.rxyz[s][r][0];
            data_t ry=par.rxyz[s][r][1];
            data_t rz=par.rxyz[s][r][2];
            data_t dip=par.rxyz[s][r][3];
            data_t az=par.rxyz[s][r][4];

            if (par.gl<=0)
            {
                check = ((rx -X.o >= -precision) && (-precision <= -rx + X.o+(X.n-1)*X.d) && 
                        (ry -Y.o >= -precision) && (-precision <= -ry + Y.o+(Y.n-1)*Y.d) && 
                        (rz -Z.o >= -precision) && (-precision <= -rz + Z.o+(Z.n-1)*Z.d));
                std::string msg = "Receiver "+std::to_string(r)+" at x="+std::to_string(rx)+" km, y="+std::to_string(ry)+" km, z="+std::to_string(rz)+" km, falls outside subsurface model boundaries\n";
                successCheck(check,msg);
            }
            else
            {
                check=
                    !((rx - par.gl/2*cos(dip)*cos(az) < X.o) 
                    || (rx - par.gl/2*cos(dip)*cos(az) > X.o+(X.n-1)*X.d)
                    || (rx + par.gl/2*cos(dip)*cos(az) < X.o)
                    || (rx + par.gl/2*cos(dip)*cos(az) > X.o+(X.n-1)*X.d)
                    || (ry - par.gl/2*cos(dip)*sin(az) < Y.o) 
                    || (ry - par.gl/2*cos(dip)*sin(az) > Y.o+(Y.n-1)*Y.d)
                    || (ry + par.gl/2*cos(dip)*sin(az) < Y.o)
                    || (ry + par.gl/2*cos(dip)*sin(az) > Y.o+(Y.n-1)*Y.d)
                    || (rz - par.gl/2*sin(dip) < Z.o) 
                    || (rz - par.gl/2*sin(dip) > Z.o+(Z.n-1)*Z.d)
                    || (rz + par.gl/2*sin(dip) < Z.o)
                    || (rz + par.gl/2*sin(dip) > Z.o+(Z.n-1)*Z.d)
                    );
                std::string msg = "DAS channel "+std::to_string(r)+" at x="+std::to_string(rx)+" km, y="+std::to_string(ry)+" km, z="+std::to_string(rz)+" km, dip="+std::to_string(dip)+" rad, azimuth="+std::to_string(az)+" rad, must be away from the subsurface model boundaries by half of the gauge length along the fiber\n";
                successCheck(check,msg);
            }
            r++;
        }
        s++;
    }
    if (verbose){
        if (par.gl<=0) {
            if (par.nmodels>=3) fprintf(stderr,"Total number of 3-components receivers to be modeled = %d\n",par.nr);
            else fprintf(stderr,"Total number of hydrophones to be modeled = %d\n",par.nr);
        }
        else fprintf(stderr,"Total number of DAS channels to be modeled = %d\n",par.nr);
    }
}

std::shared_ptr<vecReg<data_t> > analyzeWavelet(std::shared_ptr<vecReg<data_t> > src, const param &par, bool verbose)
{
    if (verbose) fprintf(stderr,"\n==========================\n Source wavelet\n==========================\n");
    if (src->getHyper()->getNdim()==1)
    {
        if (par.nscomp==1)
        {
            if (verbose) {
                fprintf(stderr,"Acoustic pressure source is assumed\n");
                fprintf(stderr,"Input wavelet will be duplicated %d times\n",par.ns);
                fprintf(stderr,"Parameters fangle, fazimuth, mt, mxx, myy, mzz, mxy, mxz, myz are ignored\n");
            }
            axis<data_t> T = src->getHyper()->getAxis(1);
            axis<data_t> S(par.ns,0,1);
            hypercube<data_t> hyp(T,S);
            std::shared_ptr<vecReg<data_t> > allsrc = std::make_shared<vecReg<data_t> > (hyp);
            const data_t * p = src->getCVals();
            data_t (* pall) [T.n] = (data_t (*) [T.n]) allsrc->getVals();
            for (int its=0; its<S.n; its++) memcpy(pall[its], p, T.n*sizeof(data_t));
            return allsrc;
        }
        else if (par.nscomp==3)
        {
            if (verbose) {
                fprintf(stderr,"Vector force is assumed with a dip of %f rad clockwise wrt to the horizontal and azimuth of %f rad anti-clockwise wrt to the x-axis\n",par.fangle,par.fazimuth);
                fprintf(stderr,"Input wavelet will be duplicated 3 x %d times with the appropriate scaling\n",par.ns);
                fprintf(stderr,"Parameters mxx, myy, mzz, mxy, mxz, myz are ignored\n");
            }
            axis<data_t> T = src->getHyper()->getAxis(1);
            axis<data_t> C(3,0,1);
            axis<data_t> S(par.ns,0,1);
            hypercube<data_t> hyp(T,C,S);
            std::shared_ptr<vecReg<data_t> > allsrc = std::make_shared<vecReg<data_t> > (hyp);
            const data_t * p = src->getCVals();
            data_t (* pall) [3][T.n] = (data_t (*) [3][T.n]) allsrc->getVals();
            for (int its=0; its<S.n; its++){
                for (int it=0; it<T.n; it++){
                    pall[its][0][it] = p[it]*cos(par.fangle)*cos(par.fazimuth); // fx component
                    pall[its][1][it] = p[it]*cos(par.fangle)*sin(par.fazimuth); // fy component
                    pall[its][2][it] = p[it]*sin(par.fangle); // fz component
                }
            }
            return allsrc;
        } 
        else 
        {
            if (verbose){
                fprintf(stderr,"Moment tensor is assumed with mxx=%f, myy=%f, mzz=%f, mxy=%f, mxz=%f, myz=%f\n",par.mxx,par.myy,par.mzz,par.mxy,par.mxz,par.myz);
                fprintf(stderr,"Input wavelet will be duplicated 6 x %d times with the appropriate scaling\n",par.ns);
                fprintf(stderr,"Parameters fangle, fazimuth are ignored\n");
            }
            axis<data_t> T = src->getHyper()->getAxis(1);
            axis<data_t> C(6,0,1);
            axis<data_t> S(par.ns,0,1);
            hypercube<data_t> hyp(T,C,S);
            std::shared_ptr<vecReg<data_t> > allsrc = std::make_shared<vecReg<data_t> > (hyp);
            const data_t * p = src->getCVals();
            data_t (* pall) [6][T.n] = (data_t (*) [6][T.n]) allsrc->getVals();
            for (int its=0; its<S.n; its++){
                for (int it=0; it<T.n; it++){
                    pall[its][0][it] = p[it]*par.mxx; // mxx component
                    pall[its][1][it] = p[it]*par.myy; // myy component
                    pall[its][2][it] = p[it]*par.mzz; // mzz component
                    pall[its][3][it] = p[it]*par.mxy; // mxy component
                    pall[its][4][it] = p[it]*par.mxz; // mxz component
                    pall[its][5][it] = p[it]*par.myz; // myz component
                }
            }
            return allsrc;
        }
    }
    else
    {
        int ntr = src->getN123()/src->getHyper()->getAxis(1).n;
        if (par.nmodels==2)
        {
            std::string msg = std::to_string(par.ns) + " wavelets are expected, " + std::to_string(ntr) + " are provided\n";
            successCheck(ntr == par.ns,msg);
            if (verbose) fprintf(stderr,"Parameters fangle, fazimuth, mt, mxx, myy, mzz, mxy, mxz, myz are ignored\n");
        }
        else if (par.mt==false)
        {
            std::string msg = std::to_string(3*par.ns) + " wavelets are expected, " + std::to_string(ntr) + " are provided\n";
            successCheck(ntr == 3*par.ns,msg);
            if (verbose) fprintf(stderr,"Parameters mxx, myy, mzz, mxy, mxz, myz are ignored\n");
        }
        else
        {
            std::string msg = std::to_string(6*par.ns) + " wavelets are expected, " + std::to_string(ntr) + " are provided\n";
            successCheck(ntr == 6*par.ns,msg);
            if (verbose) fprintf(stderr,"Parameters fangle, fazimuth are ignored\n");
        }
        return src;
    }
}

void analyzeModel(const hypercube<data_t> &domain, std::shared_ptr<vecReg<data_t> > model, param &par)
{
    successCheck(model->getHyper()->getNdim()==4,"The model must have four axes\n");
    successCheck(par.nmodels==2 || par.nmodels==3 || par.nmodels==6,"The model must contain two (Vp, rho), three (Vp, Vs, rho) or six (Vp, Vs, rho, delta, eps, gamma) parameters\n");

    axis<data_t> Z = model->getHyper()->getAxis(1);
    axis<data_t> X = model->getHyper()->getAxis(2);
    axis<data_t> Y = model->getHyper()->getAxis(3);
    int nx=X.n;
    int ny=Y.n;
    int nz=Z.n;

    if (par.verbose>1) {
        fprintf(stderr,"\n==========================\n Subsurface model bounds\n==========================\n");
        if (par.nmodels==2) fprintf(stderr,"Model is assumed to contain Vp (km/s) and density Rho (g/cc) in that order\n");
        else if (par.nmodels==3) fprintf(stderr,"Model is assumed to contain Vp (km/s), Vs (km/s) and density Rho (g/cc) in that order\n");
        else if (par.nmodels==6) fprintf(stderr,"Model is assumed to contain Vp (km/s), Vs (km/s), density Rho (g/cc), Delta, Epsilon, Gamma in that order\n");
    }

    if (par.nmodels>=3)
    {
        data_t vpmin = model->min(0,ny*nx*nz);
        data_t vpmax = model->max(0,ny*nx*nz);
        data_t vsmin = model->min(ny*nx*nz,2*ny*nx*nz);
        data_t vsmax = model->max(ny*nx*nz,2*ny*nx*nz);
        data_t rhomin = model->min(2*ny*nx*nz,3*ny*nx*nz);
        data_t rhomax = model->max(2*ny*nx*nz,3*ny*nx*nz);
        data_t delmin, delmax, epsmin, epsmax, gammamin, gammamax;

        if (par.nmodels==6)
        {
            delmin = model->min(3*ny*nx*nz,4*ny*nx*nz);
            delmax = model->max(3*ny*nx*nz,4*ny*nx*nz);
            epsmin = model->min(4*ny*nx*nz,5*ny*nx*nz);
            epsmax = model->max(4*ny*nx*nz,5*ny*nx*nz);
            gammamin = model->min(5*ny*nx*nz,6*ny*nx*nz);
            gammamax = model->max(5*ny*nx*nz,6*ny*nx*nz);
        }

        if (par.verbose>1) {
            fprintf(stderr,"Vp bounds are %.2f - %.2f km/s\n",vpmin,vpmax);
            fprintf(stderr,"Vs bounds are %.2f - %.2f km/s\n",vsmin,vsmax);
            fprintf(stderr,"Rho bounds are %.2f - %.2f g/cc\n",rhomin,rhomax);
            if (par.nmodels==6){
                fprintf(stderr,"Delta bounds are %.2f - %.2f\n",delmin,delmax);
                fprintf(stderr,"Epsilon bounds are %.2f - %.2f\n",epsmin,epsmax);
                fprintf(stderr,"Gamma bounds are %.2f - %.2f\n",gammamin,gammamax);
            }
        }

        if (!par.soft_clip)
        {
            model->clip(par.vpmin,par.vpmax,0,ny*nx*nz);
            model->clip(par.vsmin,par.vsmax,ny*nx*nz,2*ny*nx*nz);
            model->clip(par.rhomin,par.rhomax,2*ny*nx*nz,3*ny*nx*nz);

            if (par.nmodels==6)
            {
                model->clip(par.deltamin,par.deltamax,3*ny*nx*nz,4*ny*nx*nz);
                model->clip(par.epsilonmin,par.epsilonmax,4*ny*nx*nz,5*ny*nx*nz);
                model->clip(par.gammamin,par.gammamax,5*ny*nx*nz,6*ny*nx*nz);
            }

            bool check=true;
            data_t (* pm) [ny][nx][nz]= (data_t (*) [ny][nx][nz]) model->getVals();
            for (int iy=0; iy<ny; iy++){
                for (int ix=0; ix<nx; ix++){
                    for (int iz=0; iz<nz; iz++){
                        if (pm[0][iy][ix][iz] < sqrt(2)*pm[1][iy][ix][iz]){
                            check = false;
                            pm[1][iy][ix][iz] = pm[0][iy][ix][iz] / sqrt(2.0001);
                        }
                    }
                }
            }
            //successCheck(check,"Vp values must always be larger than sqrt(2)*Vs\n");
            if (par.verbose>1 && !check) fprintf(stderr,"WARNING: Vs exceeds Vp/sqrt(2) at some locations and will be clipped accordingly\n");

            vpmin = model->min(0,ny*nx*nz);
            vpmax = model->max(0,ny*nx*nz);
            vsmin = model->min(ny*nx*nz,2*ny*nx*nz);
            vsmax = model->max(ny*nx*nz,2*ny*nx*nz);
            rhomin = model->min(2*ny*nx*nz,3*ny*nx*nz);
            rhomax = model->max(2*ny*nx*nz,3*ny*nx*nz);

            if (par.verbose>1) {
                fprintf(stderr,"Vp bounds after hard clipping are %.2f - %.2f km/s\n",vpmin,vpmax);
                fprintf(stderr,"Vs bounds after hard clipping are %.2f - %.2f km/s\n",vsmin,vsmax);
                fprintf(stderr,"Rho bounds after hard clipping are %.2f - %.2f g/cc\n",rhomin,rhomax);
            }
        }
        par.vmax=vpmax;
        par.vmin=vsmin;
    }

    else
    {
        data_t vpmin = model->min(0,ny*nx*nz);
        data_t vpmax = model->max(0,ny*nx*nz);
        data_t rhomin = model->min(ny*nx*nz,2*ny*nx*nz);
        data_t rhomax = model->max(ny*nx*nz,2*ny*nx*nz);

        if (par.verbose>1) {
            fprintf(stderr,"Vp bounds are %.2f - %.2f km/s\n",vpmin,vpmax);
            fprintf(stderr,"Rho bounds are %.2f - %.2f g/cc\n",rhomin,rhomax);
        }

        model->clip(par.vpmin,par.vpmax,0,ny*nx*nz);
        model->clip(par.rhomin,par.rhomax,ny*nx*nz,2*ny*nx*nz);

        vpmin = model->min(0,ny*nx*nz);
        vpmax = model->max(0,ny*nx*nz);
        rhomin = model->min(ny*nx*nz,2*ny*nx*nz);
        rhomax = model->max(ny*nx*nz,2*ny*nx*nz);

        if (par.verbose>1) {
            fprintf(stderr,"Vp bounds after hard clipping are %.2f - %.2f km/s\n",vpmin,vpmax);
            fprintf(stderr,"Rho bounds after hard clipping are %.2f - %.2f g/cc\n",rhomin,rhomax);
        }

        par.vmax=vpmax;
        par.vmin=vpmin;
    }    

    if (par.verbose>2) {
        fprintf(stderr,"\n==========================\n Dispersion analysis\n==========================\n");
        fprintf(stderr,"Maximum frequency assumed by the user = %.2f Hz\n",par.fmax);
        fprintf(stderr,"Corresponding minimum number of grid points per wavelength according to minimum velocity = %.1f\n",par.vmin/(par.fmax*std::max(std::max(X.d,Y.d),Z.d)));
    }
    
    if (par.verbose>2) fprintf(stderr,"\n==========================\n Time analysis\n==========================\n");
    axis<data_t> T = domain.getAxis(1);
    data_t tmax = (T.n-1)*T.d;
    if (par.verbose>2){
        fprintf(stderr,"Maximum duration of the input = %.3f sec, sampling = %.5f sec, number of samples = %d\n",tmax, T.d, T.n);
        fprintf(stderr,"Original Courant number = %.2f\n",par.courant);
        fprintf(stderr,"Corresponding propagation time step based on CFL condition = %.5f sec\n",par.courant*std::min(std::min(X.d,Y.d),Z.d)/par.vmax);
    }

    if (par.dt==0) par.dt=par.courant*std::min(std::min(X.d,Y.d),Z.d)/par.vmax;
    else if (par.dt > 0) {
        par.dt = std::min(par.dt,T.d);
        par.courant = par.dt * par.vmax / std::min(std::min(X.d,Y.d),Z.d);
    }
    else {
        par.dt=par.courant*std::min(std::min(X.d,Y.d),Z.d)/par.vmax;
        if (T.d > par.dt)
        {
            par.dt = T.d/ceil((T.d/par.dt));
        }
        else 
        {
            par.dt = T.d;
        }
        par.courant = par.dt * par.vmax / std::min(std::min(X.d,Y.d),Z.d);
    }
    if (par.verbose>2) fprintf(stderr,"Updated Courant number = %.2f\n",par.courant);
    par.nt = round(tmax/par.dt) + 1;
    par.tmax=(par.nt-1)*par.dt;
    if (par.verbose>2){
        fprintf(stderr,"Maximum duration of the propagation = %.3f sec, sampling = %.5f sec, number of samples = %d\n", par.tmax, par.dt, par.nt);
        if (par.resampling!="linear") par.resampling="sinc";
        fprintf(stderr,"Time resampling method = %s interpolation\n", par.resampling.c_str());
    }

    if (par.sub < 0) par.sub = round(T.d/par.dt);
    if (par.sub > 0 && par.verbose>2) fprintf(stderr,"Full wavefield (gradient) will be saved (computed) every %d propagation time steps\n", par.sub);
}

void nl_we_op_e::convert_model(data_t * m, int n, bool forward) const
{
    data_t (* pm) [n] = (data_t (*) [n]) m;

    if (forward)
    {
        for (int i=0; i<n; i++){
            pm[0][i] = pm[2][i]*(pm[0][i]*pm[0][i]-2*pm[1][i]*pm[1][i]);
            pm[1][i] = pm[2][i]*pm[1][i]*pm[1][i];
        }
    }
    else
    {
        for (int i=0; i<n; i++){
            pm[0][i] = sqrt((pm[0][i]+2*pm[1][i])/pm[2][i]);
            pm[1][i] = sqrt(pm[1][i]/pm[2][i]);
        }
    }
}

// variable coefficients expressions for 2nd derivative SBP and SAT operators
static inline data_t zero(const data_t ** par, int i){return 0;} // e.g. = zero function
static inline data_t one(const data_t ** par, int i){return 1;} // e.g. = one function
static inline data_t lam(const data_t ** par, int i){return par[0][i];} // e.g. = lambda
static inline data_t mu(const data_t ** par, int i){return par[1][i];} // e.g. = mu
static inline data_t rho(const data_t ** par, int i){return par[2][i];} // e.g. = rho
static inline data_t c13(const data_t ** par, int i){return par[3][i];} // e.g. = c13
static inline data_t eps(const data_t ** par, int i){return par[4][i];} // e.g. = eps
static inline data_t c66(const data_t ** par, int i){return par[5][i];} // e.g. = c66
static inline data_t mu2(const data_t ** par, int i){return 2*par[1][i];} // e.g. = 2 mu
static inline data_t c11b(const data_t ** par, int i){return 2*(1+2*par[4][i])*par[1][i];} // e.g. = c11b=2.(1+2.eps).mu
static inline data_t c11c(const data_t ** par, int i){return (1+2*par[4][i])*par[0][i];} // e.g. = c11c=(1+2.eps).lambda
static inline data_t c12(const data_t ** par, int i){return (1+2*par[4][i])*(par[0][i]+2*par[1][i])-2*par[5][i];} // e.g. = c12=c11-2.c66 
static inline data_t simp(const data_t ** par, int i){return 0.5*sqrt(par[1][i]*par[2][i]);} // e.g. = 1/2 sqrt(rho.mu)
static inline data_t pimp(const data_t ** par, int i){return 0.5*sqrt((par[0][i]+2*par[1][i])*par[2][i]);} // e.g. = 1/2 sqrt(rho.(lambda+2 mu))
static inline data_t pimp2(const data_t ** par, int i){return 0.5*sqrt((1+2*par[4][i])*(par[0][i]+2*par[1][i])*par[2][i]);} // e.g. = 1/2 sqrt(rho.(1+2.eps).(lambda+2 mu))
static inline data_t aimp(const data_t ** par, int i){return 0.5/sqrt(par[0][i]/par[1][i]);} // e.g. = 1/(2.sqrt(rho.K))
static inline data_t eps2(const data_t ** par, int i){return 1+2*par[4][i];} // e.g. = 1+2.eps

void nl_we_op_e::compute_gradients(const data_t * model, const data_t * u_full, const data_t * curr, const data_t * u_x,  const data_t * u_y, const data_t * u_z, data_t * tmp, data_t * grad, const param &par, int nx, int ny, int nz, int it, data_t dx, data_t dy, data_t dz, data_t dt) const
{
    // Grad_lambda for wide stencil = div(adjoint).H.Ht.div(forward)
    // Grad_mu for wide stencil = (adjointx_z+adjointz_x).H.Ht.(forwardx_z+forwardz_x)
    //      + (adjointx_y+adjointy_x).H.Ht.(forwardx_y+forwardy_x)
    //      + (adjointy_z+adjointz_y).H.Ht.(forwardy_z+forwardz_y)
    //      + 2.adjointx_x.H.Ht.forwardx_x + 2.adjointy_y.H.Ht.forwardy_y + 2.adjointz_z.H.Ht.forwardz_z
    // Grad_rho = adjointx.H.Ht.forwardx_tt + adjointy.H.Ht.forwardy_tt + adjointz.H.Ht.forwardz_tt
    // The H quadrature will be applied elsewhere to the final gradients (all shots included)

    int nxyz = nx*ny*nz;
    int itlocal = par.nt/par.sub-it;
    data_t (* __restrict pm) [nxyz] = (data_t (*)[nxyz]) model;
    data_t (* __restrict pfor) [3][nxyz] = (data_t (*) [3][nxyz]) u_full;
    data_t (* __restrict padj) [nxyz] = (data_t (*)[nxyz]) curr;
    data_t (* __restrict padj_x) [nxyz] = (data_t (*)[nxyz]) u_x;
    data_t (* __restrict padj_y) [nxyz] = (data_t (*)[nxyz]) u_y;
    data_t (* __restrict padj_z) [nxyz] = (data_t (*)[nxyz]) u_z;
    data_t (* __restrict temp) [nxyz] = (data_t (*)[nxyz]) tmp;
    data_t (* __restrict g) [nxyz] = (data_t (*)[nxyz]) grad;

    Dx(false, pfor[itlocal][0], temp[0], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz); // forwardx_x
    Dy(false, pfor[itlocal][1], temp[1], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz); // forwardy_y
    Dz(false, pfor[itlocal][2], temp[2], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz); // forwardz_z
    Dx(false, pfor[itlocal][2], temp[3], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz); // forwardz_x
    Dz(false, pfor[itlocal][0], temp[4], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz); // forwardx_z
    Dx(false, pfor[itlocal][1], temp[5], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz); // forwardy_x
    Dy(false, pfor[itlocal][0], temp[6], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz); // forwardx_y
    Dy(false, pfor[itlocal][2], temp[7], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz); // forwardz_y
    Dz(false, pfor[itlocal][1], temp[8], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz); // forwardy_z

    // different from time zero
    if (itlocal>0){
        #pragma omp parallel for
        for (int i=0; i<nxyz; i++) 
        {
            g[0][i] += dt*(padj_x[0][i] + padj_y[1][i] + padj_z[2][i])*(temp[0][i] + temp[1][i] + temp[2][i]); // lambda gradient
            g[1][i] += dt*((padj_z[0][i] + padj_x[2][i])*(temp[3][i] + temp[4][i])
                            + (padj_y[0][i] + padj_x[1][i])*(temp[5][i] + temp[6][i])
                            + (padj_z[1][i] + padj_y[2][i])*(temp[7][i] + temp[8][i])
                            + 2*padj_x[0][i]*temp[0][i] + 2*padj_y[1][i]*temp[1][i] + 2*padj_z[2][i]*temp[2][i]); // mu gradient
            g[2][i] += 1.0/dt*(padj[0][i]*(pfor[itlocal+1][0][i]-2*pfor[itlocal][0][i]+pfor[itlocal-1][0][i]) 
                            + padj[1][i]*(pfor[itlocal+1][1][i]-2*pfor[itlocal][1][i]+pfor[itlocal-1][1][i])
                            + padj[2][i]*(pfor[itlocal+1][2][i]-2*pfor[itlocal][2][i]+pfor[itlocal-1][2][i]) ); // rho gradient
        }
    }
    // time zero
    else{
        #pragma omp parallel for
        for (int i=0; i<nxyz; i++) 
        {
            g[0][i] += 0.5*dt*(padj_x[0][i] + padj_y[1][i] + padj_z[2][i])*(temp[0][i] + temp[1][i] + temp[2][i]); // lambda gradient
            g[1][i] += 0.5*dt*((padj_z[0][i] + padj_x[2][i])*(temp[3][i] + temp[4][i])
                            + (padj_y[0][i] + padj_x[1][i])*(temp[5][i] + temp[6][i])
                            + (padj_z[1][i] + padj_y[2][i])*(temp[7][i] + temp[8][i])
                            + 2*padj_x[0][i]*temp[0][i] + 2*padj_y[1][i]*temp[1][i] + 2*padj_z[2][i]*temp[2][i]); // mu gradient
            g[2][i] += 1.0/dt*( padj[0][i]*(pfor[itlocal+1][0][i]-pfor[itlocal][0][i]) 
                                + padj[1][i]*(pfor[itlocal+1][1][i]-pfor[itlocal][1][i])
                                + padj[2][i]*(pfor[itlocal+1][2][i]-pfor[itlocal][2][i]) ); // rho gradient
        }
    }
}

void nl_we_op_e::propagate(bool adj, const data_t * model, const data_t * src, data_t * rcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int ny, int nz, data_t dx, data_t dy, data_t dz, data_t ox, data_t oy, data_t oz) const
{
    // wavefields allocations and pointers
    int nxyz = nx*ny*nz;
    data_t * u  = new data_t[9*nxyz];
    data_t * dux  = new data_t[3*nxyz];
    data_t * duy  = new data_t[3*nxyz];
    data_t * duz  = new data_t[3*nxyz];
    memset(u, 0, 9*nxyz*sizeof(data_t));
    memset(dux, 0, 3*nxyz*sizeof(data_t));
    memset(duy, 0, 3*nxyz*sizeof(data_t));
    memset(duz, 0, 3*nxyz*sizeof(data_t));
    data_t (* __restrict prev) [nxyz] = (data_t (*)[nxyz]) u;
    data_t (* __restrict curr) [nxyz] = (data_t (*)[nxyz]) (u + 3*nxyz);
    data_t (* __restrict next) [nxyz] = (data_t (*)[nxyz]) (u + 6*nxyz);
    data_t (* __restrict bucket) [nxyz];
    data_t (* __restrict u_x) [nxyz] = (data_t (*)[nxyz]) dux;
    data_t (* __restrict u_y) [nxyz] = (data_t (*)[nxyz]) duy;
    data_t (* __restrict u_z) [nxyz] = (data_t (*)[nxyz]) duz;
    data_t (* __restrict u_full) [3][nxyz];
    data_t * __restrict tmp;

    if (par.sub>0) u_full = (data_t (*) [3][nxyz]) full_wfld;
    if (grad != nullptr) 
    {
        tmp = new data_t[9*nxyz];
        memset(tmp, 0, 9*nxyz*sizeof(data_t));
    }
    else {
        tmp = new data_t[nxyz];
        memset(tmp, 0, nxyz*sizeof(data_t));
    }
	const data_t* mod[3] = {model, model+nxyz, model+2*nxyz};
    
    
    // source and receiver components for injection/extraction
    int nscomp = par.nscomp;
    int nrcomp = par.nrcomp;
    if (adj) {nscomp = par.nrcomp; nrcomp=par.nscomp;}

    // prev = 1/(2 * rho) * dt2 * src
    inj->inject(false, src, prev[0], nx, ny, nz, par.nt, inj->_npts, 0, 0, 0, inj->_npts, inj->_xind.data(), inj->_yind.data(), inj->_zind.data(), inj->_xw.data(), inj->_yw.data(), inj->_zw.data());
    inj->inject(false, src, prev[1], nx, ny, nz, par.nt, inj->_npts, 1, 0, 0, inj->_npts, inj->_xind.data(), inj->_yind.data(), inj->_zind.data(), inj->_xw.data(), inj->_yw.data(), inj->_zw.data());
    inj->inject(false, src, prev[2], nx, ny, nz, par.nt, inj->_npts, 2, 0, 0, inj->_npts, inj->_xind.data(), inj->_yind.data(), inj->_zind.data(), inj->_xw.data(), inj->_yw.data(), inj->_zw.data());
    for (int i=0; i<nxyz; i++)
    {
        prev[0][i]  *= 0.5 * par.dt*par.dt/ mod[2][i];
        prev[1][i]  *= 0.5 * par.dt*par.dt/ mod[2][i];
        prev[2][i]  *= 0.5 * par.dt*par.dt/ mod[2][i];
    }    

    int pct10 = round(par.nt/10);

    for (int it=0; it<par.nt-1; it++)
    {
        // copy the current wfld to the full wfld vector
        if ((par.sub>0) && (it%par.sub==0) && (grad==nullptr)) memcpy(u_full[it/par.sub], curr, 3*nxyz*sizeof(data_t));

        // extract receivers
        if (grad == nullptr)
        {
            ext->extract(true, curr[0], rcv, nx, ny, nz, par.nt, ext->_npts, 0, it, 0, ext->_npts, ext->_xind.data(), ext->_yind.data(), ext->_zind.data(), ext->_xw.data(), ext->_yw.data(), ext->_zw.data());
            ext->extract(true, curr[1], rcv, nx, ny, nz, par.nt, ext->_npts, 1, it, 0, ext->_npts, ext->_xind.data(), ext->_yind.data(), ext->_zind.data(), ext->_xw.data(), ext->_yw.data(), ext->_zw.data());
            ext->extract(true, curr[2], rcv, nx, ny, nz, par.nt, ext->_npts, 2, it, 0, ext->_npts, ext->_xind.data(), ext->_yind.data(), ext->_zind.data(), ext->_xw.data(), ext->_yw.data(), ext->_zw.data());
        }

        // compute FWI gradients except for first and last time samples
        if ((grad != nullptr) && (par.sub>0) && (it%par.sub==0) && it!=0) compute_gradients(model, full_wfld, curr[0], u_x[0], u_y[0], u_z[0], tmp, grad, par, nx, ny, nz, it/par.sub, dx, dy, dz, par.sub*par.dt);

        // apply spatial SBP operators
        Dxx_var<mu>(false, curr[0], next[0], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz, mod, 2.0);
        Dyy_var<mu>(true,  curr[0], next[0], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz, mod, 1.0);
        Dzz_var<mu>(true,  curr[0], next[0], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz, mod, 1.0);

        Dxx_var<mu>(false, curr[1], next[1], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz, mod, 1.0);
        Dyy_var<mu>(true,  curr[1], next[1], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz, mod, 2.0);
        Dzz_var<mu>(true,  curr[1], next[1], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz, mod, 1.0);
        
        Dxx_var<mu>(false, curr[2], next[2], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz, mod, 1.0);
        Dyy_var<mu>(true,  curr[2], next[2], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz, mod, 1.0);
        Dzz_var<mu>(true,  curr[2], next[2], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz, mod, 2.0);

        
       
        Dx(false, curr[0], u_x[0], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz);
        Dy(false, curr[1], u_y[1], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz);
        Dz(false, curr[2], u_z[2], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz);
        #pragma omp parallel for
        for (int i=0; i<nxyz; i++)
        {
            tmp[i] = mod[0][i]*(u_x[0][i] + u_y[1][i] + u_z[2][i]);
        }
        Dx(true, tmp, next[0], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz);
        Dy(true, tmp, next[1], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz);
        Dz(true, tmp, next[2], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz);
    


        Dx(false, curr[1], u_x[1], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz);
        Dx(false, curr[2], u_x[2], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz);     
        mult_Dy(true, u_x[1], next[0], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz, mod[1]);
        mult_Dz(true, u_x[2], next[0], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz, mod[1]);

        Dy(false, curr[0], u_y[0], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz);
        Dy(false, curr[2], u_y[2], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz);     
        mult_Dx(true, u_y[0], next[1], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz, mod[1]);
        mult_Dz(true, u_y[2], next[1], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz, mod[1]);

        Dz(false, curr[0], u_z[0], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz);
        Dz(false, curr[1], u_z[1], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz);     
        mult_Dx(true, u_z[0], next[2], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz, mod[1]);
        mult_Dy(true, u_z[1], next[2], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz, mod[1]);

        // inject sources
        inj->inject(true, src, next[0], nx, ny, nz, par.nt, inj->_npts, 0, it, 0, inj->_npts, inj->_xind.data(), inj->_yind.data(), inj->_zind.data(), inj->_xw.data(), inj->_yw.data(), inj->_zw.data());
        inj->inject(true, src, next[1], nx, ny, nz, par.nt, inj->_npts, 1, it, 0, inj->_npts, inj->_xind.data(), inj->_yind.data(), inj->_zind.data(), inj->_xw.data(), inj->_yw.data(), inj->_zw.data());
        inj->inject(true, src, next[2], nx, ny, nz, par.nt, inj->_npts, 2, it, 0, inj->_npts, inj->_xind.data(), inj->_yind.data(), inj->_zind.data(), inj->_xw.data(), inj->_yw.data(), inj->_zw.data());


        // apply boundary conditions
        if (par.bc_top==1)
        {
            const data_t * in[4] = {curr[2],curr[1],curr[0],prev[0]};
            esat_neumann_top<mu,zero,zero,mu,zero>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
            in[2] = curr[1]; in[1] = curr[2];
            esat_neumann_top<zero,mu,zero,mu,zero>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
            in[0] = curr[0]; in[1] = curr[1]; in[2] = curr[2];
            esat_neumann_top<lam,lam,lam,mu2,zero>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
        }
        else if (par.bc_top==2)
        {
            const data_t * in[4] = {curr[2],curr[1],curr[0],prev[0]};
            esat_neumann_top<mu,zero,zero,mu,simp>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
            in[2] = curr[1]; in[1] = curr[2]; in[3]=prev[1];
            esat_neumann_top<zero,mu,zero,mu,simp>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
            in[0] = curr[0]; in[1] = curr[1]; in[2] = curr[2]; in[3]=prev[2];
            esat_neumann_top<lam,lam,lam,mu2,pimp>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
        }
        if (par.bc_bottom==1)
        {
            const data_t * in[4] = {curr[2],curr[1],curr[0],prev[0]};
            esat_neumann_bottom<mu,zero,zero,mu,zero>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
            in[2] = curr[1]; in[1] = curr[2];
            esat_neumann_bottom<zero,mu,zero,mu,zero>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
            in[0] = curr[0]; in[1] = curr[1]; in[2] = curr[2];
            esat_neumann_bottom<lam,lam,lam,mu2,zero>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
        }
        else if (par.bc_bottom==2)
        {
            const data_t * in[4] = {curr[2],curr[1],curr[0],prev[0]};
            esat_neumann_bottom<mu,zero,zero,mu,simp>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
            in[2] = curr[1]; in[1] = curr[2]; in[3]=prev[1];
            esat_neumann_bottom<zero,mu,zero,mu,simp>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
            in[0] = curr[0]; in[1] = curr[1]; in[2] = curr[2]; in[3]=prev[2];
            esat_neumann_bottom<lam,lam,lam,mu2,pimp>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
        }
        if (par.bc_left==1)
        {
            const data_t * in[4] = {curr[0],curr[2],curr[1],prev[0]};
            esat_neumann_left<mu,zero,zero,mu,zero>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
            in[2] = curr[2]; in[1] = curr[0];
            esat_neumann_left<zero,mu,zero,mu,zero>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
            in[0] = curr[1]; in[1] = curr[2]; in[2] = curr[0];
            esat_neumann_left<lam,lam,lam,mu2,zero>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
        }
        else if (par.bc_left==2)
        {
            const data_t * in[4] = {curr[0],curr[2],curr[1],prev[1]};
            esat_neumann_left<mu,zero,zero,mu,simp>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
            in[2] = curr[2]; in[1] = curr[0]; in[3]=prev[2];
            esat_neumann_left<zero,mu,zero,mu,simp>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
            in[0] = curr[1]; in[1] = curr[2]; in[2] = curr[0]; in[3]=prev[0];
            esat_neumann_left<lam,lam,lam,mu2,pimp>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
        }
        if (par.bc_right==1)
        {
            const data_t * in[4] = {curr[0],curr[2],curr[1],prev[0]};
            esat_neumann_right<mu,zero,zero,mu,zero>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
            in[2] = curr[2]; in[1] = curr[0];
            esat_neumann_right<zero,mu,zero,mu,zero>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
            in[0] = curr[1]; in[1] = curr[2]; in[2] = curr[0];
            esat_neumann_right<lam,lam,lam,mu2,zero>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
        }
        else if (par.bc_right==2)
        {
            const data_t * in[4] = {curr[0],curr[2],curr[1],prev[1]};
            esat_neumann_right<mu,zero,zero,mu,simp>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
            in[2] = curr[2]; in[1] = curr[0]; in[3]=prev[2];
            esat_neumann_right<zero,mu,zero,mu,simp>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
            in[0] = curr[1]; in[1] = curr[2]; in[2] = curr[0]; in[3]=prev[0];
            esat_neumann_right<lam,lam,lam,mu2,pimp>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
        }
        if (par.bc_front==1)
        {
            const data_t * in[4] = {curr[1],curr[2],curr[0],prev[0]};
            esat_neumann_front<mu,zero,zero,mu,zero>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
            in[2] = curr[2]; in[1] = curr[1];
            esat_neumann_front<zero,mu,zero,mu,zero>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
            in[0] = curr[0]; in[1] = curr[2]; in[2] = curr[1];
            esat_neumann_front<lam,lam,lam,mu2,zero>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
        }
        else if (par.bc_front==2)
        {
            const data_t * in[4] = {curr[1],curr[2],curr[0],prev[0]};
            esat_neumann_front<mu,zero,zero,mu,simp>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
            in[2] = curr[2]; in[1] = curr[1]; in[3]=prev[2];
            esat_neumann_front<zero,mu,zero,mu,simp>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
            in[0] = curr[0]; in[1] = curr[2]; in[2] = curr[1]; in[3]=prev[1];
            esat_neumann_front<lam,lam,lam,mu2,pimp>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
        }
        if (par.bc_back==1)
        {
            const data_t * in[4] = {curr[1],curr[2],curr[0],prev[0]};
            esat_neumann_back<mu,zero,zero,mu,zero>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
            in[2] = curr[2]; in[1] = curr[1];
            esat_neumann_back<zero,mu,zero,mu,zero>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
            in[0] = curr[0]; in[1] = curr[2]; in[2] = curr[1];
            esat_neumann_back<lam,lam,lam,mu2,zero>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
        }
        else if (par.bc_back==2)
        {
            const data_t * in[4] = {curr[1],curr[2],curr[0],prev[0]};
            esat_neumann_back<mu,zero,zero,mu,simp>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
            in[2] = curr[2]; in[1] = curr[1]; in[3]=prev[2];
            esat_neumann_back<zero,mu,zero,mu,simp>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
            in[0] = curr[0]; in[1] = curr[2]; in[2] = curr[1]; in[3]=prev[1];
            esat_neumann_back<lam,lam,lam,mu2,pimp>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
        }

        // update wfld with the 2 steps time recursion
        #pragma omp parallel for
        for (int i=0; i<nxyz; i++)
        {
            next[0][i] = par.dt*par.dt*next[0][i]/mod[2][i] + 2*curr[0][i] - prev[0][i];
            next[1][i] = par.dt*par.dt*next[1][i]/mod[2][i] + 2*curr[1][i] - prev[1][i];
            next[2][i] = par.dt*par.dt*next[2][i]/mod[2][i] + 2*curr[2][i] - prev[2][i];
        }

        // scale boundaries when relevant (for locally absorbing BC only)
        data_t * in[3] = {next[0],next[1],next[2]};
        esat_scale_boundaries<one>(in, nx, ny, nz, dx, dy, dz, 0, nx, 0, ny, 0, nz, mod, par.dt, par.bc_top==2, par.bc_bottom==2, par.bc_left==2, par.bc_right==2, par.bc_front==2, par.bc_back==2);
        
        // apply taper 
        taperz(curr[0], nx, ny, nz, 3, 0, nx, 0, ny, par.taper_top, 0, 0, 3, par.taper_strength);
        taperz(next[0], nx, ny, nz, 3, 0, nx, 0, ny, par.taper_top, 0, 0, 3, par.taper_strength);
        taperz(curr[0], nx, ny, nz, 3, 0, nx, 0, ny, nz-par.taper_bottom, nz, 0, 3, par.taper_strength);
        taperz(next[0], nx, ny, nz, 3, 0, nx, 0, ny, nz-par.taper_bottom, nz, 0, 3, par.taper_strength);
        taperx(curr[0], nx, ny, nz, 3, par.taper_left, 0, 0, ny, 0, nz, 0, 3, par.taper_strength);
        taperx(next[0], nx, ny, nz, 3, par.taper_left, 0, 0, ny, 0, nz, 0, 3, par.taper_strength);
        taperx(curr[0], nx, ny, nz, 3, nx-par.taper_right, nx, 0, ny, 0, nz, 0, 3, par.taper_strength);
        taperx(next[0], nx, ny, nz, 3, nx-par.taper_right, nx, 0, ny, 0, nz, 0, 3, par.taper_strength);
        tapery(curr[0], nx, ny, nz, 3, 0, nx, par.taper_front, 0, 0, nz, 0, 3, par.taper_strength);
        tapery(next[0], nx, ny, nz, 3, 0, nx, par.taper_front, 0, 0, nz, 0, 3, par.taper_strength);
        tapery(curr[0], nx, ny, nz, 3, 0, nx, ny-par.taper_back, ny, 0, nz, 0, 3, par.taper_strength);
        tapery(next[0], nx, ny, nz, 3, 0, nx, ny-par.taper_back, ny, 0, nz, 0, 3, par.taper_strength);
       
        bucket=prev;
        prev=curr;
        curr=next;
        next=bucket;

        if ((it+1) % pct10 == 0 && par.verbose>2) fprintf(stderr,"Propagation progress = %d\%\n",10*(it+1)/pct10);
    }

    // copy the last wfld to the full wfld vector
    if ((par.sub>0) && ((par.nt-1)%par.sub==0) && (grad==nullptr)) memcpy(u_full[(par.nt-1)/par.sub], curr, 3*nxyz*sizeof(data_t));

    // extract receivers last sample
    if (grad == nullptr)
    {
        ext->extract(true, curr[0], rcv, nx, ny, nz, par.nt, ext->_npts, 0, par.nt-1, 0, ext->_npts, ext->_xind.data(), ext->_yind.data(), ext->_zind.data(), ext->_xw.data(), ext->_yw.data(), ext->_zw.data());
        ext->extract(true, curr[1], rcv, nx, ny, nz, par.nt, ext->_npts, 1, par.nt-1, 0, ext->_npts, ext->_xind.data(), ext->_yind.data(), ext->_zind.data(), ext->_xw.data(), ext->_yw.data(), ext->_zw.data());
        ext->extract(true, curr[2], rcv, nx, ny, nz, par.nt, ext->_npts, 2, par.nt-1, 0, ext->_npts, ext->_xind.data(), ext->_yind.data(), ext->_zind.data(), ext->_xw.data(), ext->_yw.data(), ext->_zw.data());
    }

    // last sample gradient
    if ((grad != nullptr) && (par.sub>0) && ((par.nt-1)%par.sub==0) ) compute_gradients(model, full_wfld, curr[0], u_x[0], u_y[0], u_z[0], tmp, grad, par, nx, ny, nz, (par.nt-1)/par.sub, dx, dy, dz, par.sub*par.dt);
    
    delete [] u;
    delete [] dux;
    delete [] duy;
    delete [] duz;
    delete [] tmp;
}

void nl_we_op_e::apply_forward(bool add, const data_t * pmod, data_t * pdat)
{
    std::shared_ptr<vecReg<data_t> > model = std::make_shared<vecReg<data_t> >(_domain);
    memcpy(model->getVals(), pmod, _domain.getN123()*sizeof(data_t));
    analyzeModel(*_src->getHyper(),model,_par);
    convert_model(model->getVals(), model->getN123()/_par.nmodels, true);
    const data_t * pm = model->getCVals();

    axis<data_t> X = _domain.getAxis(2);
    axis<data_t> Y = _domain.getAxis(3);
    axis<data_t> Z = _domain.getAxis(1);
    hypercube<data_t> domain = *_src->getHyper();

    if (_par.sub>0){
        if (_par.nmodels>=3) _full_wfld=std::make_shared<vecReg<data_t> > (hypercube<data_t>(Z,X,Y,axis<data_t>(3,0,1), axis<data_t>(1+_par.nt/_par.sub,0,_par.dt*_par.sub)));
        else _full_wfld=std::make_shared<vecReg<data_t> > (hypercube<data_t>(Z,X,Y, axis<data_t>(1+_par.nt/_par.sub,0,_par.dt*_par.sub)));
    }

    // setting up and apply the time resampling operator
    resampler * resamp;

    axis<data_t> T = domain.getAxis(1);
    axis<data_t> Xr = _range.getAxis(2);
    axis<data_t> Cr0 (_par.nrcomp,0,1);
    axis<data_t> Cr (3,0,1);
    axis<data_t> Cs(_par.nscomp,0,1);
    if (_par.nmodels==2) Cr.n = 1;
    data_t alpha = _par.dt/T.d; // ratio between sampling rates
    T.n = _par.nt;
    T.d = _par.dt;
    hypercube<data_t> hyper_s(T,Cs);
    hypercube<data_t> hyper_r0(T,Xr,Cr0);
    hypercube<data_t> hyper_r(T,Xr,Cr);
    std::shared_ptr<vecReg<data_t> > src = std::make_shared<vecReg<data_t> >(hyper_s);
    std::shared_ptr<vecReg<data_t> > rcv = std::make_shared<vecReg<data_t> >(hyper_r);
    src->zero();
    rcv->zero();

    if (_par.resampling == "linear") resamp = new linear_resampler(domain, hyper_s);
    else resamp = new sinc_resampler(domain, hyper_s, _par.sinc_half_length);

    // resample in time (interpolation)
    resamp->apply_forward(false,_src->getCVals(),src->getVals());

    // setting up the injection and extraction operators
    injector * inj;
    if (_par.mt==true) inj = new ddelta_m3(_domain,{_par.sxyz[0]});
    else inj = new delta_m3(_domain,{_par.sxyz[0]});

    injector * ext;
    if (_par.gl<=0) ext = new delta_m3(_domain,_par.rxyz[0]);
    else ext = new dipole_m3(_domain,_par.rxyz[0],_par.gl);

    // perform the wave propagation
    data_t * full = nullptr;
    if (_par.sub>0) full = _full_wfld->getVals();
    propagate(false, pm, src->getCVals(), rcv->getVals(), inj, ext, full, nullptr, _par, X.n, Y.n, Z.n, X.d, Y.d, Z.d, X.o, Y.o, Z.o);

    delete inj;
    delete ext;

    // convert to DAS (strain) data when relevant
    if (_par.gl>0) dipole_to_strain(false, rcv->getVals(), _par.rdip.data(), _par.raz.data(), Xr.n, T.n, 0, Xr.n, _par.gl);

    // convert to particle velocity or strain rate when relevant: forward  mode
    if (_par.seismotype==1)
    {
        std::shared_ptr<vecReg<data_t> > rcv_dt = std::make_shared<vecReg<data_t> >(hyper_r0);
        rcv_dt->zero();
        Dt(false, false, rcv->getCVals(), rcv_dt->getVals(), Xr.n*Cr0.n, T.n, T.d, 0, Xr.n*Cr0.n);
        rcv = rcv_dt;
    }

    // resample in time (decimation)
    rcv->scale(alpha);
    resamp->setDomainRange(_range, hyper_r0);
    resamp->apply_adjoint(add, pdat, rcv->getCVals());
        
    delete resamp;
}

void nl_we_op_e::apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat)
{

    std::shared_ptr<vecReg<data_t> > model0 = std::make_shared<vecReg<data_t> >(_domain);
    memcpy(model0->getVals(), pmod0, _domain.getN123()*sizeof(data_t));
    int verbose = _par.verbose;
    _par.verbose=0;
    analyzeModel(*_src->getHyper(),model0,_par);
    _par.verbose=verbose;
    convert_model(model0->getVals(), model0->getN123()/_par.nmodels, true);
    const data_t * pm0 = model0->getCVals();

    axis<data_t> X = _domain.getAxis(2);
    axis<data_t> Y = _domain.getAxis(3);
    axis<data_t> Z = _domain.getAxis(1);
    hypercube<data_t> domain = *_src->getHyper();
    int nxyz = X.n*Y.n*Z.n;

    data_t * temp;
    int nm=std::min(3,_par.nmodels);
    if (!add) memset(pmod, 0, nm*nxyz*sizeof(data_t));
    else { // copy pre-existant gradient to a temporary container
        temp = new data_t[nm*nxyz];
        memcpy(temp, pmod, nm*nxyz*sizeof(data_t));
    }

    // setting up and apply the time resampling operator
    resampler * resamp;

    axis<data_t> T = domain.getAxis(1);
    axis<data_t> Xr = _range.getAxis(2);
    axis<data_t> Cr0 (_par.nrcomp,0,1);
    axis<data_t> Cr (3,0,1);
    axis<data_t> Cs(_par.nscomp,0,1);
    if (_par.nmodels==2) Cr.n = 1;
    data_t alpha = _par.dt/T.d; // ratio between sampling rates
    T.n = _par.nt;
    T.d = _par.dt;
    hypercube<data_t> hyper_s(T,Cs);
    hypercube<data_t> hyper_r0(T,Xr,Cr0);
    hypercube<data_t> hyper_r(T,Xr,Cr);
    std::shared_ptr<vecReg<data_t> > src = std::make_shared<vecReg<data_t> >(hyper_s);
    std::shared_ptr<vecReg<data_t> > rcv = std::make_shared<vecReg<data_t> >(hyper_r);
    src->zero();
    rcv->zero();
    
    if (_par.resampling == "linear") resamp = new linear_resampler(_range, hyper_r0);
    else resamp = new sinc_resampler(_range, hyper_r0, _par.sinc_half_length);

    // resample in time (interpolation)
    resamp->apply_forward(false,pdat,rcv->getVals());

    // apply inverse time quadrature, revert in time and multiply by -1
    applyHt(true, false, rcv->getCVals(), rcv->getVals(), Xr.n*Cr.n, T.n, T.d, 0, Xr.n*Cr.n);
    rcv->revert(1);
    rcv->scale(-alpha);

    // convert from particle velocity or strain rate when relevant: adjoint mode
    if (_par.seismotype==1)
    {
        std::shared_ptr<vecReg<data_t> > rcv_dt = std::make_shared<vecReg<data_t> >(hyper_r);
        rcv_dt->zero();
        Dt(true, false, rcv->getCVals(), rcv_dt->getVals(), Xr.n*Cr.n, T.n, T.d, 0, Xr.n*Cr.n);
        rcv = rcv_dt;
    }

    // convert from DAS (strain) to dipole data when relevant
    if (_par.gl>0) dipole_to_strain(true, rcv->getVals(), _par.rdip.data(), _par.raz.data(), Xr.n, T.n, 0, Xr.n, _par.gl);

    // setting up the injection and extraction operators
    injector * ext;
    if (_par.mt==true) ext = new ddelta_m3(_domain,{_par.sxyz[0]});
    else ext = new delta_m3(_domain,{_par.sxyz[0]});

    injector * inj;
    if (_par.gl<=0) inj = new delta_m3(_domain,_par.rxyz[0]);
    else inj = new dipole_m3(_domain,_par.rxyz[0],_par.gl);

    // perform the wave propagation
    data_t * full = nullptr;
    if (_par.sub>0) full = _full_wfld->getVals();
    propagate(true, pm0, rcv->getCVals(), src->getVals(), inj, ext, full, pmod, _par, X.n, Y.n, Z.n, X.d, Y.d, Z.d, X.o, Y.o, Z.o);

    delete inj;
    delete ext;

    // apply H quadrature to final gradient
    for (int p=0; p<nm; p++)
    {
        applyHx(false, false, pmod+p*nxyz, pmod+p*nxyz, X.n, Y.n, Z.n, X.d, 0, X.n, 0, Y.n, 0, Z.n);
        applyHy(false, false, pmod+p*nxyz, pmod+p*nxyz, X.n, Y.n, Z.n, Y.d, 0, X.n, 0, Y.n, 0, Z.n);
        applyHz(false, false, pmod+p*nxyz, pmod+p*nxyz, X.n, Y.n, Z.n, Z.d, 0, X.n, 0, Y.n, 0, Z.n);
    }

    // convert gradients from lambda, mu, rho to vp, vs, rho    or   from K, rho-1 to vp, rho
    // Grad_vp = 2.rho.vp.Grad_lambda = 2.sqrt(rho(lambda+2mu)).Grad_lambda
    // Grad_vs = 2.rho.vs.(Grad_mu - 2.Grad_lambda) = 2.sqrt(rho.mu).(Grad_mu - 2.Grad_lambda)
    // Grad_rho = (vp2 - 2.vs2).Grad_lambda + vs2.Grad_mu + Grad_rho = (lambda/rho).Grad_lambda + (mu/rho).Grad_mu + Grad_rho
    // or
    // Grad_vp = 2.rho.vp.Grad_K = 2.sqrt(K/rho-1).Grad_K
    // Grad_rho = vp^2.Grad_K - 1/rho^2.Grad_rho-1 = K.rho-1.Grad_K - (rho-1)^2.Grad_rho-1
    data_t (*g) [nxyz] = (data_t (*)[nxyz]) pmod;
    data_t (*pm) [nxyz] = (data_t (*)[nxyz]) pm0;
    if (nm>=3)
    {
        #pragma omp parallel for
        for(int i=0; i<nxyz; i++)
        {
            g[2][i] = (pm[0][i]/pm[2][i])*g[0][i] + (pm[1][i]/pm[2][i])*g[1][i] + g[2][i]; // Rho Gradient
            g[1][i] = 2*sqrt(pm[1][i]*pm[2][i])*(g[1][i] - 2*g[0][i]); // Vs gradient
            g[0][i] = 2*sqrt(pm[2][i]*(pm[0][i]+2*pm[1][i]))*g[0][i]; // Vp gradient
        }
    }
    else
    {
        #pragma omp parallel for
        for(int i=0; i<nxyz; i++)
        {
            g[1][i] = pm[1][i]*pm[0][i]*g[0][i] - pm[1][i]*pm[1][i]*g[1][i]; // rho gradient
            g[0][i] = 2*sqrt(pm[0][i]/pm[1][i])*g[0][i]; // Vp gradient
        }
    }

    if(add)
    {
        #pragma omp parallel for
        for(int i=0; i<nm*nxyz; i++) pmod[i] += temp[i];

        delete [] temp;
    }

    delete resamp;
}




void nl_we_op_vti::convert_model(data_t * m, int n, bool forward) const
{
    data_t (* pm) [n] = (data_t (*) [n]) m;
    data_t val=0;

    if (forward)
    {
        for (int i=0; i<n; i++){
            pm[0][i] = pm[2][i]*(pm[0][i]*pm[0][i]-2*pm[1][i]*pm[1][i]); // lambda
            pm[1][i] = pm[2][i]*pm[1][i]*pm[1][i]; // mu
            pm[3][i] = sqrt(2*(pm[0][i]+2*pm[1][i])*(pm[0][i]+pm[1][i])*pm[3][i] + (pm[0][i]+pm[1][i])*(pm[0][i]+pm[1][i])) - pm[1][i]; // c13
            pm[5][i] = (1+2*pm[5][i])*pm[1][i]; // c66
        }
    }
    else
    {
        for (int i=0; i<n; i++){
            pm[3][i] = ((pm[3][i]+pm[1][i])*(pm[3][i]+pm[1][i]) - (pm[0][i]+pm[1][i])*(pm[0][i]+pm[1][i])) / (2*(pm[0][i]+pm[1][i])*(pm[0][i]+2*pm[1][i])); // delta
            pm[5][i] = 0.5*(pm[5][i]/pm[1][i] - 1); // gamma
            pm[0][i] = sqrt((pm[0][i]+2*pm[1][i])/pm[2][i]); // vp
            pm[1][i] = sqrt(pm[1][i]/pm[2][i]); // vs
        }
    }
}

void nl_we_op_vti::compute_gradients(const data_t * model, const data_t * u_full, const data_t * curr, const data_t * u_x,  const data_t * u_y, const data_t * u_z, data_t * tmp, data_t * grad, const param &par, int nx, int ny, int nz, int it, data_t dx, data_t dy, data_t dz, data_t dt) const
{
    // see notes for gradients expression
    // The H quadrature will be applied elsewhere to the final gradients (all shots included)

    int nxyz = nx*ny*nz;
    int itlocal = par.nt/par.sub-it;
    data_t (* __restrict pm) [nxyz] = (data_t (*)[nxyz]) model;
    data_t (* __restrict pfor) [3][nxyz] = (data_t (*) [3][nxyz]) u_full;
    data_t (* __restrict padj) [nxyz] = (data_t (*)[nxyz]) curr;
    data_t (* __restrict padj_x) [nxyz] = (data_t (*)[nxyz]) u_x;
    data_t (* __restrict padj_y) [nxyz] = (data_t (*)[nxyz]) u_y;
    data_t (* __restrict padj_z) [nxyz] = (data_t (*)[nxyz]) u_z;
    data_t (* __restrict temp) [nxyz] = (data_t (*)[nxyz]) tmp;
    data_t (* __restrict g) [nxyz] = (data_t (*)[nxyz]) grad;

    Dx(false, pfor[itlocal][0], temp[0], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz); // forwardx_x
    Dy(false, pfor[itlocal][1], temp[1], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz); // forwardy_y
    Dz(false, pfor[itlocal][2], temp[2], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz); // forwardz_z
    Dx(false, pfor[itlocal][2], temp[3], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz); // forwardz_x
    Dz(false, pfor[itlocal][0], temp[4], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz); // forwardx_z
    Dx(false, pfor[itlocal][1], temp[5], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz); // forwardy_x
    Dy(false, pfor[itlocal][0], temp[6], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz); // forwardx_y
    Dy(false, pfor[itlocal][2], temp[7], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz); // forwardz_y
    Dz(false, pfor[itlocal][1], temp[8], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz); // forwardy_z

    data_t val1=0, val2=0, val3=0, val4=0, del=0, gamma=0;

    // different from time zero
    if (itlocal>0){
        #pragma omp parallel for private(val1,val2,val3,val4,del,gamma)
        for (int i=0; i<nxyz; i++) 
        {
            del = ((pm[3][i]+pm[1][i])*(pm[3][i]+pm[1][i]) - (pm[0][i]+pm[1][i])*(pm[0][i]+pm[1][i])) / (2*(pm[0][i]+pm[1][i])*(pm[0][i]+2*pm[1][i]));((pm[3][i]+pm[1][i])*(pm[3][i]+pm[1][i]) - (pm[0][i]+pm[1][i])*(pm[0][i]+pm[1][i])) / (2*(pm[0][i]+pm[1][i])*(pm[0][i]+2*pm[1][i]));
            gamma = (pm[5][i] - pm[1][i]) / (2*pm[1][i]);
            val1 = sqrt(2*(pm[0][i]+2*pm[1][i])*(pm[0][i]+pm[1][i])*del + (pm[0][i]+pm[1][i])*(pm[0][i]+pm[1][i]));
            val2 = ((1+2*del)*pm[0][i] + (1+3*del)*pm[1][i])/val1; // d(C13)/d(lambda)
            val3 = ((1+3*del)*pm[0][i] + (1+4*del)*pm[1][i])/val1 - 1; // d(C13)/d(mu)
            val4 = (padj_x[0][i] + padj_y[1][i])*temp[2][i] + padj_z[2][i]*(temp[0][i] + temp[1][i]);

            g[0][i] += dt*( (1+2*pm[4][i])*(padj_x[0][i] + padj_y[1][i])*(temp[0][i] + temp[1][i])
                            + padj_z[2][i]*temp[2][i] + val2*val4 ); // lambda gradient
            g[1][i] += dt*( 2*(1+2*pm[4][i])*(padj_x[0][i]*temp[0][i] + padj_y[1][i]*temp[1][i])
                            + 2*padj_z[2][i]*temp[2][i]
                            + 4*(pm[4][i] - gamma)*(padj_x[0][i]*temp[1][i] + padj_x[0][i]*temp[0][i])
                            + (padj_y[2][i] + padj_z[1][i])*(temp[7][i] + temp[8][i])
                            + (padj_x[2][i] + padj_z[0][i])*(temp[3][i] + temp[4][i])
                            + (1+2*gamma)*(padj_x[1][i] + padj_y[0][i])*(temp[5][i] + temp[6][i])
                            + val3*val4 ); // mu gradient
            g[2][i] += 1.0/dt*(padj[0][i]*(pfor[itlocal+1][0][i]-2*pfor[itlocal][0][i]+pfor[itlocal-1][0][i]) 
                            + padj[1][i]*(pfor[itlocal+1][1][i]-2*pfor[itlocal][1][i]+pfor[itlocal-1][1][i])
                            + padj[2][i]*(pfor[itlocal+1][2][i]-2*pfor[itlocal][2][i]+pfor[itlocal-1][2][i]) ); // rho gradient
        }
    }
    // time zero
    else{
        #pragma omp parallel for private(val1,val2,val3,val4,del,gamma)
        for (int i=0; i<nxyz; i++) 
        {
            del = ((pm[3][i]+pm[1][i])*(pm[3][i]+pm[1][i]) - (pm[0][i]+pm[1][i])*(pm[0][i]+pm[1][i])) / (2*(pm[0][i]+pm[1][i])*(pm[0][i]+2*pm[1][i]));((pm[3][i]+pm[1][i])*(pm[3][i]+pm[1][i]) - (pm[0][i]+pm[1][i])*(pm[0][i]+pm[1][i])) / (2*(pm[0][i]+pm[1][i])*(pm[0][i]+2*pm[1][i]));
            gamma = (pm[5][i] - pm[1][i]) / (2*pm[1][i]);
            val1 = sqrt(2*(pm[0][i]+2*pm[1][i])*(pm[0][i]+pm[1][i])*del + (pm[0][i]+pm[1][i])*(pm[0][i]+pm[1][i]));
            val2 = ((1+2*del)*pm[0][i] + (1+3*del)*pm[1][i])/val1; // d(C13)/d(lambda)
            val3 = ((1+3*del)*pm[0][i] + (1+4*del)*pm[1][i])/val1 - 1; // d(C13)/d(mu)
            val4 = (padj_x[0][i] + padj_y[1][i])*temp[2][i] + padj_z[2][i]*(temp[0][i] + temp[1][i]);

            g[0][i] += 0.5*dt*( (1+2*pm[4][i])*(padj_x[0][i] + padj_y[1][i])*(temp[0][i] + temp[1][i])
                            + padj_z[2][i]*temp[2][i] + val2*val4 ); // lambda gradient
            g[1][i] += 0.5*dt*( 2*(1+2*pm[4][i])*(padj_x[0][i]*temp[0][i] + padj_y[1][i]*temp[1][i])
                            + 2*padj_z[2][i]*temp[2][i]
                            + 4*(pm[4][i] - gamma)*(padj_x[0][i]*temp[1][i] + padj_x[0][i]*temp[0][i])
                            + (padj_y[2][i] + padj_z[1][i])*(temp[7][i] + temp[8][i])
                            + (padj_x[2][i] + padj_z[0][i])*(temp[3][i] + temp[4][i])
                            + (1+2*gamma)*(padj_x[1][i] + padj_y[0][i])*(temp[5][i] + temp[6][i])
                            + val3*val4 ); // mu gradient
            g[2][i] += 1.0/dt*( padj[0][i]*(pfor[itlocal+1][0][i]-pfor[itlocal][0][i]) 
                                + padj[1][i]*(pfor[itlocal+1][1][i]-pfor[itlocal][1][i])
                                + padj[2][i]*(pfor[itlocal+1][2][i]-pfor[itlocal][2][i]) ); // rho gradient
        }
    }
}

void nl_we_op_vti::propagate(bool adj, const data_t * model, const data_t * src, data_t * rcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int ny, int nz, data_t dx, data_t dy, data_t dz, data_t ox, data_t oy, data_t oz) const
{
    // wavefields allocations and pointers
    int nxyz = nx*ny*nz;
    data_t * u  = new data_t[9*nxyz];
    data_t * dux  = new data_t[3*nxyz];
    data_t * duy  = new data_t[3*nxyz];
    data_t * duz  = new data_t[3*nxyz];
    memset(u, 0, 9*nxyz*sizeof(data_t));
    memset(dux, 0, 3*nxyz*sizeof(data_t));
    memset(duy, 0, 3*nxyz*sizeof(data_t));
    memset(duz, 0, 3*nxyz*sizeof(data_t));
    data_t (* __restrict prev) [nxyz] = (data_t (*)[nxyz]) u;
    data_t (* __restrict curr) [nxyz] = (data_t (*)[nxyz]) (u + 3*nxyz);
    data_t (* __restrict next) [nxyz] = (data_t (*)[nxyz]) (u + 6*nxyz);
    data_t (* __restrict bucket) [nxyz];
    data_t (* __restrict u_x) [nxyz] = (data_t (*)[nxyz]) dux;
    data_t (* __restrict u_y) [nxyz] = (data_t (*)[nxyz]) duy;
    data_t (* __restrict u_z) [nxyz] = (data_t (*)[nxyz]) duz;
    data_t (* __restrict u_full) [3][nxyz];
    data_t * __restrict tmp;

    if (par.sub>0) u_full = (data_t (*) [3][nxyz]) full_wfld;
    if (grad != nullptr) 
    {
        tmp = new data_t[9*nxyz];
        memset(tmp, 0, 9*nxyz*sizeof(data_t));
    }
    else {
        tmp = new data_t[nxyz];
        memset(tmp, 0, nxyz*sizeof(data_t));
    }
	const data_t* mod[6] = {model, model+nxyz, model+2*nxyz, model+3*nxyz, model+4*nxyz, model+5*nxyz};
    
    
    // source and receiver components for injection/extraction
    int nscomp = par.nscomp;
    int nrcomp = par.nrcomp;
    if (adj) {nscomp = par.nrcomp; nrcomp=par.nscomp;}

    // prev = 1/(2 * rho) * dt2 * src
    inj->inject(false, src, prev[0], nx, ny, nz, par.nt, inj->_npts, 0, 0, 0, inj->_npts, inj->_xind.data(), inj->_yind.data(), inj->_zind.data(), inj->_xw.data(), inj->_yw.data(), inj->_zw.data());
    inj->inject(false, src, prev[1], nx, ny, nz, par.nt, inj->_npts, 1, 0, 0, inj->_npts, inj->_xind.data(), inj->_yind.data(), inj->_zind.data(), inj->_xw.data(), inj->_yw.data(), inj->_zw.data());
    inj->inject(false, src, prev[2], nx, ny, nz, par.nt, inj->_npts, 2, 0, 0, inj->_npts, inj->_xind.data(), inj->_yind.data(), inj->_zind.data(), inj->_xw.data(), inj->_yw.data(), inj->_zw.data());
    for (int i=0; i<nxyz; i++)
    {
        prev[0][i]  *= 0.5 * par.dt*par.dt/ mod[2][i];
        prev[1][i]  *= 0.5 * par.dt*par.dt/ mod[2][i];
        prev[2][i]  *= 0.5 * par.dt*par.dt/ mod[2][i];
    }    

    int pct10 = round(par.nt/10);

    for (int it=0; it<par.nt-1; it++)
    {
        // copy the current wfld to the full wfld vector
        if ((par.sub>0) && (it%par.sub==0) && (grad==nullptr)) memcpy(u_full[it/par.sub], curr, 3*nxyz*sizeof(data_t));

        // extract receivers
        if (grad == nullptr)
        {
            ext->extract(true, curr[0], rcv, nx, ny, nz, par.nt, ext->_npts, 0, it, 0, ext->_npts, ext->_xind.data(), ext->_yind.data(), ext->_zind.data(), ext->_xw.data(), ext->_yw.data(), ext->_zw.data());
            ext->extract(true, curr[1], rcv, nx, ny, nz, par.nt, ext->_npts, 1, it, 0, ext->_npts, ext->_xind.data(), ext->_yind.data(), ext->_zind.data(), ext->_xw.data(), ext->_yw.data(), ext->_zw.data());
            ext->extract(true, curr[2], rcv, nx, ny, nz, par.nt, ext->_npts, 2, it, 0, ext->_npts, ext->_xind.data(), ext->_yind.data(), ext->_zind.data(), ext->_xw.data(), ext->_yw.data(), ext->_zw.data());
        }

        // compute FWI gradients except for first and last time samples
        if ((grad != nullptr) && (par.sub>0) && (it%par.sub==0) && it!=0) compute_gradients(model, full_wfld, curr[0], u_x[0], u_y[0], u_z[0], tmp, grad, par, nx, ny, nz, it/par.sub, dx, dy, dz, par.sub*par.dt);

        // apply spatial SBP operators
        Dxx_var<c11b>(false, curr[0], next[0], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz, mod, 1.0);
        Dyy_var<c66>(true,  curr[0], next[0], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz, mod, 1.0);
        Dzz_var<mu>(true,  curr[0], next[0], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz, mod, 1.0);

        Dxx_var<c66>(false, curr[1], next[1], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz, mod, 1.0);
        Dyy_var<c11b>(true,  curr[1], next[1], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz, mod, 1.0);
        Dzz_var<mu>(true,  curr[1], next[1], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz, mod, 1.0);
        
        Dxx_var<mu>(false, curr[2], next[2], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz, mod, 1.0);
        Dyy_var<mu>(true,  curr[2], next[2], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz, mod, 1.0);
        Dzz_var<mu>(true,  curr[2], next[2], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz, mod, 2.0);

        
       
        Dx(false, curr[0], u_x[0], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz);
        Dy(false, curr[1], u_y[1], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz);
        Dz(false, curr[2], u_z[2], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz);
        #pragma omp parallel for
        for (int i=0; i<nxyz; i++)
        {
            tmp[i] = (1+2*mod[4][i])*mod[0][i]*u_x[0][i] + ((1+2*mod[4][i])*(mod[0][i]+2*mod[1][i])-2*mod[5][i])*u_y[1][i] + mod[3][i]*u_z[2][i];
        }
        Dx(true, tmp, next[0], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz);
        #pragma omp parallel for
        for (int i=0; i<nxyz; i++)
        {
            tmp[i] = (1+2*mod[4][i])*mod[0][i]*u_y[1][i] + ((1+2*mod[4][i])*(mod[0][i]+2*mod[1][i])-2*mod[5][i])*u_x[0][i] + mod[3][i]*u_z[2][i];
        }
        Dy(true, tmp, next[1], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz);
        #pragma omp parallel for
        for (int i=0; i<nxyz; i++)
        {
            tmp[i] = mod[3][i]*(u_x[0][i]+u_y[1][i]) + mod[0][i]*u_z[2][i];
        }
        Dz(true, tmp, next[2], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz);
    


        Dx(false, curr[1], u_x[1], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz);
        Dx(false, curr[2], u_x[2], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz);     
        mult_Dy(true, u_x[1], next[0], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz, mod[5]);
        mult_Dz(true, u_x[2], next[0], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz, mod[1]);

        Dy(false, curr[0], u_y[0], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz);
        Dy(false, curr[2], u_y[2], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz);     
        mult_Dx(true, u_y[0], next[1], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz, mod[5]);
        mult_Dz(true, u_y[2], next[1], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz, mod[1]);

        Dz(false, curr[0], u_z[0], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz);
        Dz(false, curr[1], u_z[1], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz);     
        mult_Dx(true, u_z[0], next[2], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz, mod[1]);
        mult_Dy(true, u_z[1], next[2], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz, mod[1]);

        // inject sources
        inj->inject(true, src, next[0], nx, ny, nz, par.nt, inj->_npts, 0, it, 0, inj->_npts, inj->_xind.data(), inj->_yind.data(), inj->_zind.data(), inj->_xw.data(), inj->_yw.data(), inj->_zw.data());
        inj->inject(true, src, next[1], nx, ny, nz, par.nt, inj->_npts, 1, it, 0, inj->_npts, inj->_xind.data(), inj->_yind.data(), inj->_zind.data(), inj->_xw.data(), inj->_yw.data(), inj->_zw.data());
        inj->inject(true, src, next[2], nx, ny, nz, par.nt, inj->_npts, 2, it, 0, inj->_npts, inj->_xind.data(), inj->_yind.data(), inj->_zind.data(), inj->_xw.data(), inj->_yw.data(), inj->_zw.data());


        // apply boundary conditions
        if (par.bc_top==1)
        {
            const data_t * in[4] = {curr[2],curr[1],curr[0],prev[0]};
            esat_neumann_top<mu,zero,zero,mu,zero>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
            in[2] = curr[1]; in[1] = curr[2];
            esat_neumann_top<zero,mu,zero,mu,zero>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
            in[0] = curr[0]; in[1] = curr[1]; in[2] = curr[2];
            esat_neumann_top<c13,c13,lam,mu2,zero>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
        }
        else if (par.bc_top==2)
        {
            const data_t * in[4] = {curr[2],curr[1],curr[0],prev[0]};
            esat_neumann_top<mu,zero,zero,mu,simp>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
            in[2] = curr[1]; in[1] = curr[2]; in[3]=prev[1];
            esat_neumann_top<zero,mu,zero,mu,simp>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
            in[0] = curr[0]; in[1] = curr[1]; in[2] = curr[2]; in[3]=prev[2];
            esat_neumann_top<c13,c13,lam,mu2,pimp>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
        }
        if (par.bc_bottom==1)
        {
            const data_t * in[4] = {curr[2],curr[1],curr[0],prev[0]};
            esat_neumann_bottom<mu,zero,zero,mu,zero>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
            in[2] = curr[1]; in[1] = curr[2];
            esat_neumann_bottom<zero,mu,zero,mu,zero>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
            in[0] = curr[0]; in[1] = curr[1]; in[2] = curr[2];
            esat_neumann_bottom<c13,c13,lam,mu2,zero>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
        }
        else if (par.bc_bottom==2)
        {
            const data_t * in[4] = {curr[2],curr[1],curr[0],prev[0]};
            esat_neumann_bottom<mu,zero,zero,mu,simp>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
            in[2] = curr[1]; in[1] = curr[2]; in[3]=prev[1];
            esat_neumann_bottom<zero,mu,zero,mu,simp>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
            in[0] = curr[0]; in[1] = curr[1]; in[2] = curr[2]; in[3]=prev[2];
            esat_neumann_bottom<c13,c13,lam,mu2,pimp>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
        }
        if (par.bc_left==1)
        {
            const data_t * in[4] = {curr[0],curr[2],curr[1],prev[0]};
            esat_neumann_left<c66,zero,zero,c66,zero>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
            in[2] = curr[2]; in[1] = curr[0];
            esat_neumann_left<zero,mu,zero,mu,zero>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
            in[0] = curr[1]; in[1] = curr[2]; in[2] = curr[0];
            esat_neumann_left<c12,c13,c11c,c11b,zero>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
        }
        else if (par.bc_left==2)
        {
            const data_t * in[4] = {curr[0],curr[2],curr[1],prev[1]};
            esat_neumann_left<c66,zero,zero,c66,simp>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
            in[2] = curr[2]; in[1] = curr[0]; in[3]=prev[2];
            esat_neumann_left<zero,mu,zero,mu,simp>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
            in[0] = curr[1]; in[1] = curr[2]; in[2] = curr[0]; in[3]=prev[0];
            esat_neumann_left<c12,c13,c11c,c11b,pimp2>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
        }
        if (par.bc_right==1)
        {
            const data_t * in[4] = {curr[0],curr[2],curr[1],prev[0]};
            esat_neumann_right<c66,zero,zero,c66,zero>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
            in[2] = curr[2]; in[1] = curr[0];
            esat_neumann_right<zero,mu,zero,mu,zero>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
            in[0] = curr[1]; in[1] = curr[2]; in[2] = curr[0];
            esat_neumann_right<c12,c13,c11c,c11b,zero>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
        }
        else if (par.bc_right==2)
        {
            const data_t * in[4] = {curr[0],curr[2],curr[1],prev[1]};
            esat_neumann_right<c66,zero,zero,c66,simp>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
            in[2] = curr[2]; in[1] = curr[0]; in[3]=prev[2];
            esat_neumann_right<zero,mu,zero,mu,simp>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
            in[0] = curr[1]; in[1] = curr[2]; in[2] = curr[0]; in[3]=prev[0];
            esat_neumann_right<c12,c13,c11c,c11b,pimp2>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
        }
        if (par.bc_front==1)
        {
            const data_t * in[4] = {curr[1],curr[2],curr[0],prev[0]};
            esat_neumann_front<c66,zero,zero,c66,zero>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
            in[2] = curr[2]; in[1] = curr[1];
            esat_neumann_front<zero,mu,zero,mu,zero>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
            in[0] = curr[0]; in[1] = curr[2]; in[2] = curr[1];
            esat_neumann_front<c12,c13,c11c,c11b,zero>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
        }
        else if (par.bc_front==2)
        {
            const data_t * in[4] = {curr[1],curr[2],curr[0],prev[0]};
            esat_neumann_front<c66,zero,zero,c66,simp>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
            in[2] = curr[2]; in[1] = curr[1]; in[3]=prev[2];
            esat_neumann_front<zero,mu,zero,mu,simp>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
            in[0] = curr[0]; in[1] = curr[2]; in[2] = curr[1]; in[3]=prev[1];
            esat_neumann_front<c12,c13,c11c,c11b,pimp2>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
        }
        if (par.bc_back==1)
        {
            const data_t * in[4] = {curr[1],curr[2],curr[0],prev[0]};
            esat_neumann_back<c66,zero,zero,c66,zero>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
            in[2] = curr[2]; in[1] = curr[1];
            esat_neumann_back<zero,mu,zero,mu,zero>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
            in[0] = curr[0]; in[1] = curr[2]; in[2] = curr[1];
            esat_neumann_back<c12,c13,c11c,c11b,zero>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
        }
        else if (par.bc_back==2)
        {
            const data_t * in[4] = {curr[1],curr[2],curr[0],prev[0]};
            esat_neumann_back<c66,zero,zero,c66,simp>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
            in[2] = curr[2]; in[1] = curr[1]; in[3]=prev[2];
            esat_neumann_back<zero,mu,zero,mu,simp>(true, in, next[2], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
            in[0] = curr[0]; in[1] = curr[2]; in[2] = curr[1]; in[3]=prev[1];
            esat_neumann_back<c12,c13,c11c,c11b,pimp2>(true, in, next[1], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
        }

        // update wfld with the 2 steps time recursion
        #pragma omp parallel for
        for (int i=0; i<nxyz; i++)
        {
            next[0][i] = par.dt*par.dt*next[0][i]/mod[2][i] + 2*curr[0][i] - prev[0][i];
            next[1][i] = par.dt*par.dt*next[1][i]/mod[2][i] + 2*curr[1][i] - prev[1][i];
            next[2][i] = par.dt*par.dt*next[2][i]/mod[2][i] + 2*curr[2][i] - prev[2][i];
        }

        // scale boundaries when relevant (for locally absorbing BC only)
        data_t * in[3] = {next[0],next[1],next[2]};
        esat_scale_boundaries<eps2>(in, nx, ny, nz, dx, dy, dz, 0, nx, 0, ny, 0, nz, mod, par.dt, par.bc_top==2, par.bc_bottom==2, par.bc_left==2, par.bc_right==2, par.bc_front==2, par.bc_back==2);
        
        // apply taper 
        taperz(curr[0], nx, ny, nz, 3, 0, nx, 0, ny, par.taper_top, 0, 0, 3, par.taper_strength);
        taperz(next[0], nx, ny, nz, 3, 0, nx, 0, ny, par.taper_top, 0, 0, 3, par.taper_strength);
        taperz(curr[0], nx, ny, nz, 3, 0, nx, 0, ny, nz-par.taper_bottom, nz, 0, 3, par.taper_strength);
        taperz(next[0], nx, ny, nz, 3, 0, nx, 0, ny, nz-par.taper_bottom, nz, 0, 3, par.taper_strength);
        taperx(curr[0], nx, ny, nz, 3, par.taper_left, 0, 0, ny, 0, nz, 0, 3, par.taper_strength);
        taperx(next[0], nx, ny, nz, 3, par.taper_left, 0, 0, ny, 0, nz, 0, 3, par.taper_strength);
        taperx(curr[0], nx, ny, nz, 3, nx-par.taper_right, nx, 0, ny, 0, nz, 0, 3, par.taper_strength);
        taperx(next[0], nx, ny, nz, 3, nx-par.taper_right, nx, 0, ny, 0, nz, 0, 3, par.taper_strength);
        tapery(curr[0], nx, ny, nz, 3, 0, nx, par.taper_front, 0, 0, nz, 0, 3, par.taper_strength);
        tapery(next[0], nx, ny, nz, 3, 0, nx, par.taper_front, 0, 0, nz, 0, 3, par.taper_strength);
        tapery(curr[0], nx, ny, nz, 3, 0, nx, ny-par.taper_back, ny, 0, nz, 0, 3, par.taper_strength);
        tapery(next[0], nx, ny, nz, 3, 0, nx, ny-par.taper_back, ny, 0, nz, 0, 3, par.taper_strength);
       
        bucket=prev;
        prev=curr;
        curr=next;
        next=bucket;

        if ((it+1) % pct10 == 0 && par.verbose>2) fprintf(stderr,"Propagation progress = %d\%\n",10*(it+1)/pct10);
    }

    // copy the last wfld to the full wfld vector
    if ((par.sub>0) && ((par.nt-1)%par.sub==0) && (grad==nullptr)) memcpy(u_full[(par.nt-1)/par.sub], curr, 3*nxyz*sizeof(data_t));

    // extract receivers last sample
    if (grad == nullptr)
    {
        ext->extract(true, curr[0], rcv, nx, ny, nz, par.nt, ext->_npts, 0, par.nt-1, 0, ext->_npts, ext->_xind.data(), ext->_yind.data(), ext->_zind.data(), ext->_xw.data(), ext->_yw.data(), ext->_zw.data());
        ext->extract(true, curr[1], rcv, nx, ny, nz, par.nt, ext->_npts, 1, par.nt-1, 0, ext->_npts, ext->_xind.data(), ext->_yind.data(), ext->_zind.data(), ext->_xw.data(), ext->_yw.data(), ext->_zw.data());
        ext->extract(true, curr[2], rcv, nx, ny, nz, par.nt, ext->_npts, 2, par.nt-1, 0, ext->_npts, ext->_xind.data(), ext->_yind.data(), ext->_zind.data(), ext->_xw.data(), ext->_yw.data(), ext->_zw.data());
    }

    // last sample gradient
    if ((grad != nullptr) && (par.sub>0) && ((par.nt-1)%par.sub==0) ) compute_gradients(model, full_wfld, curr[0], u_x[0], u_y[0], u_z[0], tmp, grad, par, nx, ny, nz, (par.nt-1)/par.sub, dx, dy, dz, par.sub*par.dt);
    
    delete [] u;
    delete [] dux;
    delete [] duy;
    delete [] duz;
    delete [] tmp;
}



void nl_we_op_a::convert_model(data_t * m, int n, bool forward) const
{
    data_t (* pm) [n] = (data_t (*) [n]) m;

    if (forward)
    {
        for (int i=0; i<n; i++){
            pm[0][i] = pm[1][i]*pm[0][i]*pm[0][i];
            pm[1][i] = 1.0/pm[1][i];
        }
    }
    else
    {
        for (int i=0; i<n; i++){
            pm[0][i] = sqrt(pm[0][i]*pm[1][i]);
            pm[1][i] = 1.0/pm[1][i];
        }
    }
}

void nl_we_op_a::compute_gradients(const data_t * model, const data_t * u_full, const data_t * curr, data_t * tmp, data_t * grad, const param &par, int nx, int ny, int nz, int it, data_t dx, data_t dy, data_t dz, data_t dt) const
{
    // Grad_K-1 = adjoint.H.Ht.forward_tt
    // Grad_K = -(1/K^2).Grad_K-1
    // Grad_rho-1 = adjoint_x.H.Ht.forward_x + adjoint_y.H.Ht.forward_y + adjoint_z.H.Ht.forward_z
    // The H quadrature will be applied elsewhere to the final gradients (all shots included)

    int nxyz = nx*ny*nz;
    int itlocal = par.nt/par.sub+1-it-1;
    data_t (*pm) [nxyz] = (data_t (*)[nxyz]) model;
    data_t (*pfor) [nxyz] = (data_t (*) [nxyz]) u_full;
    const data_t *padj = curr;
    data_t (*temp) [nxyz] = (data_t (*)[nxyz]) tmp;
    data_t (*g) [nxyz] = (data_t (*)[nxyz]) grad;

    Dx(false, pfor[itlocal], temp[0], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz); // forward_x
    Dy(false, pfor[itlocal], temp[1], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz); // forward_y
    Dz(false, pfor[itlocal], temp[2], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz); // forward_z
    Dx(false, curr, temp[3], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz); // adjoint_x
    Dy(false, curr, temp[4], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz); // adjoint_y
    Dz(false, curr, temp[5], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz); // adjoint_z

    // different from time zero
    if (par.nt/par.sub-it>0){
        #pragma omp parallel for
        for (int i=0; i<nxyz; i++) 
        {
            g[0][i] -= 1.0/(dt*pm[0][i]*pm[0][i])*padj[i]*(pfor[itlocal+1][i]-2*pfor[itlocal][i]+pfor[itlocal-1][i]); // K gradient
            g[1][i] += dt*(temp[3][i]*temp[0][i] + temp[4][i]*temp[1][i] + temp[5][i]*temp[2][i]); // rho-1 gradient
        }
    }
    // time zero
    else{
        #pragma omp parallel for
        for (int i=0; i<nxyz; i++) 
        {
            g[0][i] -= 1.0/(dt*pm[0][i]*pm[0][i])*padj[i]*(pfor[itlocal+1][i]-2*pfor[itlocal][i]+pfor[itlocal-1][i]); // K gradient
            g[1][i] += 0.5*dt*(temp[3][i]*temp[0][i] + temp[4][i]*temp[1][i] + temp[5][i]*temp[2][i]); // rho-1 gradient
        }
    }
}

void nl_we_op_a::propagate(bool adj, const data_t * model, const data_t * src, data_t * rcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int ny, int nz, data_t dx, data_t dy, data_t dz, data_t ox, data_t oy, data_t oz) const
{
    // wavefields allocations and pointers
    int nxyz=nx*ny*nz;
    data_t * u  = new data_t[3*nxyz];
    memset(u, 0, 3*nxyz*sizeof(data_t));
    data_t (* __restrict prev) [nxyz] = (data_t (*)[nxyz]) u;
    data_t (* __restrict curr) [nxyz] = (data_t (*)[nxyz]) (u + nxyz);
    data_t (* __restrict next) [nxyz] = (data_t (*)[nxyz]) (u + 2*nxyz);
    data_t (* __restrict bucket) [nxyz];
    data_t (* __restrict u_full) [nxyz];
    data_t * __restrict tmp;

    if (par.sub>0) u_full = (data_t (*) [nxyz]) full_wfld;
    if (grad != nullptr) 
    {
        tmp = new data_t[6*nxyz];
        memset(tmp, 0, 6*nxyz*sizeof(data_t));
    }
	const data_t* mod[2] = {model, model+nxyz};

    // free surface stiffness
    data_t alpha = 0.2508560249 / par.free_surface_stiffness;

    // prev = K/2 * dt2 * src
    inj->inject(false, src, prev[0], nx, ny, nz, par.nt, inj->_npts, 0, 0, 0, inj->_npts, inj->_xind.data(), inj->_yind.data(), inj->_zind.data(), inj->_xw.data(), inj->_yw.data(), inj->_zw.data());
    for (int i=0; i<nxyz; i++) prev[0][i]  *= 0.5 * mod[0][i]*par.dt*par.dt;  

    int pct10 = round(par.nt/10);

    for (int it=0; it<par.nt-1; it++)
    {
        // copy the current wfld to the full wfld vector
        if ((par.sub>0) && (it%par.sub==0) && (grad==nullptr)) memcpy(u_full[it/par.sub], curr, nxyz*sizeof(data_t));

        // extract receivers
        if (grad == nullptr)
        {
            ext->extract(true, curr[0], rcv, nx, ny, nz, par.nt, ext->_npts, 0, it, 0, ext->_npts, ext->_xind.data(), ext->_yind.data(), ext->_zind.data(), ext->_xw.data(), ext->_yw.data(), ext->_zw.data());
        }

        // compute FWI gradients except for first and last time samples
        if ((grad != nullptr) && (par.sub>0) && (it%par.sub==0) && it!=0) compute_gradients(model, full_wfld, curr[0], tmp, grad, par, nx, ny, nz, it/par.sub, dx, dy, dz, par.sub*par.dt);

        // apply spatial SBP operators
        // mu here designates the second model parameter, i.e. reciprocal of density 1/rho 
        Dxx_var<mu>(false, curr[0], next[0], nx, ny, nz, dx, 0, nx, 0, ny, 0, nz, mod, 1.0);
        Dyy_var<mu>(true, curr[0], next[0], nx, ny, nz, dy, 0, nx, 0, ny, 0, nz, mod, 1.0);
        Dzz_var<mu>(true, curr[0], next[0], nx, ny, nz, dz, 0, nx, 0, ny, 0, nz, mod, 1.0);

        // inject sources
        inj->inject(true, src, next[0], nx, ny, nz, par.nt, inj->_npts, 0, it, 0, inj->_npts, inj->_xind.data(), inj->_yind.data(), inj->_zind.data(), inj->_xw.data(), inj->_yw.data(), inj->_zw.data());

        // apply boundary conditions
        if (par.bc_top==1)
        {
            const data_t * in[1] = {curr[0]};
            asat_dirichlet_top(true, in, next[0], nx, ny, nz, dx, dy, dz, 0, nx, 0, ny, mod, alpha);
        }
        else if (par.bc_top==2)
        {
            const data_t * in[2] = {curr[0],prev[0]};
            asat_neumann_top<mu,aimp>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
        }
        else if (par.bc_top==3)
        {
            const data_t * in[2] = {curr[0],prev[0]};
            asat_neumann_top<mu,zero>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
        }
        //else if (par.bc_top==4)
        //{
        //    const data_t * in[2] = {curr[0],prev[0]};
        //    asat_neumann_absorbing_top(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, _transmission->getCVals());
        //}
        if (par.bc_bottom==1)
        {
            const data_t * in[1] = {curr[0]};
            asat_dirichlet_bottom(true, in, next[0], nx, ny, nz, dx, dy, dz, 0, nx, 0, ny, mod, alpha);
        }
        else if (par.bc_bottom==2)
        {
            const data_t * in[2] = {curr[0],prev[0]};
            asat_neumann_bottom<mu,aimp>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
        }
        else if (par.bc_bottom==3)
        {
            const data_t * in[2] = {curr[0],prev[0]};
            asat_neumann_bottom<mu,zero>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, ny, mod);
        }
        if (par.bc_left==1)
        {
            const data_t * in[1] = {curr[0]};
            asat_dirichlet_left(true, in, next[0], nx, ny, nz, dx, dy, dz, 0, ny, 0, nz, mod, alpha);
        }
        else if (par.bc_left==2)
        {
            const data_t * in[2] = {curr[0],prev[0]};
            asat_neumann_left<mu,aimp>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
        }
        else if (par.bc_left==3)
        {
            const data_t * in[2] = {curr[0],prev[0]};
            asat_neumann_left<mu,zero>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
        }
        if (par.bc_right==1)
        {
            const data_t * in[1] = {curr[0]};
            asat_dirichlet_right(true, in, next[0], nx, ny, nz, dx, dy, dz, 0, ny, 0, nz, mod, alpha);
        }
        else if (par.bc_right==2)
        {
            const data_t * in[2] = {curr[0],prev[0]};
            asat_neumann_right<mu,aimp>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
        }
        else if (par.bc_right==3)
        {
            const data_t * in[2] = {curr[0],prev[0]};
            asat_neumann_right<mu,zero>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, ny, 0, nz, mod);
        }
        if (par.bc_front==1)
        {
            const data_t * in[1] = {curr[0]};
            asat_dirichlet_front(true, in, next[0], nx, ny, nz, dx, dy, dz, 0, nx, 0, nz, mod, alpha);
        }
        else if (par.bc_front==2)
        {
            const data_t * in[2] = {curr[0],prev[0]};
            asat_neumann_front<mu,aimp>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
        }
        else if (par.bc_front==3)
        {
            const data_t * in[2] = {curr[0],prev[0]};
            asat_neumann_front<mu,zero>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
        }
        if (par.bc_back==1)
        {
            const data_t * in[1] = {curr[0]};
            asat_dirichlet_back(true, in, next[0], nx, ny, nz, dx, dy, dz, 0, nx, 0, nz, mod, alpha);
        }
        else if (par.bc_back==2)
        {
            const data_t * in[2] = {curr[0],prev[0]};
            asat_neumann_back<mu,aimp>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
        }
        else if (par.bc_back==3)
        {
            const data_t * in[2] = {curr[0],prev[0]};
            asat_neumann_back<mu,zero>(true, in, next[0], nx, ny, nz, dx, dy, dz, par.dt, 0, nx, 0, nz, mod);
        }

        // update wfld with the 2 steps time recursion
        #pragma omp parallel for
        for (int i=0; i<nxyz; i++)
        {
            next[0][i] = par.dt*par.dt*next[0][i]*mod[0][i] + 2*curr[0][i] - prev[0][i];
        }

        // scale boundaries when relevant (for locally absorbing BC only)
        data_t * in[1] = {next[0]};
        if (par.bc_top!=4) asat_scale_boundaries(in, nx, ny, nz, dx, dy, dz, 0, nx, 0, ny, 0, nz, mod, par.dt, par.bc_top==2, par.bc_bottom==2, par.bc_left==2, par.bc_right==2, par.bc_front==2, par.bc_back==2);
        // else asat_scale_boundaries_bis(in, nx, ny, nz, dx, dy, dz, 0, nx, 0, ny, 0, nz, mod, _transmission->getCVals(), par.dt, par.bc_top==4, par.bc_bottom==2, par.bc_left==2, par.bc_right==2, par.bc_front==2, par.bc_back==2);
        
        taperz(curr[0], nx, ny, nz, 1, 0, nx, 0, ny, par.taper_top, 0, 0, 1, par.taper_strength);
        taperz(next[0], nx, ny, nz, 1, 0, nx, 0, ny, par.taper_top, 0, 0, 1, par.taper_strength);
        taperz(curr[0], nx, ny, nz, 1, 0, nx, 0, ny, nz-par.taper_bottom, nz, 0, 1, par.taper_strength);
        taperz(next[0], nx, ny, nz, 1, 0, nx, 0, ny, nz-par.taper_bottom, nz, 0, 1, par.taper_strength);
        taperx(curr[0], nx, ny, nz, 1, par.taper_left, 0, 0, ny, 0, nz, 0, 1, par.taper_strength);
        taperx(next[0], nx, ny, nz, 1, par.taper_left, 0, 0, ny, 0, nz, 0, 1, par.taper_strength);
        taperx(curr[0], nx, ny, nz, 1, nx-par.taper_right, nx, 0, ny, 0, nz, 0, 1, par.taper_strength);
        taperx(next[0], nx, ny, nz, 1, nx-par.taper_right, nx, 0, ny, 0, nz, 0, 1, par.taper_strength);
        tapery(curr[0], nx, ny, nz, 1, 0, nx, par.taper_front, 0, 0, nz, 0, 1, par.taper_strength);
        tapery(next[0], nx, ny, nz, 1, 0, nx, par.taper_front, 0, 0, nz, 0, 1, par.taper_strength);
        tapery(curr[0], nx, ny, nz, 1, 0, nx, ny-par.taper_back, ny, 0, nz, 0, 1, par.taper_strength);
        tapery(next[0], nx, ny, nz, 1, 0, nx, ny-par.taper_back, ny, 0, nz, 0, 1, par.taper_strength);

        bucket=prev;
        prev=curr;
        curr=next;
        next=bucket;

        if ((it+1) % pct10 == 0 && par.verbose>2) fprintf(stderr,"Propagation progress = %d\%\n",10*(it+1)/pct10);
    }
}