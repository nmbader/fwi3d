#include "bsplines.hpp"

#define ZERO 1e-16

// B-splines functions of order 0, 1, 2, 3
// N0
data_t N0(int i, data_t u, const std::vector<data_t> &uk){
    if ((u<uk[i]) || (u>=uk[i+1])) return 0;
    else return 1;
}
// N1
data_t N1(int i, data_t u, const std::vector<data_t> &uk){
    if (std::abs(uk[i+1]-uk[i])<ZERO && std::abs(uk[i+2]-uk[i+1])<ZERO) return 0;
    else if  (std::abs(uk[i+1]-uk[i])<ZERO && std::abs(uk[i+2]-uk[i+1])>=ZERO) return (uk[i+2]-u)/(uk[i+2]-uk[i+1])*N0(i+1,u,uk);
    else if (std::abs(uk[i+1]-uk[i])>=ZERO && std::abs(uk[i+2]-uk[i+1])<ZERO) return (u-uk[i])/(uk[i+1]-uk[i])*N0(i,u,uk);
    else return (u-uk[i])/(uk[i+1]-uk[i])*N0(i,u,uk) + (uk[i+2]-u)/(uk[i+2]-uk[i+1])*N0(i+1,u,uk);
}

// N2
data_t N2(int i, data_t u, const std::vector<data_t> &uk){
    if (std::abs(uk[i+2]-uk[i])<ZERO && std::abs(uk[i+3]-uk[i+1])<ZERO) return 0;
    else if  (std::abs(uk[i+2]-uk[i])<ZERO && std::abs(uk[i+3]-uk[i+1])>=ZERO) return (uk[i+3]-u)/(uk[i+3]-uk[i+1])*N1(i+1,u,uk);
    else if (std::abs(uk[i+2]-uk[i])>=ZERO && std::abs(uk[i+3]-uk[i+1])<ZERO) return (u-uk[i])/(uk[i+2]-uk[i])*N1(i,u,uk);
    else return (u-uk[i])/(uk[i+2]-uk[i])*N1(i,u,uk) + (uk[i+3]-u)/(uk[i+3]-uk[i+1])*N1(i+1,u,uk);
}

// N3
data_t N3(int i, data_t u, const std::vector<data_t> &uk){
    if (std::abs(uk[i+3]-uk[i])<ZERO && std::abs(uk[i+4]-uk[i+1])<ZERO) return 0;
    else if  (std::abs(uk[i+3]-uk[i])<ZERO && std::abs(uk[i+4]-uk[i+1])>=ZERO) return (uk[i+4]-u)/(uk[i+4]-uk[i+1])*N2(i+1,u,uk);
    else if (std::abs(uk[i+3]-uk[i])>=ZERO && std::abs(uk[i+4]-uk[i+1])<ZERO) return (u-uk[i])/(uk[i+3]-uk[i])*N2(i,u,uk);
    else return (u-uk[i])/(uk[i+3]-uk[i])*N2(i,u,uk) + (uk[i+4]-u)/(uk[i+4]-uk[i+1])*N2(i+1,u,uk);
}


void setKnot(std::vector<data_t> &u, std::vector<data_t> &c, std::vector<int> &m){
    successCheck(c.size() == m.size(),"Control and multiplicity vectors must have the same size\n");
    u.clear();
    u.push_back(c[0]);
    u.push_back(c[0]);
    for (int i=0; i<m.size(); i++){
        for (int j=0; j<m[i]; j++) u.push_back(c[i]);
    }
    u.push_back(c[c.size()-1]);
    u.push_back(c[c.size()-1]);
}

void fillin3d(std::shared_ptr<vecReg<data_t> > c, const std::shared_ptr<vecReg<data_t> > v, const std::vector<data_t> &cx, const std::vector<data_t> &cy, const std::vector<data_t> &cz){
    axis<data_t> Z = c->getHyper()->getAxis(1);
    axis<data_t> Z2 = v->getHyper()->getAxis(1);
    axis<data_t> X = c->getHyper()->getAxis(2);
    axis<data_t> X2 = v->getHyper()->getAxis(2);
    axis<data_t> Y = c->getHyper()->getAxis(3);
    axis<data_t> Y2 = v->getHyper()->getAxis(3);
    successCheck(Z.n == cz.size(),"The coarse spline vector must have the same z size as the control vector\n");
    successCheck(X.n == cx.size(),"The coarse spline vector must have the same x size as the control vector\n");
    successCheck(Y.n == cy.size(),"The coarse spline vector must have the same y size as the control vector\n");
    int nc = c->getN123()/(Y.n*X.n*Z.n);
    data_t (* __restrict pc) [Y.n][X.n][Z.n] = (data_t (*) [Y.n][X.n][Z.n]) c->getVals();
    const data_t (* pv) [Y2.n][X2.n][Z2.n] = (const data_t (*) [Y2.n][X2.n][Z2.n]) v->getCVals();
    int ix2, iy2, iz2, ix3, iy3, iz3;
    data_t wx, wy, wz;
    for (int ic=0; ic<nc; ic++){
        for (int iy=0; iy<Y.n; iy++){
            iy2 = floor((cy[iy]-Y2.o)/Y2.d);
            iy3 = std::min(Y2.n-1,iy2+1);
            wy = (cy[iy] - Y2.o - iy2*Y2.d)/Y2.d;
            for (int ix=0; ix<X.n; ix++){
                ix2 = floor((cx[ix]-X2.o)/X2.d);
                ix3 = std::min(X2.n-1,ix2+1);
                wx = (cx[ix] - X2.o - ix2*X2.d)/X2.d;
                for (int iz=0; iz<Z.n; iz++){
                    iz2 = floor((cz[iz]-Z2.o)/Z2.d);
                    iz3 = std::min(Z2.n-1,iz2+1);
                    wz = (cz[iz] - Z2.o - iz2*Z2.d)/Z2.d;
                    pc[ic][iy][ix][iz] = (1-wx)*(1-wy)*(1-wz)*pv[ic][iy2][ix2][iz2]  
                    + (wx)*(1-wy)*(1-wz)*pv[ic][iy2][ix3][iz2] 
                    + (1-wx)*(wy)*(1-wz)*pv[ic][iy3][ix2][iz2]
                    + (1-wx)*(1-wy)*(wz)*pv[ic][iy2][ix2][iz3]
                    + (wx)*(wy)*(1-wz)*pv[ic][iy3][ix3][iz2]
                    + (wx)*(1-wy)*(wz)*pv[ic][iy2][ix3][iz3]
                    + (1-wx)*(wy)*(wz)*pv[ic][iy3][ix2][iz3]
                    + (wx)*(wy)*(wz)*pv[ic][iy3][ix3][iz3];
                }
            }
        }
    }
}

void fillin2d(std::shared_ptr<vecReg<data_t> > c, const std::shared_ptr<vecReg<data_t> > v, const std::vector<data_t> &cx, const std::vector<data_t> &cz){
    axis<data_t> Z = c->getHyper()->getAxis(1);
    axis<data_t> Z2 = v->getHyper()->getAxis(1);
    axis<data_t> X = c->getHyper()->getAxis(2);
    axis<data_t> X2 = v->getHyper()->getAxis(2);
    successCheck(Z.n == cz.size(),"The coarse spline vector must have the same z size as the control vector\n");
    successCheck(X.n == cx.size(),"The coarse spline vector must have the same x size as the control vector\n");
    int ny = c->getN123()/(X.n*Z.n);
    data_t (*pc) [X.n][Z.n] = (data_t (*) [X.n][Z.n]) c->getVals();
    const data_t (*pv) [X2.n][Z2.n] = (const data_t (*) [X2.n][Z2.n]) v->getCVals();
    int ix2, iz2, ix3, iz3;
    data_t wx, wz;
    for (int iy=0; iy<ny; iy++){
        for (int ix=0; ix<X.n; ix++){
            ix2 = floor((cx[ix]-X2.o)/X2.d);
            ix3 = std::min(X2.n-1,ix2+1);
            wx = (cx[ix] - X2.o - ix2*X2.d)/X2.d;
            for (int iz=0; iz<Z.n; iz++){
                iz2 = floor((cz[iz]-Z2.o)/Z2.d);
                iz3 = std::min(Z2.n-1,iz2+1);
                wz = (cz[iz] - Z2.o - iz2*Z2.d)/Z2.d;
                pc[iy][ix][iz] = (1-wx)*(1-wz)*pv[iy][ix2][iz2]  + wx*(1-wz)*pv[iy][ix3][iz2] + (1-wx)*wz*pv[iy][ix2][iz3] + wx*wz*pv[iy][ix3][iz3];
            }
        }
    }
}

void fillin1d(std::shared_ptr<vecReg<data_t> > c, const std::shared_ptr<vecReg<data_t> > v, const std::vector<data_t> &cz){
    axis<data_t> Z = c->getHyper()->getAxis(1);
    axis<data_t> Z2 = v->getHyper()->getAxis(1);
    successCheck(Z.n == cz.size(),"The coarse spline vector must have the same z size as the control vector\n");
    int ny = c->getN123()/(Z.n);
    data_t * pc = c->getVals();
    const data_t * pv = v->getCVals();
    int iz2, iz3;
    data_t wz;
    for (int iy=0; iy<ny; iy++){
        for (int iz=0; iz<Z.n; iz++){
            iz2 = floor((cz[iz]-Z2.o)/Z2.d);
            iz3 = std::min(Z2.n-1,iz2+1);
            wz = (cz[iz] - Z2.o - iz2*Z2.d)/Z2.d;
            pc[iy*Z.n+iz] = (1-wz)*pv[iy*Z2.n+iz2] + wz*pv[iy*Z2.n+iz3];
        }
    }
}

void duplicate3d::apply_forward(bool add, const data_t * pmod, data_t * pdat) {
    
    int nz = _domain.getAxis(1).n;
    int nz2 = _range.getAxis(1).n;
    int nx = _domain.getAxis(2).n;
    int nx2 = _range.getAxis(2).n;
    int ny = _domain.getAxis(3).n;
    int ny2 = _range.getAxis(3).n;
    int nc = _domain.getN123()/(nx*ny*nz);

    const data_t (* pm) [ny][nx][nz] = (const data_t (*) [ny][nx][nz]) pmod; 
    data_t (* __restrict pd) [ny2][nx2][nz2] = (data_t (*) [ny2][nx2][nz2]) pdat;

    #pragma omp parallel for
    for (int ic=0; ic<nc; ic++){
        int county=0;
        for (int iy=0; iy<ny; county+=_my[iy], iy++){
            for (int jy=county; jy<county+_my[iy]; jy++){
                int countx=0;
                for (int ix=0; ix<nx; countx+=_mx[ix], ix++){
                    for (int jx=countx; jx<countx+_mx[ix]; jx++){
                        int countz=0;
                        for (int iz=0; iz<nz; countz+=_mz[iz], iz++){
                            for (int jz=countz; jz<countz+_mz[iz]; jz++){
                                pd[ic][jy][jx][jz] = add*pd[ic][jy][jx][jz] + pm[ic][iy][ix][iz];
                            }
                        }
                    }
                }
            }
        }
    }
}

void duplicate3d::apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

    int nz = _domain.getAxis(1).n;
    int nz2 = _range.getAxis(1).n;
    int nx = _domain.getAxis(2).n;
    int nx2 = _range.getAxis(2).n;
    int ny = _domain.getAxis(3).n;
    int ny2 = _range.getAxis(3).n;
    int nc = _domain.getN123()/(nx*ny*nz);

    if (!add) memset(pmod, 0, _domain.getN123()*sizeof(data_t));

    data_t (* __restrict pm) [ny][nx][nz] = (data_t (*) [ny][nx][nz]) pmod; 
    const data_t (* pd) [ny2][nx2][nz2] = (const data_t (*) [ny2][nx2][nz2]) pdat;

    #pragma omp parallel for
    for (int ic=0; ic<nc; ic++){
        int county=0;
        for (int iy=0; iy<ny; county+=_my[iy], iy++){
            for (int jy=county; jy<county+_my[iy]; jy++){
                int countx=0;
                for (int ix=0; ix<nx; countx+=_mx[ix], ix++){
                    for (int jx=countx; jx<countx+_mx[ix]; jx++){
                        int countz=0;
                        for (int iz=0; iz<nz; countz+=_mz[iz], iz++){
                            for (int jz=countz; jz<countz+_mz[iz]; jz++){
                                pm[ic][iy][ix][iz] += pd[ic][jy][jx][jz];
                            }
                        }
                    }
                }
            }
        }
    }
}

void duplicate2d::apply_forward(bool add, const data_t * pmod, data_t * pdat) {
    
    int nz = _domain.getAxis(1).n;
    int nz2 = _range.getAxis(1).n;
    int nx = _domain.getAxis(2).n;
    int nx2 = _range.getAxis(2).n;
    int ny = _domain.getN123()/(nx*nz);

    data_t * p = new data_t[ny*nx*nz2];

    int count;
    #pragma omp parallel for private(count)
    for (int i=0; i<ny*nx; i++){
        int count = 0;
        for (int iz=0; iz<nz; iz++){
            for (int j=count; j<count+_mz[iz]; j++) p[i*nz2+j] = pmod[i*nz+iz];
            count += _mz[iz];
        }
    }
    #pragma omp parallel for private(count)
    for (int iy=0; iy<ny; iy++){
        count = 0;
        for (int ix=0; ix<nx; ix++){
            for (int j=count; j<count+_mx[ix]; j++){
                for (int iz=0; iz<nz2; iz++) pdat[iy*nx2*nz2+j*nz2+iz] = add*pdat[iy*nx2*nz2+j*nz2+iz] + p[iy*nx*nz2+ix*nz2+iz];
            }
            count += _mx[ix];
        }
    }
    delete [] p;
}
void duplicate2d::apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

    int nz = _domain.getAxis(1).n;
    int nz2 = _range.getAxis(1).n;
    int nx = _domain.getAxis(2).n;
    int nx2 = _range.getAxis(2).n;
    int ny = _domain.getN123()/(nx*nz);

    if (add == false) memset(pmod, 0, _domain.getN123()*sizeof(data_t));

    data_t * p = new data_t[ny*nx2*nz];
    int count;
    #pragma omp parallel for private(count)
    for (int i=0; i<ny*nx2; i++){
        count = 0;
        for (int iz=0; iz<nz; iz++){
            p[i*nz+iz] = 0;
            for (int j=count; j<count+_mz[iz]; j++) p[i*nz+iz] += pdat[i*nz2+j];
            count += _mz[iz];
        }
    }
    #pragma omp parallel for private(count)
    for (int iy=0; iy<ny; iy++){
        count = 0;
        for (int ix=0; ix<nx; ix++){
            for (int j=count; j<count+_mx[ix]; j++){
                for (int iz=0; iz<nz; iz++) pmod[iy*nx*nz+ix*nz+iz] += p[iy*nx2*nz+j*nz+iz];
            }
            count += _mx[ix];
        }
    }
    delete [] p;
}

void duplicate1d::apply_forward(bool add, const data_t * pmod, data_t * pdat) {
    
    int nz = _domain.getAxis(1).n;
    int nz2 = _range.getAxis(1).n;
    int ny = _domain.getN123()/(nz);

    int count;
    #pragma omp parallel for private(count)
    for (int i=0; i<ny; i++){
        int count = 0;
        for (int iz=0; iz<nz; iz++){
            for (int j=count; j<count+_mz[iz]; j++) pdat[i*nz2+j] = add*pdat[i*nz2+j] + pmod[i*nz+iz];
            count += _mz[iz];
        }
    }
}
void duplicate1d::apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

    int nz = _domain.getAxis(1).n;
    int nz2 = _range.getAxis(1).n;
    int ny = _domain.getN123()/(nz);

    if (add == false) memset(pmod, 0, _domain.getN123()*sizeof(data_t));

    int count;
    #pragma omp parallel for private(count)
    for (int i=0; i<ny; i++){
        count = 0;
        for (int iz=0; iz<nz; iz++){
            for (int j=count; j<count+_mz[iz]; j++) pmod[i*nz+iz] += pdat[i*nz2+j];
            count += _mz[iz];
        }
    }
}

void bsplines3d::apply_forward(bool add, const data_t * pmod, data_t * pdat){

    axis<data_t> Z = _range.getAxis(1);
    axis<data_t> X = _range.getAxis(2);
    axis<data_t> Y = _range.getAxis(3);
    int nc = _range.getN123()/(Y.n*X.n*Z.n);
    int kx = _kx.size()-4;
    int ky = _ky.size()-4;
    int kz = _kz.size()-4;
    data_t x, y, z;

    const data_t (* pm) [ky][kx][kz] = (const data_t (*) [ky][kx][kz]) pmod; 
    data_t (* __restrict pd) [Y.n][X.n][Z.n] = (data_t (*) [Y.n][X.n][Z.n]) pdat;

    if (!add) memset(pdat, 0, _range.getN123()*sizeof(data_t));

    for (int ic=0; ic<nc; ic++){
        #pragma omp parallel for
        for (int iy=0; iy<Y.n; iy++){
            for (int ix=0; ix<X.n; ix++){
                for (int iz=0; iz<Z.n; iz++){
                    for (int k=0; k<4; k++){
                        for (int i=0; i<4; i++){
                            for (int j=0; j<4; j++){
                                pd[ic][iy][ix][iz] += _N3y[k][iy]*_N3x[i][ix]*_N3z[j][iz]*pm[ic][k+_kymin[iy]][i+_kxmin[ix]][j+_kzmin[iz]];
                            }
                        }
                    }
                }
            }
        }
    }
}

void bsplines3d::apply_adjoint(bool add, data_t * pmod, const data_t * pdat){
    
    axis<data_t> Z = _range.getAxis(1);
    axis<data_t> X = _range.getAxis(2);
    axis<data_t> Y = _range.getAxis(3);
    int nc = _range.getN123()/(Y.n*X.n*Z.n);
    int kx = _kx.size()-4;
    int ky = _ky.size()-4;
    int kz = _kz.size()-4;
    data_t x, y, z;

    data_t (* __restrict pm) [ky][kx][kz] = (data_t (*) [ky][kx][kz]) pmod; 
    const data_t (* pd) [Y.n][X.n][Z.n] = (const data_t (*) [Y.n][X.n][Z.n]) pdat;

    if (!add) memset(pmod, 0, _domain.getN123()*sizeof(data_t));

    #pragma omp parallel for
    for (int ic=0; ic<nc; ic++){
        for (int iy=0; iy<Y.n; iy++){
            for (int ix=0; ix<X.n; ix++){
                for (int iz=0; iz<Z.n; iz++){
                    for (int k=0; k<4; k++){
                        for (int i=0; i<4; i++){
                            for (int j=0; j<4; j++){
                                pm[ic][k+_kymin[iy]][i+_kxmin[ix]][j+_kzmin[iz]] += _N3y[k][iy]*_N3x[i][ix]*_N3z[j][iz]*pd[ic][iy][ix][iz];
                            }
                        }
                    }
                }
            }
        }
    }
}

void bsplines2d::apply_forward(bool add, const data_t * pmod, data_t * pdat){

    axis<data_t> Z = _range.getAxis(1);
    axis<data_t> X = _range.getAxis(2);
    int ny = _range.getN123()/(X.n*Z.n);
    int kx = _kx.size()-4;
    int kz = _kz.size()-4;
    data_t x, z;

    if (add == false) memset(pdat, 0, _range.getN123()*sizeof(data_t));

    for (int iy=0; iy<ny; iy++){
        #pragma omp parallel for private(x,z)
        for (int ix=0; ix<X.n; ix++){
            x = ix*X.d+X.o;
            for (int iz=0; iz<Z.n; iz++){
                z = iz*Z.d+Z.o;
                //pdat[iy*X.n*Z.n+ix*Z.n+iz] = add*pdat[iy*X.n*Z.n+ix*Z.n+iz];
                for (int i=_kxmin[ix]; i<_kxmin[ix]+4; i++){
                    for (int j=_kzmin[iz]; j<_kzmin[iz]+4; j++){
                        pdat[iy*X.n*Z.n+ix*Z.n+iz] += N3(i,x,_kx)*N3(j,z,_kz)*pmod[iy*kx*kz+i*kz+j];
                    }
                }
            }
        }
    }
}
void bsplines2d::apply_adjoint(bool add, data_t * pmod, const data_t * pdat){
    
    axis<data_t> Z = _range.getAxis(1);
    axis<data_t> X = _range.getAxis(2);
    int ny = _range.getN123()/(X.n*Z.n);
    int kx = _kx.size()-4;
    int kz = _kz.size()-4;
    data_t x, z;

    if (add == false) memset(pmod, 0, _domain.getN123()*sizeof(data_t));

    for (int iy=0; iy<ny; iy++){
        //#pragma omp parallel for private(x,z)
        for (int ix=0; ix<X.n; ix++){
            x = ix*X.d+X.o;
            for (int iz=0; iz<Z.n; iz++){
                z = iz*Z.d+Z.o;
                for (int i=_kxmin[ix]; i<_kxmin[ix]+4; i++){
                    for (int j=_kzmin[iz]; j<_kzmin[iz]+4; j++){
                        pmod[iy*kx*kz+i*kz+j] += N3(i,x,_kx)*N3(j,z,_kz)*pdat[iy*X.n*Z.n+ix*Z.n+iz];
                    }
                }
            }
        }
    }
}

void bsplines1d::apply_forward(bool add, const data_t * pmod, data_t * pdat){

    axis<data_t> Z = _range.getAxis(1);
    int ny = _range.getN123()/(Z.n);
    int kz = _kz.size()-4;
    data_t z;

    if (add == false) memset(pdat, 0, _range.getN123()*sizeof(data_t));

    #pragma omp parallel for private(z)
    for (int iy=0; iy<ny; iy++){
        for (int iz=0; iz<Z.n; iz++){
            z = iz*Z.d+Z.o;
            for (int j=_kzmin[iz]; j<_kzmin[iz]+4; j++){
                pdat[iy*Z.n+iz] += N3(j,z,_kz)*pmod[iy*kz+j];
            }
        }
    }
}
void bsplines1d::apply_adjoint(bool add, data_t * pmod, const data_t * pdat){
    
    axis<data_t> Z = _range.getAxis(1);
    int ny = _range.getN123()/(Z.n);
    int kz = _kz.size()-4;
    data_t z;

    if (add == false) memset(pmod, 0, _domain.getN123()*sizeof(data_t));

    for (int iy=0; iy<ny; iy++){
        for (int iz=0; iz<Z.n; iz++){
            z = iz*Z.d+Z.o;
            for (int j=_kzmin[iz]; j<_kzmin[iz]+4; j++){
                pmod[iy*kz+j] += N3(j,z,_kz)*pdat[iy*Z.n+iz];
            }
        }
    }
}

void bsfillin3d::apply_forward(bool add, const data_t * pmod, data_t * pdat){

    axis<data_t> Z = _range.getAxis(1);
    axis<data_t> Z2 = _domain.getAxis(1);
    axis<data_t> X = _range.getAxis(2);
    axis<data_t> X2 = _domain.getAxis(2);
    axis<data_t> Y = _range.getAxis(3);
    axis<data_t> Y2 = _domain.getAxis(3);
    int nc = _range.getN123()/(Y.n*X.n*Z.n);
    data_t (* __restrict pc) [Y.n][X.n][Z.n] = (data_t (*) [Y.n][X.n][Z.n]) pdat;
    const data_t (*pv) [Y2.n][X2.n][Z2.n] = (const data_t (*) [Y2.n][X2.n][Z2.n]) pmod;
    for (int ic=0; ic<nc; ic++){
        #pragma omp parallel for
        for (int iy=0; iy<Y.n; iy++){
            int iy2 = floor((_controly[iy]-Y2.o)/Y2.d);
            int iy3 = std::min(Y2.n-1,iy2+1);
            data_t wy = (_controly[iy] - Y2.o - iy2*Y2.d)/Y2.d;
            for (int ix=0; ix<X.n; ix++){
                int ix2 = floor((_controlx[ix]-X2.o)/X2.d);
                int ix3 = std::min(X2.n-1,ix2+1);
                data_t wx = (_controlx[ix] - X2.o - ix2*X2.d)/X2.d;
                for (int iz=0; iz<Z.n; iz++){
                    int iz2 = floor((_controlz[iz]-Z2.o)/Z2.d);
                    int iz3 = std::min(Z2.n-1,iz2+1);
                    data_t wz = (_controlz[iz] - Z2.o - iz2*Z2.d)/Z2.d;
                    pc[ic][iy][ix][iz] = add*pc[ic][iy][ix][iz] + (1-wx)*(1-wy)*(1-wz)*pv[ic][iy2][ix2][iz2]  
                    + (wx)*(1-wy)*(1-wz)*pv[ic][iy2][ix3][iz2] 
                    + (1-wx)*(wy)*(1-wz)*pv[ic][iy3][ix2][iz2]
                    + (1-wx)*(1-wy)*(wz)*pv[ic][iy2][ix2][iz3]
                    + (wx)*(wy)*(1-wz)*pv[ic][iy3][ix3][iz2]
                    + (wx)*(1-wy)*(wz)*pv[ic][iy2][ix3][iz3]
                    + (1-wx)*(wy)*(wz)*pv[ic][iy3][ix2][iz3]
                    + (wx)*(wy)*(wz)*pv[ic][iy3][ix3][iz3];
                }
            }
        }
    }
}
void bsfillin3d::apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

    if (!add) memset(pmod, 0, _domain.getN123()*sizeof(data_t));
    
    axis<data_t> Z = _range.getAxis(1);
    axis<data_t> Z2 = _domain.getAxis(1);
    axis<data_t> X = _range.getAxis(2);
    axis<data_t> X2 = _domain.getAxis(2);
    axis<data_t> Y = _range.getAxis(3);
    axis<data_t> Y2 = _domain.getAxis(3);
    int nc = _range.getN123()/(Y.n*X.n*Z.n);
    const data_t (*pc) [Y.n][X.n][Z.n] = (const data_t (*) [Y.n][X.n][Z.n]) pdat;
    data_t (* __restrict pv) [Y2.n][X2.n][Z2.n] = (data_t (*) [Y2.n][X2.n][Z2.n]) pmod;

    #pragma omp parallel for
    for (int ic=0; ic<nc; ic++){
        for (int iy=0; iy<Y.n; iy++){
            int iy2 = floor((_controly[iy]-Y2.o)/Y2.d);
            int iy3 = std::min(Y2.n-1,iy2+1);
            data_t wy = (_controly[iy] - Y2.o - iy2*Y2.d)/Y2.d;
            for (int ix=0; ix<X.n; ix++){
                int ix2 = floor((_controlx[ix]-X2.o)/X2.d);
                int ix3 = std::min(X2.n-1,ix2+1);
                data_t wx = (_controlx[ix] - X2.o - ix2*X2.d)/X2.d;
                for (int iz=0; iz<Z.n; iz++){
                    int iz2 = floor((_controlz[iz]-Z2.o)/Z2.d);
                    int iz3 = std::min(Z2.n-1,iz2+1);
                    data_t wz = (_controlz[iz] - Z2.o - iz2*Z2.d)/Z2.d;

                    pv[ic][iy2][ix2][iz2] += (1-wx)*(1-wy)*(1-wz)*pc[ic][iy][ix][iz];
                    pv[ic][iy2][ix3][iz2] += (wx)*(1-wy)*(1-wz)*pc[ic][iy][ix][iz];
                    pv[ic][iy3][ix2][iz2] += (1-wx)*(wy)*(1-wz)*pc[ic][iy][ix][iz];
                    pv[ic][iy2][ix2][iz3] += (1-wx)*(1-wy)*(wz)*pc[ic][iy][ix][iz];
                    pv[ic][iy3][ix3][iz2] += (wx)*(wy)*(1-wz)*pc[ic][iy][ix][iz];
                    pv[ic][iy2][ix3][iz3] += (wx)*(1-wy)*(wz)*pc[ic][iy][ix][iz];
                    pv[ic][iy3][ix2][iz3] += (1-wx)*(wy)*(wz)*pc[ic][iy][ix][iz];
                    pv[ic][iy3][ix3][iz3] += (wx)*(wy)*(wz)*pc[ic][iy][ix][iz];
                }
            }
        }
    }
}

#undef ZERO