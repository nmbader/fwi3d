#include <string.h>
#include <algorithm>
#include <omp.h>
#include "spatial_operators.hpp"

#define IZ(iz) (iy*nx*nz+ix*nz+iz)
#define IX(ix) (iy*nx*nz+(ix)*nz+iz)
#define IY(iy) ((iy)*nx*nz+ix*nz+iz)
#define IXYZ(ix,iy,iz) ((iy)*nx*nz+(ix)*nz+iz)

void applyHz(bool inv, bool add, const data_t* in, data_t* out, int nx, int ny, int nz, data_t d, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax){

    int izminb=std::max(4,izmin);
    int izmaxb=std::min(nz - 4, izmax);
    data_t coef[4] = {17.0/48, 59.0/48, 43.0/48, 49.0/48};
    const data_t (* pin) [nx][nz] = (const data_t (*) [nx][nz]) in; 
    data_t (* __restrict pout) [nx][nz] = (data_t (*) [nx][nz]) out; 

    if (!inv){
        #pragma omp parallel for
        for (int iy=iymin; iy<iymax; iy++){
            for (int ix=ixmin; ix<ixmax; ix++){
                for (int iz = izmin; iz<4; iz++){
                    pout[iy][ix][iz] = add*pout[iy][ix][iz] + pin[iy][ix][iz]*coef[iz]*d;
                }
                for (int iz = nz-izmax; iz<4; iz++){
                    pout[iy][ix][nz-1-iz] = add*pout[iy][ix][nz-1-iz] + pin[iy][ix][nz-1-iz]*coef[iz]*d;
                }
                #pragma ivdep
                for (int iz=izminb; iz<izmaxb; iz++){
                    pout[iy][ix][iz] = add*pout[iy][ix][iz] + pin[iy][ix][iz]*d;
                }
            }
        }
    }
    else{
        #pragma omp parallel for
        for (int iy=iymin; iy<iymax; iy++){
            for (int ix=ixmin; ix<ixmax; ix++){
                for (int iz = izmin; iz<4; iz++){
                    pout[iy][ix][iz] = add*pout[iy][ix][iz] + pin[iy][ix][iz]/(coef[iz]*d);
                }
                for (int iz = nz-izmax; iz<4; iz++){
                    pout[iy][ix][nz-1-iz] = add*pout[iy][ix][nz-1-iz] + pin[iy][ix][nz-1-iz]/(coef[iz]*d);
                }
                #pragma ivdep
                for (int iz=izminb; iz<izmaxb; iz++){
                    pout[iy][ix][iz] = add*pout[iy][ix][iz] + pin[iy][ix][iz]/d;
                }
            }
        }
    }
}

void applyHx(bool inv, bool add, const data_t* in, data_t* out, int nx, int ny, int nz, data_t d, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax){

    int ixminb=std::max(4,ixmin);
    int ixmaxb=std::min(nx - 4, ixmax);
    data_t coef[4] = {17.0/48, 59.0/48, 43.0/48, 49.0/48};
    const data_t (* pin) [nx][nz] = (const data_t (*) [nx][nz]) in; 
    data_t (* __restrict pout) [nx][nz] = (data_t (*) [nx][nz]) out; 

    if (!inv){
        #pragma omp parallel for
        for (int iy=iymin; iy<iymax; iy++){
            for (int ix=ixmin; ix<4; ix++){
                #pragma ivdep
                for (int iz=izmin; iz<izmax; iz++){
                    pout[iy][ix][iz] = add*pout[iy][ix][iz] + pin[iy][ix][iz]*coef[ix]*d;
                }
            }
            for (int ix=nx-ixmax; ix<4; ix++){
                #pragma ivdep
                for (int iz=izmin; iz<izmax; iz++){
                    pout[iy][nx-1-ix][iz] = add*pout[iy][nx-1-ix][iz] + pin[iy][nx-1-ix][iz]*coef[ix]*d;
                }
            }
            for (int ix=ixminb; ix<ixmaxb; ix++){
                #pragma ivdep
                for (int iz=izmin; iz<izmax; iz++){
                    pout[iy][ix][iz] = add*pout[iy][ix][iz] + pin[iy][ix][iz]*d;
                }
            }
        }
    }
    else{
        #pragma omp parallel for
        for (int iy=iymin; iy<iymax; iy++){
            for (int ix=ixmin; ix<4; ix++){
                #pragma ivdep
                for (int iz=izmin; iz<izmax; iz++){
                    pout[iy][ix][iz] = add*pout[iy][ix][iz] + pin[iy][ix][iz]/(coef[ix]*d);
                }
            }
            for (int ix=nx-ixmax; ix<4; ix++){
                #pragma ivdep
                for (int iz=izmin; iz<izmax; iz++){
                    pout[iy][nx-1-ix][iz] = add*pout[iy][nx-1-ix][iz] + pin[iy][nx-1-ix][iz]/(coef[ix]*d);
                }
            }
            for (int ix=ixminb; ix<ixmaxb; ix++){
                #pragma ivdep
                for (int iz=izmin; iz<izmax; iz++){
                    pout[iy][ix][iz] = add*pout[iy][ix][iz] + pin[iy][ix][iz]/d;
                }
            }
        }
    }
}

void applyHy(bool inv, bool add, const data_t* in, data_t* out, int nx, int ny, int nz, data_t d, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax){

    int iyminb=std::max(4,iymin);
    int iymaxb=std::min(ny - 4, iymax);
    data_t coef[4] = {17.0/48, 59.0/48, 43.0/48, 49.0/48};
    const data_t (* pin) [nx][nz] = (const data_t (*) [nx][nz]) in; 
    data_t (* __restrict pout) [nx][nz] = (data_t (*) [nx][nz]) out; 

    if (!inv){
        #pragma omp parallel for
        for (int iy=iymin; iy<4; iy++){
            for (int ix=ixmin; ix<ixmax; ix++){
                #pragma ivdep
                for (int iz=izmin; iz<izmax; iz++){
                    pout[iy][ix][iz] = add*pout[iy][ix][iz] + pin[iy][ix][iz]*coef[iy]*d;
                }
            }
        }
        #pragma omp parallel for
        for (int iy=ny-iymax; iy<4; iy++){
            for (int ix=ixmin; ix<ixmax; ix++){
                #pragma ivdep
                for (int iz=izmin; iz<izmax; iz++){
                    pout[ny-1-iy][ix][iz] = add*pout[ny-1-iy][ix][iz] + pin[ny-1-iy][ix][iz]*coef[iy]*d;
                }
            }
        }
        #pragma omp parallel for
        for (int iy=iyminb; iy<iymaxb; iy++){
            for (int ix=ixmin; ix<ixmax; ix++){
                #pragma ivdep
                for (int iz=izmin; iz<izmax; iz++){
                    pout[iy][ix][iz] = add*pout[iy][ix][iz] + pin[iy][ix][iz]*d;
                }
            }
        }
    }
    else{
        #pragma omp parallel for
        for (int iy=iymin; iy<4; iy++){
            for (int ix=ixmin; ix<ixmax; ix++){
                #pragma ivdep
                for (int iz=izmin; iz<izmax; iz++){
                    pout[iy][ix][iz] = add*pout[iy][ix][iz] + pin[iy][ix][iz]/(coef[iy]*d);
                }
            }
        }
        #pragma omp parallel for
        for (int iy=ny-iymax; iy<4; iy++){
            for (int ix=ixmin; ix<ixmax; ix++){
                #pragma ivdep
                for (int iz=izmin; iz<izmax; iz++){
                    pout[ny-1-iy][ix][iz] = add*pout[ny-1-iy][ix][iz] + pin[ny-1-iy][ix][iz]/(coef[iy]*d);
                }
            }
        }
        #pragma omp parallel for
        for (int iy=iyminb; iy<iymaxb; iy++){
            for (int ix=ixmin; ix<ixmax; ix++){
                #pragma ivdep
                for (int iz=izmin; iz<izmax; iz++){
                    pout[iy][ix][iz] = add*pout[iy][ix][iz] + pin[iy][ix][iz]/d;
                }
            }
        }
    }
}

void Dz(bool add, const data_t* in, data_t* out, int nx, int ny, int nz, data_t d, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax){
    
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t coef0=2.0/3;
    data_t coef1=-1.0/12;
    int nc1=4, nc2=6;
    int izminb=std::max(nc1,izmin);
    int izmaxb=std::min(nz - nc1, izmax);
    const data_t (* pin) [nx][nz] = (const data_t (*) [nx][nz]) in; 
    data_t (* __restrict pout) [nx][nz] = (data_t (*) [nx][nz]) out;

    #pragma omp parallel for
    for (int iy=iymin; iy<iymax; iy++){
        for (int ix = ixmin; ix < ixmax; ix++){
            // apply the operator near the top boundary if included
            for (int iz=izmin; iz<nc1; iz++){
                pout[iy][ix][iz] = add*pout[iy][ix][iz] + (bnd_coef[iz*nc2] * pin[iy][ix][0] + bnd_coef[iz*nc2+1] * pin[iy][ix][1] + bnd_coef[iz*nc2+2] * pin[iy][ix][2] + bnd_coef[iz*nc2+3] * pin[iy][ix][3] + bnd_coef[iz*nc2+4] * pin[iy][ix][4] + bnd_coef[iz*nc2+5] * pin[iy][ix][5]) / d;
            }

            // apply the operator near the bottom boundary if included
            for (int iz=nz-izmax; iz<nc1; iz++){
                pout[iy][ix][nz-1-iz] = add*pout[iy][ix][nz-1-iz] - (bnd_coef[iz*nc2] * pin[iy][ix][nz-1] + bnd_coef[iz*nc2+1] * pin[iy][ix][nz-2] + bnd_coef[iz*nc2+2] * pin[iy][ix][nz-3] + bnd_coef[iz*nc2+3] * pin[iy][ix][nz-4] + bnd_coef[iz*nc2+4] * pin[iy][ix][nz-5] + bnd_coef[iz*nc2+5] * pin[iy][ix][nz-6]) / d;
            }

            // apply the operator in the interior
            #pragma omp simd
            for (int iz=izminb; iz<izmaxb; iz++){
                pout[iy][ix][iz] = add*pout[iy][ix][iz] + (coef0 * (pin[iy][ix][iz+1] - pin[iy][ix][iz-1]) + coef1 * (pin[iy][ix][iz+2] - pin[iy][ix][iz-2]))/d;
            }
        }
    }
}

void Dx(bool add, const data_t* in, data_t* out, int nx, int ny, int nz, data_t d, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax){
    
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t coef0=2.0/3;
    data_t coef1=-1.0/12;
    int nc1=4, nc2=6;
    int ixminb=std::max(nc1,ixmin);
    int ixmaxb=std::min(nx - nc1, ixmax);
    const data_t (* pin) [nx][nz] = (const data_t (*) [nx][nz]) in; 
    data_t (* __restrict pout) [nx][nz] = (data_t (*) [nx][nz]) out;

    #pragma omp parallel for
    for (int iy=iymin; iy<iymax; iy++){
        // apply the operator near the left boundary if included
        for (int ix=ixmin; ix<nc1; ix++){
            #pragma omp simd
            for (int iz = izmin; iz < izmax; iz++){
                pout[iy][ix][iz] = add*pout[iy][ix][iz] + (bnd_coef[ix*nc2] * pin[iy][0][iz] + bnd_coef[ix*nc2+1] * pin[iy][1][iz] + bnd_coef[ix*nc2+2] * pin[iy][2][iz] + bnd_coef[ix*nc2+3] * pin[iy][3][iz] + bnd_coef[ix*nc2+4] * pin[iy][4][iz] + bnd_coef[ix*nc2+5] * pin[iy][5][iz]) / d;
            }
        }
        // apply the operator near the right boundary if included
        for (int ix=nx-ixmax; ix<nc1; ix++){
            #pragma omp simd
            for (int iz = izmin; iz < izmax; iz++){
                pout[iy][nx-1-ix][iz] = add*pout[iy][nx-1-ix][iz] - (bnd_coef[ix*nc2] * pin[iy][nx-1][iz] + bnd_coef[ix*nc2+1] * pin[iy][nx-2][iz] + bnd_coef[ix*nc2+2] * pin[iy][nx-3][iz] + bnd_coef[ix*nc2+3] * pin[iy][nx-4][iz] + bnd_coef[ix*nc2+4] * pin[iy][nx-5][iz] + bnd_coef[ix*nc2+5] * pin[iy][nx-6][iz]) / d;
            }
        }
        // apply the operator in the interior
        for (int ix=ixminb; ix<ixmaxb; ix++){
            #pragma omp simd
            for (int iz = izmin; iz < izmax; iz++){
                pout[iy][ix][iz] = add*pout[iy][ix][iz] + (coef0 * (pin[iy][ix+1][iz] - pin[iy][ix-1][iz]) + coef1 * (pin[iy][ix+2][iz] - pin[iy][ix-2][iz]))/d;
            }
        }
    }
}

void Dy(bool add, const data_t* in, data_t* out, int nx, int ny, int nz, data_t d, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax){
    
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t coef0=2.0/3;
    data_t coef1=-1.0/12;
    int nc1=4, nc2=6;
    int iyminb=std::max(nc1,iymin);
    int iymaxb=std::min(ny - nc1, iymax);
    const data_t (* pin) [nx][nz] = (const data_t (*) [nx][nz]) in; 
    data_t (* __restrict pout) [nx][nz] = (data_t (*) [nx][nz]) out;

    // apply the operator near the front boundary if included
    #pragma omp parallel for
    for (int iy=iymin; iy<nc1; iy++){
        for (int ix=ixmin; ix<ixmax; ix++){
            #pragma omp simd
            for (int iz=izmin; iz<izmax; iz++){
                pout[iy][ix][iz] = add*pout[iy][ix][iz] + (bnd_coef[iy*nc2] * pin[0][ix][iz] + bnd_coef[iy*nc2+1] * pin[1][ix][iz] + bnd_coef[iy*nc2+2] * pin[2][ix][iz] + bnd_coef[iy*nc2+3] * pin[3][ix][iz] + bnd_coef[iy*nc2+4] * pin[4][ix][iz] + bnd_coef[iy*nc2+5] * pin[5][ix][iz]) / d;
            }
        }
    }
    // apply the operator near the back boundary if included
    #pragma omp parallel for
    for (int iy=ny-iymax; iy<nc1; iy++){
        for (int ix=ixmin; ix<ixmax; ix++){
            #pragma omp simd
            for (int iz=izmin; iz<izmax; iz++){
                pout[ny-1-iy][ix][iz] = add*pout[ny-1-iy][ix][iz] - (bnd_coef[iy*nc2] * pin[ny-1][ix][iz] + bnd_coef[iy*nc2+1] * pin[ny-2][ix][iz] + bnd_coef[iy*nc2+2] * pin[ny-3][ix][iz] + bnd_coef[iy*nc2+3] * pin[ny-4][ix][iz] + bnd_coef[iy*nc2+4] * pin[ny-5][ix][iz] + bnd_coef[iy*nc2+5] * pin[ny-6][ix][iz]) / d;
            }
        }
    }
    // apply the operator in the interior
    #pragma omp parallel for
    for (int iy=iyminb; iy<iymaxb; iy++){
        for (int ix=ixmin; ix<ixmax; ix++){
            #pragma omp simd
            for (int iz=izmin; iz<izmax; iz++){
                pout[iy][ix][iz] = add*pout[iy][ix][iz] + (coef0 * (pin[iy+1][ix][iz] - pin[iy-1][ix][iz]) + coef1 * (pin[iy+2][ix][iz] - pin[iy-2][ix][iz]))/d;
            }
        }
    }
}

void mult_Dz(bool add, const data_t* in, data_t* out, int nx, int ny, int nz, data_t d, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax, const data_t * par){
    
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t coef0=2.0/3;
    data_t coef1=-1.0/12;
    int nc1=4, nc2=6;
    int izminb=std::max(nc1,izmin);
    int izmaxb=std::min(nz - nc1, izmax);
    const data_t (* pp) [nx][nz] = (const data_t (*) [nx][nz]) par; 
    const data_t (* pin) [nx][nz] = (const data_t (*) [nx][nz]) in; 
    data_t (* __restrict pout) [nx][nz] = (data_t (*) [nx][nz]) out;

    #pragma omp parallel for
    for (int iy=iymin; iy<iymax; iy++){
        for (int ix = ixmin; ix < ixmax; ix++){
            // apply the operator near the top boundary if included
            for (int iz=izmin; iz<nc1; iz++){
                pout[iy][ix][iz] = add*pout[iy][ix][iz] + (bnd_coef[iz*nc2] * pin[iy][ix][0]*pp[iy][ix][0] + bnd_coef[iz*nc2+1] * pin[iy][ix][1]*pp[iy][ix][1] + bnd_coef[iz*nc2+2] * pin[iy][ix][2]*pp[iy][ix][2] + bnd_coef[iz*nc2+3] * pin[iy][ix][3]*pp[iy][ix][3] + bnd_coef[iz*nc2+4] * pin[iy][ix][4]*pp[iy][ix][4] + bnd_coef[iz*nc2+5] * pin[iy][ix][5]*pp[iy][ix][5])/d;
            }

            // apply the operator near the bottom boundary if included
            for (int iz=nz-izmax; iz<nc1; iz++){
                pout[iy][ix][nz-1-iz] = add*pout[iy][ix][nz-1-iz] - (bnd_coef[iz*nc2] * pin[iy][ix][nz-1]*pp[iy][ix][nz-1] + bnd_coef[iz*nc2+1] * pin[iy][ix][nz-2]*pp[iy][ix][nz-2] + bnd_coef[iz*nc2+2] * pin[iy][ix][nz-3]*pp[iy][ix][nz-3] + bnd_coef[iz*nc2+3] * pin[iy][ix][nz-4]*pp[iy][ix][nz-4] + bnd_coef[iz*nc2+4] * pin[iy][ix][nz-5]*pp[iy][ix][nz-5] + bnd_coef[iz*nc2+5] * pin[iy][ix][nz-6]*pp[iy][ix][nz-6])/d;
            }

            // apply the operator in the interior
            #pragma omp simd
            for (int iz=izminb; iz<izmaxb; iz++){
                pout[iy][ix][iz] = add*pout[iy][ix][iz] + (coef0 * (pin[iy][ix][iz+1]*pp[iy][ix][iz+1] - pin[iy][ix][iz-1]*pp[iy][ix][iz-1]) + coef1 * (pin[iy][ix][iz+2]*pp[iy][ix][iz+2] - pin[iy][ix][iz-2]*pp[iy][ix][iz-2]))/d;
            }
        }
    }
}

void mult_Dx(bool add, const data_t* in, data_t* out, int nx, int ny, int nz, data_t d, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax, const data_t * par){
    
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t coef0=2.0/3;
    data_t coef1=-1.0/12;
    int nc1=4, nc2=6;
    int ixminb=std::max(nc1,ixmin);
    int ixmaxb=std::min(nx - nc1, ixmax);
    const data_t (* pp) [nx][nz] = (const data_t (*) [nx][nz]) par; 
    const data_t (* pin) [nx][nz] = (const data_t (*) [nx][nz]) in; 
    data_t (* __restrict pout) [nx][nz] = (data_t (*) [nx][nz]) out;

    #pragma omp parallel for
    for (int iy=iymin; iy<iymax; iy++){
        // apply the operator near the left boundary if included
        for (int ix=ixmin; ix<nc1; ix++){
            #pragma omp simd
            for (int iz = izmin; iz < izmax; iz++){
                pout[iy][ix][iz] = add*pout[iy][ix][iz] + (bnd_coef[ix*nc2] * pin[iy][0][iz]*pp[iy][0][iz] + bnd_coef[ix*nc2+1] * pin[iy][1][iz]*pp[iy][1][iz] + bnd_coef[ix*nc2+2] * pin[iy][2][iz]*pp[iy][2][iz] + bnd_coef[ix*nc2+3] * pin[iy][3][iz]*pp[iy][3][iz] + bnd_coef[ix*nc2+4] * pin[iy][4][iz]*pp[iy][4][iz] + bnd_coef[ix*nc2+5] * pin[iy][5][iz]*pp[iy][5][iz])/d;
            }
        }
        // apply the operator near the right boundary if included
        for (int ix=nx-ixmax; ix<nc1; ix++){
            #pragma omp simd
            for (int iz = izmin; iz < izmax; iz++){
                pout[iy][nx-1-ix][iz] = add*pout[iy][nx-1-ix][iz] - (bnd_coef[ix*nc2] * pin[iy][nx-1][iz]*pp[iy][nx-1][iz] + bnd_coef[ix*nc2+1] * pin[iy][nx-2][iz]*pp[iy][nx-2][iz] + bnd_coef[ix*nc2+2] * pin[iy][nx-3][iz]*pp[iy][nx-3][iz] + bnd_coef[ix*nc2+3] * pin[iy][nx-4][iz]*pp[iy][nx-4][iz] + bnd_coef[ix*nc2+4] * pin[iy][nx-5][iz]*pp[iy][nx-5][iz] + bnd_coef[ix*nc2+5] * pin[iy][nx-6][iz]*pp[iy][nx-6][iz])/d;
            }
        }
        // apply the operator in the interior
        for (int ix=ixminb; ix<ixmaxb; ix++){
            #pragma omp simd
            for (int iz = izmin; iz < izmax; iz++){
                pout[iy][ix][iz] = add*pout[iy][ix][iz] + (coef0 * (pin[iy][ix+1][iz]*pp[iy][ix+1][iz] - pin[iy][ix-1][iz]*pp[iy][ix-1][iz]) + coef1 * (pin[iy][ix+2][iz]*pp[iy][ix+2][iz] - pin[iy][ix-2][iz]*pp[iy][ix-2][iz]))/d;
            }
        }
    }
}

void mult_Dy(bool add, const data_t* in, data_t* out, int nx, int ny, int nz, data_t d, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax, const data_t * par){
    
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t coef0=2.0/3;
    data_t coef1=-1.0/12;
    int nc1=4, nc2=6;
    int iyminb=std::max(nc1,iymin);
    int iymaxb=std::min(ny - nc1, iymax);
    const data_t (* pp) [nx][nz] = (const data_t (*) [nx][nz]) par; 
    const data_t (* pin) [nx][nz] = (const data_t (*) [nx][nz]) in; 
    data_t (* __restrict pout) [nx][nz] = (data_t (*) [nx][nz]) out;

    // apply the operator near the front boundary if included
    #pragma omp parallel for
    for (int iy=iymin; iy<nc1; iy++){
        for (int ix=ixmin; ix<ixmax; ix++){
            #pragma omp simd
            for (int iz=izmin; iz<izmax; iz++){
                pout[iy][ix][iz] = add*pout[iy][ix][iz] + (bnd_coef[iy*nc2] * pin[0][ix][iz]*pp[0][ix][iz] + bnd_coef[iy*nc2+1] * pin[1][ix][iz]*pp[1][ix][iz] + bnd_coef[iy*nc2+2] * pin[2][ix][iz]*pp[2][ix][iz] + bnd_coef[iy*nc2+3] * pin[3][ix][iz]*pp[3][ix][iz] + bnd_coef[iy*nc2+4] * pin[4][ix][iz]*pp[4][ix][iz] + bnd_coef[iy*nc2+5] * pin[5][ix][iz]*pp[5][ix][iz])/d;
            }
        }
    }
    // apply the operator near the back boundary if included
    #pragma omp parallel for
    for (int iy=ny-iymax; iy<nc1; iy++){
        for (int ix=ixmin; ix<ixmax; ix++){
            #pragma omp simd
            for (int iz=izmin; iz<izmax; iz++){
                pout[ny-1-iy][ix][iz] = add*pout[ny-1-iy][ix][iz] - (bnd_coef[iy*nc2] * pin[ny-1][ix][iz]*pp[ny-1][ix][iz] + bnd_coef[iy*nc2+1] * pin[ny-2][ix][iz]*pp[ny-2][ix][iz] + bnd_coef[iy*nc2+2] * pin[ny-3][ix][iz]*pp[ny-3][ix][iz] + bnd_coef[iy*nc2+3] * pin[ny-4][ix][iz]*pp[ny-4][ix][iz] + bnd_coef[iy*nc2+4] * pin[ny-5][ix][iz]*pp[ny-5][ix][iz] + bnd_coef[iy*nc2+5] * pin[ny-6][ix][iz]*pp[ny-6][ix][iz])/d;
            }
        }
    }
    // apply the operator in the interior
    #pragma omp parallel for
    for (int iy=iyminb; iy<iymaxb; iy++){
        for (int ix=ixmin; ix<ixmax; ix++){
            #pragma omp simd
            for (int iz=izmin; iz<izmax; iz++){
                pout[iy][ix][iz] = add*pout[iy][ix][iz] + (coef0 * (pin[iy+1][ix][iz]*pp[iy+1][ix][iz] - pin[iy-1][ix][iz]*pp[iy-1][ix][iz]) + coef1 * (pin[iy+2][ix][iz]*pp[iy+2][ix][iz] - pin[iy-2][ix][iz]*pp[iy-2][ix][iz]))/d;
            }
        }
    }
}

void esat_scale_boundaries(data_t** in, int nx, int ny, int nz, data_t dx, data_t dy, data_t dz, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax, const data_t** par, data_t dt, bool top, bool bottom, bool left, bool right, bool front, bool back){
    
    // par must at least contain lambda, mu, rho in that order
    data_t h0 = 17.0/48;
    data_t scalerx, scalery, scalerz; // scalers for faces without the corners
    std::vector<data_t> sc(24,0); // scalers for the corners: top front left,right, top back right, left, bottom front left, right, bottom back right, left
    data_t (* psc) [3] = (data_t (*) [3]) sc.data();

    // scale the top boundary without the corners
    if (top){
        psc[0][0] += sqrt(par[1][0]/par[2][0])*dt / (2  * dz * h0); // top front left
        psc[0][1] += sqrt(par[1][0]/par[2][0])*dt / (2  * dz * h0);
        psc[0][2] += sqrt((par[0][0]+2*par[1][0])/par[2][0])*dt / (2  * dz * h0);
        psc[1][0] += sqrt(par[1][IXYZ(nx-1,0,0)]/par[2][IXYZ(nx-1,0,0)])*dt / (2  * dz * h0); // top front left
        psc[1][1] += sqrt(par[1][IXYZ(nx-1,0,0)]/par[2][IXYZ(nx-1,0,0)])*dt / (2  * dz * h0);
        psc[1][2] += sqrt((par[0][IXYZ(nx-1,0,0)]+2*par[1][IXYZ(nx-1,0,0)])/par[2][IXYZ(nx-1,0,0)])*dt / (2  * dz * h0);
        psc[2][0] += sqrt(par[1][IXYZ(nx-1,ny-1,0)]/par[2][IXYZ(nx-1,ny-1,0)])*dt / (2  * dz * h0); // top back right
        psc[2][1] += sqrt(par[1][IXYZ(nx-1,ny-1,0)]/par[2][IXYZ(nx-1,ny-1,0)])*dt / (2  * dz * h0);
        psc[2][2] += sqrt((par[0][IXYZ(nx-1,ny-1,0)]+2*par[1][IXYZ(nx-1,ny-1,0)])/par[2][IXYZ(nx-1,ny-1,0)])*dt / (2  * dz * h0);
        psc[3][0] += sqrt(par[1][IXYZ(0,ny-1,0)]/par[2][IXYZ(0,ny-1,0)])*dt / (2  * dz * h0); // top back left
        psc[3][1] += sqrt(par[1][IXYZ(0,ny-1,0)]/par[2][IXYZ(0,ny-1,0)])*dt / (2  * dz * h0);
        psc[3][2] += sqrt((par[0][IXYZ(0,ny-1,0)]+2*par[1][IXYZ(0,ny-1,0)])/par[2][IXYZ(0,ny-1,0)])*dt / (2  * dz * h0);
        
        if (izmin==0){
            int iz=0;
            #pragma omp parallel for private(scalerx,scalery,scalerz)
            for (int iy = std::max(1,iymin); iy < std::min(ny-1,iymax); iy++){
                for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0);
                    scalery = scalerx;
                    scalerz = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            // treat the edges without the corners
            if (left){
                int ix=0;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iy = std::max(1,iymin); iy < std::min(ny-1,iymax); iy++){
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    scalery = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    scalerz = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }
            else{
                int ix=0;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iy = std::max(1,iymin); iy < std::min(ny-1,iymax); iy++){
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0);
                    scalery = scalerx;
                    scalerz = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            if (right){
                int ix=nx-1;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iy = std::max(1,iymin); iy < std::min(ny-1,iymax); iy++){
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    scalery = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    scalerz = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }
            else{
                int ix=nx-1;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iy = std::max(1,iymin); iy < std::min(ny-1,iymax); iy++){
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0);
                    scalery = scalerx;
                    scalerz = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            if (front){
                int iy=0;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    scalery = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    scalerz = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }
            else{
                int iy=0;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0);
                    scalery = scalerx;
                    scalerz = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            if (back){
                int iy=ny-1;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    scalery = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    scalerz = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }
            else{
                int iy=ny-1;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0);
                    scalery = scalerx;
                    scalerz = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }
        }
    }

    // scale the bottom boundary without the corners
    if (bottom){
        psc[4][0] += sqrt(par[1][nz-1]/par[2][nz-1])*dt / (2  * dz * h0); // bottom front left
        psc[4][1] += sqrt(par[1][nz-1]/par[2][nz-1])*dt / (2  * dz * h0);
        psc[4][2] += sqrt((par[0][nz-1]+2*par[1][nz-1])/par[2][nz-1])*dt / (2  * dz * h0);
        psc[5][0] += sqrt(par[1][IXYZ(nx-1,0,nz-1)]/par[2][IXYZ(nx-1,0,nz-1)])*dt / (2  * dz * h0); // bottom front left
        psc[5][1] += sqrt(par[1][IXYZ(nx-1,0,nz-1)]/par[2][IXYZ(nx-1,0,nz-1)])*dt / (2  * dz * h0);
        psc[5][2] += sqrt((par[0][IXYZ(nx-1,0,nz-1)]+2*par[1][IXYZ(nx-1,0,nz-1)])/par[2][IXYZ(nx-1,0,nz-1)])*dt / (2  * dz * h0);
        psc[6][0] += sqrt(par[1][IXYZ(nx-1,ny-1,nz-1)]/par[2][IXYZ(nx-1,ny-1,nz-1)])*dt / (2  * dz * h0); // bottom back right
        psc[6][1] += sqrt(par[1][IXYZ(nx-1,ny-1,nz-1)]/par[2][IXYZ(nx-1,ny-1,nz-1)])*dt / (2  * dz * h0);
        psc[6][2] += sqrt((par[0][IXYZ(nx-1,ny-1,nz-1)]+2*par[1][IXYZ(nx-1,ny-1,nz-1)])/par[2][IXYZ(nx-1,ny-1,nz-1)])*dt / (2  * dz * h0);
        psc[7][0] += sqrt(par[1][IXYZ(0,ny-1,nz-1)]/par[2][IXYZ(0,ny-1,nz-1)])*dt / (2  * dz * h0); // bottom back left
        psc[7][1] += sqrt(par[1][IXYZ(0,ny-1,nz-1)]/par[2][IXYZ(0,ny-1,nz-1)])*dt / (2  * dz * h0);
        psc[7][2] += sqrt((par[0][IXYZ(0,ny-1,nz-1)]+2*par[1][IXYZ(0,ny-1,nz-1)])/par[2][IXYZ(0,ny-1,nz-1)])*dt / (2  * dz * h0);
        
        if (izmax==nz){
            int iz=nz-1;
            #pragma omp parallel for private(scalerx,scalery,scalerz)
            for (int iy = std::max(1,iymin); iy < std::min(ny-1,iymax); iy++){
                for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0);
                    scalery = scalerx;
                    scalerz = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            // treat the edges without the corners
            if (left){
                int ix=0;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iy = std::max(1,iymin); iy < std::min(ny-1,iymax); iy++){
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    scalery = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    scalerz = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }
            else{
                int ix=0;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iy = std::max(1,iymin); iy < std::min(ny-1,iymax); iy++){
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0);
                    scalery = scalerx;
                    scalerz = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            if (right){
                int ix=nx-1;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iy = std::max(1,iymin); iy < std::min(ny-1,iymax); iy++){
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    scalery = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    scalerz = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }
            else{
                int ix=nx-1;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iy = std::max(1,iymin); iy < std::min(ny-1,iymax); iy++){
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0);
                    scalery = scalerx;
                    scalerz = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            if (front){
                int iy=0;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    scalery = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    scalerz = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }
            else{
                int iy=0;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0);
                    scalery = scalerx;
                    scalerz = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            if (back){
                int iy=ny-1;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    scalery = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    scalerz = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }
            else{
                int iy=ny-1;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0);
                    scalery = scalerx;
                    scalerz = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dz * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }
        }
    }

    
    // scale the left boundary without the corners
    if (left){
        psc[0][0] += sqrt((par[0][0]+2*par[1][0])/par[2][0])*dt / (2 * dx * h0); // left front top
        psc[0][1] += sqrt(par[1][0]/par[2][0])*dt / (2 * dx * h0);
        psc[0][2] += sqrt(par[1][0]/par[2][0])*dt / (2 * dx * h0);
        psc[3][0] += sqrt((par[0][IXYZ(0,ny-1,0)]+2*par[1][IXYZ(0,ny-1,0)])/par[2][IXYZ(0,ny-1,0)])*dt / (2 * dx * h0); // left back top
        psc[3][1] += sqrt(par[1][IXYZ(0,ny-1,0)]/par[2][IXYZ(0,ny-1,0)])*dt / (2 * dx * h0);
        psc[3][2] += sqrt(par[1][IXYZ(0,ny-1,0)]/par[2][IXYZ(0,ny-1,0)])*dt / (2 * dx * h0);
        psc[7][0] += sqrt((par[0][IXYZ(0,ny-1,nz-1)]+2*par[1][IXYZ(0,ny-1,nz-1)])/par[2][IXYZ(0,ny-1,nz-1)])*dt / (2 * dx * h0); // left back bottom
        psc[7][1] += sqrt(par[1][IXYZ(0,ny-1,nz-1)]/par[2][IXYZ(0,ny-1,nz-1)])*dt / (2 * dx * h0);
        psc[7][2] += sqrt(par[1][IXYZ(0,ny-1,nz-1)]/par[2][IXYZ(0,ny-1,nz-1)])*dt / (2 * dx * h0);
        psc[4][0] += sqrt((par[0][IXYZ(0,0,nz-1)]+2*par[1][IXYZ(0,0,nz-1)])/par[2][IXYZ(0,0,nz-1)])*dt / (2 * dx * h0); // left front bottom
        psc[4][1] += sqrt(par[1][IXYZ(0,0,nz-1)]/par[2][IXYZ(0,0,nz-1)])*dt / (2 * dx * h0);
        psc[4][2] += sqrt(par[1][IXYZ(0,0,nz-1)]/par[2][IXYZ(0,0,nz-1)])*dt / (2 * dx * h0);

        if (ixmin==0){
            int ix=0;
            #pragma omp parallel for private(scalerx,scalery,scalerz)
            for (int iy = std::max(1,iymin); iy < std::min(ny-1,iymax); iy++){
                for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                    scalery = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    scalerz = scalery;
                    scalerx = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            // treat the edges without the corners
            if (!top){
                int iz=0;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iy = std::max(1,iymin); iy < std::min(ny-1,iymax); iy++){
                    scalerx = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    scalery = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx* h0);
                    scalerz = scalery;
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            if (!bottom){
                int iz=nz-1;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iy = std::max(1,iymin); iy < std::min(ny-1,iymax); iy++){
                    scalerx = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    scalery = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx* h0);
                    scalerz = scalery;
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            if (front){
                int iy=0;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                    scalerz = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0) + sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    scalery = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }
            else{
                int iy=0;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                    scalery = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    scalerz = scalery;
                    scalerx = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            if (back){
                int iy=ny-1;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                    scalerz = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0) + sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    scalery = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }
            else{
                int iy=ny-1;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                    scalery = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    scalerz = scalery;
                    scalerx = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }
        }
    }

    // scale the right boundary without the corners
    if (right){
        psc[1][0] += sqrt((par[0][IXYZ(nx-1,0,0)]+2*par[1][IXYZ(nx-1,0,0)])/par[2][IXYZ(nx-1,0,0)])*dt / (2 * dx * h0); // right front top
        psc[1][1] += sqrt(par[1][IXYZ(nx-1,0,0)]/par[2][IXYZ(nx-1,0,0)])*dt / (2 * dx * h0);
        psc[1][2] += sqrt(par[1][IXYZ(nx-1,0,0)]/par[2][IXYZ(nx-1,0,0)])*dt / (2 * dx * h0);
        psc[2][0] += sqrt((par[0][IXYZ(nx-1,ny-1,0)]+2*par[1][IXYZ(nx-1,ny-1,0)])/par[2][IXYZ(nx-1,ny-1,0)])*dt / (2 * dx * h0); // right back top
        psc[2][1] += sqrt(par[1][IXYZ(nx-1,ny-1,0)]/par[2][IXYZ(nx-1,ny-1,0)])*dt / (2 * dx * h0);
        psc[2][2] += sqrt(par[1][IXYZ(nx-1,ny-1,0)]/par[2][IXYZ(nx-1,ny-1,0)])*dt / (2 * dx * h0);
        psc[6][0] += sqrt((par[0][IXYZ(nx-1,ny-1,nz-1)]+2*par[1][IXYZ(nx-1,ny-1,nz-1)])/par[2][IXYZ(nx-1,ny-1,nz-1)])*dt / (2 * dx * h0); // right back bottom
        psc[6][1] += sqrt(par[1][IXYZ(nx-1,ny-1,nz-1)]/par[2][IXYZ(nx-1,ny-1,nz-1)])*dt / (2 * dx * h0);
        psc[6][2] += sqrt(par[1][IXYZ(nx-1,ny-1,nz-1)]/par[2][IXYZ(nx-1,ny-1,nz-1)])*dt / (2 * dx * h0);
        psc[5][0] += sqrt((par[0][IXYZ(nx-1,0,nz-1)]+2*par[1][IXYZ(nx-1,0,nz-1)])/par[2][IXYZ(nx-1,0,nz-1)])*dt / (2 * dx * h0); // right front bottom
        psc[5][1] += sqrt(par[1][IXYZ(nx-1,0,nz-1)]/par[2][IXYZ(nx-1,0,nz-1)])*dt / (2 * dx * h0);
        psc[5][2] += sqrt(par[1][IXYZ(nx-1,0,nz-1)]/par[2][IXYZ(nx-1,0,nz-1)])*dt / (2 * dx * h0);

        if (ixmax==nx){
            int ix=nx-1;
            #pragma omp parallel for private(scalerx,scalery,scalerz)
            for (int iy = std::max(1,iymin); iy < std::min(ny-1,iymax); iy++){
                for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                    scalery = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    scalerz = scalery;
                    scalerx = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            // treat the edges without the corners
            if (!top){
                int iz=0;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iy = std::max(1,iymin); iy < std::min(ny-1,iymax); iy++){
                    scalerx = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    scalery = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx* h0);
                    scalerz = scalery;
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            if (!bottom){
                int iz=nz-1;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iy = std::max(1,iymin); iy < std::min(ny-1,iymax); iy++){
                    scalerx = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    scalery = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx* h0);
                    scalerz = scalery;
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            if (front){
                int iy=0;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                    scalerz = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0) + sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    scalery = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }
            else{
                int iy=0;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                    scalery = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    scalerz = scalery;
                    scalerx = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            if (back){
                int iy=ny-1;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                    scalerz = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0) + sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    scalery = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0) + sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }
            else{
                int iy=ny-1;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                    scalery = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    scalerz = scalery;
                    scalerx = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dx * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }
        }
    }

    // scale the front boundary without the corners
    if (front){
        psc[0][0] += sqrt(par[1][0]/par[2][0])*dt / (2 * dy * h0); // front left top
        psc[0][1] += sqrt((par[0][0]+2*par[1][0])/par[2][0])*dt / (2 * dy * h0);
        psc[0][2] += sqrt(par[1][0]/par[2][0])*dt / (2 * dy * h0);
        psc[1][0] += sqrt(par[1][IXYZ(nx-1,0,0)]/par[2][IXYZ(nx-1,0,0)])*dt / (2 * dy * h0); // front right top
        psc[1][1] += sqrt((par[0][IXYZ(nx-1,0,0)]+2*par[1][IXYZ(nx-1,0,0)])/par[2][IXYZ(nx-1,0,0)])*dt / (2 * dy * h0);
        psc[1][2] += sqrt(par[1][IXYZ(nx-1,0,0)]/par[2][IXYZ(nx-1,0,0)])*dt / (2 * dy * h0);
        psc[5][0] += sqrt(par[1][IXYZ(nx-1,0,nz-1)]/par[2][IXYZ(nx-1,0,nz-1)])*dt / (2 * dy * h0); // front right bottom
        psc[5][1] += sqrt((par[0][IXYZ(nx-1,0,nz-1)]+2*par[1][IXYZ(nx-1,0,nz-1)])/par[2][IXYZ(nx-1,0,nz-1)])*dt / (2 * dy * h0);
        psc[5][2] += sqrt(par[1][IXYZ(nx-1,0,nz-1)]/par[2][IXYZ(nx-1,0,nz-1)])*dt / (2 * dy * h0);
        psc[4][0] += sqrt(par[1][IXYZ(0,0,nz-1)]/par[2][IXYZ(0,0,nz-1)])*dt / (2 * dy * h0); // front left bottom
        psc[4][1] += sqrt((par[0][IXYZ(0,0,nz-1)]+2*par[1][IXYZ(0,0,nz-1)])/par[2][IXYZ(0,0,nz-1)])*dt / (2 * dy * h0);
        psc[4][2] += sqrt(par[1][IXYZ(0,0,nz-1)]/par[2][IXYZ(0,0,nz-1)])*dt / (2 * dy * h0);

        if (iymin==0){
            int iy=0;
            #pragma omp parallel for private(scalerx,scalery,scalerz)
            for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    scalerz = scalerx;
                    scalery = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            // treat the edges without the corners
            if (!top){
                int iz=0;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                    scalery = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy* h0);
                    scalerz = scalerx;
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            if (!bottom){
                int iz=nz-1;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                    scalery = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy* h0);
                    scalerz = scalerx;
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            if (!left){
                int ix=0;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                    scalery = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy* h0);
                    scalerz = scalerx;
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            if (!right){
                int ix=nx-1;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                    scalery = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy* h0);
                    scalerz = scalerx;
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }
        }
    }

    // scale the back boundary without the corners
    if (back){
        psc[3][0] += sqrt(par[1][IXYZ(0,ny-1,0)]/par[2][IXYZ(0,ny-1,0)])*dt / (2 * dy * h0); // back left top
        psc[3][1] += sqrt((par[0][IXYZ(0,ny-1,0)]+2*par[1][IXYZ(0,ny-1,0)])/par[2][IXYZ(0,ny-1,0)])*dt / (2 * dy * h0);
        psc[3][2] += sqrt(par[1][IXYZ(0,ny-1,0)]/par[2][IXYZ(0,ny-1,0)])*dt / (2 * dy * h0);
        psc[2][0] += sqrt(par[1][IXYZ(nx-1,ny-1,0)]/par[2][IXYZ(nx-1,ny-1,0)])*dt / (2 * dy * h0); // back right top
        psc[2][1] += sqrt((par[0][IXYZ(nx-1,ny-1,0)]+2*par[1][IXYZ(nx-1,ny-1,0)])/par[2][IXYZ(nx-1,ny-1,0)])*dt / (2 * dy * h0);
        psc[2][2] += sqrt(par[1][IXYZ(nx-1,ny-1,0)]/par[2][IXYZ(nx-1,ny-1,0)])*dt / (2 * dy * h0);
        psc[6][0] += sqrt(par[1][IXYZ(nx-1,ny-1,nz-1)]/par[2][IXYZ(nx-1,ny-1,nz-1)])*dt / (2 * dy * h0); // back right bottom
        psc[6][1] += sqrt((par[0][IXYZ(nx-1,ny-1,nz-1)]+2*par[1][IXYZ(nx-1,ny-1,nz-1)])/par[2][IXYZ(nx-1,ny-1,nz-1)])*dt / (2 * dy * h0);
        psc[6][2] += sqrt(par[1][IXYZ(nx-1,ny-1,nz-1)]/par[2][IXYZ(nx-1,ny-1,nz-1)])*dt / (2 * dy * h0);
        psc[7][0] += sqrt(par[1][IXYZ(0,ny-1,nz-1)]/par[2][IXYZ(0,ny-1,nz-1)])*dt / (2 * dy * h0); // back left bottom
        psc[7][1] += sqrt((par[0][IXYZ(0,ny-1,nz-1)]+2*par[1][IXYZ(0,ny-1,nz-1)])/par[2][IXYZ(0,ny-1,nz-1)])*dt / (2 * dy * h0);
        psc[7][2] += sqrt(par[1][IXYZ(0,ny-1,nz-1)]/par[2][IXYZ(0,ny-1,nz-1)])*dt / (2 * dy * h0);

        if (iymax==ny){
            int iy=ny-1;
            #pragma omp parallel for private(scalerx,scalery,scalerz)
            for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    scalerz = scalerx;
                    scalery = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            // treat the edges without the corners
            if (!top){
                int iz=0;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                    scalery = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy* h0);
                    scalerz = scalerx;
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            if (!bottom){
                int iz=nz-1;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                    scalery = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy* h0);
                    scalerz = scalerx;
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            if (!left){
                int ix=0;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                    scalery = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy* h0);
                    scalerz = scalerx;
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }

            if (!right){
                int ix=nx-1;
                #pragma omp parallel for private(scalerx,scalery,scalerz)
                for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                    scalery = sqrt((par[0][IXYZ(ix,iy,iz)]+2*par[1][IXYZ(ix,iy,iz)])/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy * h0);
                    scalerx = sqrt(par[1][IXYZ(ix,iy,iz)]/par[2][IXYZ(ix,iy,iz)])*dt / (2 * dy* h0);
                    scalerz = scalerx;
                    in[0][IXYZ(ix,iy,iz)] /= (1+scalerx);
                    in[1][IXYZ(ix,iy,iz)] /= (1+scalery);
                    in[2][IXYZ(ix,iy,iz)] /= (1+scalerz);
                }
            }
        }
    }

    // scale the 8 corners
    // top front left
    in[0][IXYZ(0,0,0)] /= (1+psc[0][0]);
    in[1][IXYZ(0,0,0)] /= (1+psc[0][1]);
    in[2][IXYZ(0,0,0)] /= (1+psc[0][2]);
    // top front right
    in[0][IXYZ(nx-1,0,0)] /= (1+psc[1][0]);
    in[1][IXYZ(nx-1,0,0)] /= (1+psc[1][1]);
    in[2][IXYZ(nx-1,0,0)] /= (1+psc[1][2]);
    // top back right
    in[0][IXYZ(nx-1,ny-1,0)] /= (1+psc[2][0]);
    in[1][IXYZ(nx-1,ny-1,0)] /= (1+psc[2][1]);
    in[2][IXYZ(nx-1,ny-1,0)] /= (1+psc[2][2]);
    // top back left
    in[0][IXYZ(0,ny-1,0)] /= (1+psc[3][0]);
    in[1][IXYZ(0,ny-1,0)] /= (1+psc[3][1]);
    in[2][IXYZ(0,ny-1,0)] /= (1+psc[3][2]);
    // bottom front left
    in[0][IXYZ(0,0,nz-1)] /= (1+psc[4][0]);
    in[1][IXYZ(0,0,nz-1)] /= (1+psc[4][1]);
    in[2][IXYZ(0,0,nz-1)] /= (1+psc[4][2]);
    // bottom front right
    in[0][IXYZ(nx-1,0,nz-1)] /= (1+psc[5][0]);
    in[1][IXYZ(nx-1,0,nz-1)] /= (1+psc[5][1]);
    in[2][IXYZ(nx-1,0,nz-1)] /= (1+psc[5][2]);
    // bottom back right
    in[0][IXYZ(nx-1,ny-1,nz-1)] /= (1+psc[6][0]);
    in[1][IXYZ(nx-1,ny-1,nz-1)] /= (1+psc[6][1]);
    in[2][IXYZ(nx-1,ny-1,nz-1)] /= (1+psc[6][2]);
    // bottom back left
    in[0][IXYZ(0,ny-1,nz-1)] /= (1+psc[7][0]);
    in[1][IXYZ(0,ny-1,nz-1)] /= (1+psc[7][1]);
    in[2][IXYZ(0,ny-1,nz-1)] /= (1+psc[7][2]);
}

void taperz(data_t* in, int nx, int ny, int nz, int ncomp, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax, int icmin, int icmax, data_t a){

}
void taperx(data_t* in, int nx, int ny, int nz, int ncomp, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax, int icmin, int icmax, data_t a){
    
}
void tapery(data_t* in, int nx, int ny, int nz, int ncomp, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax, int icmin, int icmax, data_t a){
    
}

#undef IZ
#undef IX
#undef IY
#undef IXYZ