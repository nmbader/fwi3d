#pragma once

#include "vecReg.hpp"

#define IZ(iz) (iy*nx*nz+ix*nz+iz)
#define IX(ix) (iy*nx*nz+(ix)*nz+iz)
#define IY(iy) ((iy)*nx*nz+ix*nz+iz)
#define IXYZ(ix,iy,iz) ((iy)*nx*nz+(ix)*nz+iz)

typedef data_t (*expr)(const data_t ** par, int i);

// spatial quadrature operator H (or its inverse)
void applyHz(bool inv, bool add, const data_t * in, data_t * out, int nx, int ny, int nz, data_t d, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax);
void applyHx(bool inv, bool add, const data_t * in, data_t * out, int nx, int ny, int nz, data_t d, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax);
void applyHy(bool inv, bool add, const data_t * in, data_t * out, int nx, int ny, int nz, data_t d, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax);

// first derivative operators
void Dz(bool add, const data_t* in, data_t* out, int nx, int ny, int nz, data_t d, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax);
void Dx(bool add, const data_t* in, data_t* out, int nx, int ny, int nz, data_t d, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax);
void Dy(bool add, const data_t* in, data_t* out, int nx, int ny, int nz, data_t d, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax);

// first derivative operators pre-multiplied with variable parameter
void mult_Dz(bool add, const data_t* in, data_t* out, int nx, int ny, int nz, data_t d, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax, const data_t * par);
void mult_Dx(bool add, const data_t* in, data_t* out, int nx, int ny, int nz, data_t d, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax, const data_t * par);
void mult_Dy(bool add, const data_t* in, data_t* out, int nx, int ny, int nz, data_t d, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax, const data_t * par);

// Scale boundaries when absorbing SAT is used. This is needed for the time recursion of wavefields
void esat_scale_boundaries(data_t** in, int nx, int ny, int nz, data_t dx, data_t dy, data_t dz, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax, const data_t** par, data_t dt, bool top, bool bottom, bool left, bool right, bool front, bool back);

// Apply cosine^2 taper to damp the wavefield (to use in conjunction with absorbing SAT)
// taper = cos[a * pi/2 * (i - istart)/(iend-istart)]^2 ; 0 <= a <= 1
void taperz(data_t* __restrict in, int nx, int ny, int nz, int ncomp, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax, int icmin, int icmax, data_t a);
void taperx(data_t* __restrict in, int nx, int ny, int nz, int ncomp, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax, int icmin, int icmax, data_t a);
void tapery(data_t* __restrict in, int nx, int ny, int nz, int ncomp, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax, int icmin, int icmax, data_t a);

// second derivative operators with variable parameters, defined as template function to accomodate different expressions of parameters
template<expr f>
void Dzz_var(bool add, const data_t* in, data_t* __restrict out, int nx, int ny, int nz, data_t d, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax, const data_t ** par, data_t a){

    data_t coef[10] = {-3.0/4, -5.0/6, -1.0/24, 1.0/6, 1.0/2, 1.0/2, 1.0/6, -1.0/8, 1.0/6, -1.0/8};
    data_t bnd_coef[384] = {920.0/289,-59.0/68,-7549318.0/34157643,-17440994.0/184112825,0,0,0,0,-1740.0/289,0,295314153.0/366719282,262545878.0/1218454145,0,0,0,0,1128.0/289,59.0/68,-19250923.0/26254840,-12192537.0/324892213,0,0,0,0,-308.0/289,0,42283069.0/254173229,-43013531.0/427521546,0,0,0,0,0,0,-18700293.0/1355757959,18700293.0/1355757959,0,0,0,0,0,0,-3.0/833,3.0/833,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                12.0/17,0,89562243.0/385991318,16914691.0/272440171,0,0,0,0,-59.0/68,0,-47979680.0/48852831,-18034913.0/120051851,0,0,0,0,2.0/17,0,299262497.0/373256703,16156647.0/200473586,0,0,0,0,3.0/68,0,-14723969.0/177744748,22633571.0/584543661,0,0,0,0,0,0,46802031.0/1628311862,-46802031.0/1628311862,0,0,0,0,0,0,441.0/181507,-441.0/181507,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                -96.0/731,59.0/172,-47632751.0/164317197,-5570723.0/375470930,0,0,0,0,118.0/731,0,54282598.0/49343777,23802793.0/215253532,0,0,0,0,-16.0/731,-59.0/172,-39119273.0/25083370,-35971870.0/61324629,-26254.0/557679,0,0,0,-6.0/731,0,360454121.0/368940022,17254963.0/80047776,1500708.0/7993399,0,0,0,0,0,-18024731.0/79673021,24178273.0/88099647,-26254.0/185893,0,0,0,0,0,-870707.0/620833782,960119.0/1147305747,13564.0/23980197,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                -36.0/833,0,54741803.0/948483020,-13602043.0/389676498,0,0,0,0,177.0/3332,0,-35820026.0/359121865,24921773.0/534548210,0,0,0,0,-6.0/833,0,813284372.0/948584067,30057666.0/158897885,1500708.0/9108757,0,0,0,-9.0/3332,0,-95056924.0/90903639,-23417695.0/47008088,-7476412.0/9108757,-2.0/49,0,0,0,0,23159719.0/99948527,110687545.0/265515812,4502124.0/9108757,8.0/49,0,0,0,0,-3671038.0/2687426923,-1063649.0/8893843,1473580.0/9108757,-6.0/49,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,-9437957.0/1931986386,9437957.0/1931986386,0,0,0,0,0,0,17289851.0/489388053,-17289851.0/489388053,0,0,0,0,0,0,-66355919.0/327412264,24178273.0/98343792,-564461.0/4461432,0,0,0,0,0,17638343.0/74566894,103749401.0/243793650,375177.0/743572,1.0/6,0,0,0,0,-19321801.0/295845927,-50677283.0/62943042,-280535.0/371786,-5.0/6,-1.0/24,0,0,0,5130905.0/5183662662,35039615.0/213452232,1118749.0/2230716,1.0/2,1.0/6,0,0,0,0,0,-1.0/8,1.0/6,-1.0/8,0,0,0,0,0,0,0,0,0,
                0,0,-1.0/784,1.0/784,0,0,0,0,0,0,8673.0/2904112,-8673.0/2904112,0,0,0,0,0,0,-403062.0/320810033,960119.0/1280713392,3391.0/6692148,0,0,0,0,0,-1920494.0/1377228165,-1063649.0/8712336,368395.0/2230716,-1.0/8,0,0,0,0,5130905.0/5183662662,35039615.0/213452232,1118749.0/2230716,1.0/2,1.0/6,0,0,0,-117380.0/2351569839,-3290636.0/80044587,-5580181.0/6692148,-3.0/4,-5.0/6,-1.0/24,0,0,0,0,1.0/6,1.0/2,1.0/2,1.0/6,0,0,0,0,0,-1.0/8,1.0/6,-1.0/8
    };
    
    int nc1=6, nc2=8, nc3=8;
    int izminb=std::max(nc1,izmin);
    int izmaxb=std::min(nz - nc1, izmax);
    data_t  val=0;
    data_t ad2 = a/(d*d);
    
    #pragma omp parallel for private(val)
    for (int iy = iymin; iy < iymax; iy++){
        for (int ix = ixmin; ix < ixmax; ix++){

            int i1=iy*nx*nz+ix*nz;
            // apply the operator near the top boundary if included
            for (int iz=izmin; iz<nc1; iz++){
                out[IZ(iz)] = add*out[IZ(iz)];
                int i2=iz*nc2*nc3;
                for (int j = 0; j<nc2; j++){
                    val=0;
                    for (int k = 0; k<nc3; k++){
                        val += bnd_coef[i2+j*nc3+k] * f(par,i1+k);
                    }
                    out[IZ(iz)] += ad2 * val * in[IZ(j)];
                }
            }

            // apply the operator near the bottom boundary if included
            for (int iz=nz-izmax; iz<nc1; iz++){
                out[IZ(nz-1-iz)] = add*out[IZ(nz-1-iz)];
                int i2=iz*nc2*nc3;
                for (int j = 0; j<nc2; j++){
                    val=0;
                    for (int k = 0; k<nc3; k++){
                        val += bnd_coef[i2+j*nc3+k] * f(par,i1+nz-1-k);
                    }
                    out[IZ(nz-1-iz)] +=  ad2 * val * in[IZ(nz-1-j)];
                }
            }

            // apply the operator in the interior
            #pragma omp simd
            for (int iz = izminb; iz<izmaxb; iz++){
                out[IZ(iz)] = add*out[IZ(iz)]
                                            + ad2 * ( (coef[0]*f(par,i1+iz)+coef[1]*(f(par,i1+iz-1)+f(par,i1+iz+1))+coef[2]*(f(par,i1+iz-2)+f(par,i1+iz+2)))*in[IZ(iz)] 
                                            + (coef[3]*f(par,i1+iz-2)+coef[4]*f(par,i1+iz-1)+coef[5]*f(par,i1+iz)+coef[6]*f(par,i1+iz+1))*in[IZ(iz-1)]
                                            + (coef[3]*f(par,i1+iz+2)+coef[4]*f(par,i1+iz+1)+coef[5]*f(par,i1+iz)+coef[6]*f(par,i1+iz-1))*in[IZ(iz+1)]
                                            + (coef[7]*f(par,i1+iz-2)+coef[8]*f(par,i1+iz-1)+coef[9]*f(par,i1+iz))*in[IZ(iz-2)]
                                            + (coef[7]*f(par,i1+iz+2)+coef[8]*f(par,i1+iz+1)+coef[9]*f(par,i1+iz))*in[IZ(iz+2)]
                                            );
            }
        }
    }
}

template<expr f>
void Dxx_var(bool add, const data_t* in, data_t* __restrict out, int nx, int ny, int nz, data_t d, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax, const data_t ** par, data_t a){

    data_t coef[10] = {-3.0/4, -5.0/6, -1.0/24, 1.0/6, 1.0/2, 1.0/2, 1.0/6, -1.0/8, 1.0/6, -1.0/8};
    data_t bnd_coef[384] = {920.0/289,-59.0/68,-7549318.0/34157643,-17440994.0/184112825,0,0,0,0,-1740.0/289,0,295314153.0/366719282,262545878.0/1218454145,0,0,0,0,1128.0/289,59.0/68,-19250923.0/26254840,-12192537.0/324892213,0,0,0,0,-308.0/289,0,42283069.0/254173229,-43013531.0/427521546,0,0,0,0,0,0,-18700293.0/1355757959,18700293.0/1355757959,0,0,0,0,0,0,-3.0/833,3.0/833,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                12.0/17,0,89562243.0/385991318,16914691.0/272440171,0,0,0,0,-59.0/68,0,-47979680.0/48852831,-18034913.0/120051851,0,0,0,0,2.0/17,0,299262497.0/373256703,16156647.0/200473586,0,0,0,0,3.0/68,0,-14723969.0/177744748,22633571.0/584543661,0,0,0,0,0,0,46802031.0/1628311862,-46802031.0/1628311862,0,0,0,0,0,0,441.0/181507,-441.0/181507,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                -96.0/731,59.0/172,-47632751.0/164317197,-5570723.0/375470930,0,0,0,0,118.0/731,0,54282598.0/49343777,23802793.0/215253532,0,0,0,0,-16.0/731,-59.0/172,-39119273.0/25083370,-35971870.0/61324629,-26254.0/557679,0,0,0,-6.0/731,0,360454121.0/368940022,17254963.0/80047776,1500708.0/7993399,0,0,0,0,0,-18024731.0/79673021,24178273.0/88099647,-26254.0/185893,0,0,0,0,0,-870707.0/620833782,960119.0/1147305747,13564.0/23980197,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                -36.0/833,0,54741803.0/948483020,-13602043.0/389676498,0,0,0,0,177.0/3332,0,-35820026.0/359121865,24921773.0/534548210,0,0,0,0,-6.0/833,0,813284372.0/948584067,30057666.0/158897885,1500708.0/9108757,0,0,0,-9.0/3332,0,-95056924.0/90903639,-23417695.0/47008088,-7476412.0/9108757,-2.0/49,0,0,0,0,23159719.0/99948527,110687545.0/265515812,4502124.0/9108757,8.0/49,0,0,0,0,-3671038.0/2687426923,-1063649.0/8893843,1473580.0/9108757,-6.0/49,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,-9437957.0/1931986386,9437957.0/1931986386,0,0,0,0,0,0,17289851.0/489388053,-17289851.0/489388053,0,0,0,0,0,0,-66355919.0/327412264,24178273.0/98343792,-564461.0/4461432,0,0,0,0,0,17638343.0/74566894,103749401.0/243793650,375177.0/743572,1.0/6,0,0,0,0,-19321801.0/295845927,-50677283.0/62943042,-280535.0/371786,-5.0/6,-1.0/24,0,0,0,5130905.0/5183662662,35039615.0/213452232,1118749.0/2230716,1.0/2,1.0/6,0,0,0,0,0,-1.0/8,1.0/6,-1.0/8,0,0,0,0,0,0,0,0,0,
                0,0,-1.0/784,1.0/784,0,0,0,0,0,0,8673.0/2904112,-8673.0/2904112,0,0,0,0,0,0,-403062.0/320810033,960119.0/1280713392,3391.0/6692148,0,0,0,0,0,-1920494.0/1377228165,-1063649.0/8712336,368395.0/2230716,-1.0/8,0,0,0,0,5130905.0/5183662662,35039615.0/213452232,1118749.0/2230716,1.0/2,1.0/6,0,0,0,-117380.0/2351569839,-3290636.0/80044587,-5580181.0/6692148,-3.0/4,-5.0/6,-1.0/24,0,0,0,0,1.0/6,1.0/2,1.0/2,1.0/6,0,0,0,0,0,-1.0/8,1.0/6,-1.0/8
    };

    int nc1=6, nc2=8, nc3=8;
    int ixminb=std::max(nc1,ixmin);
    int ixmaxb=std::min(nx - nc1, ixmax);
    data_t  val=0;
    data_t ad2 = a/(d*d);
    
    #pragma omp parallel for private( val)
    for (int iy = iymin; iy < iymax; iy++){
        // apply the operator near the left boundary if included
        for (int ix = ixmin; ix < nc1; ix++){
            int i2=ix*nc2*nc3;
            for (int iz=izmin; iz<izmax; iz++){
                out[IX(ix)] = add*out[IX(ix)];
                for (int j = 0; j<nc2; j++){
                    val=0;
                    for (int k = 0; k<nc3; k++){
                        val += bnd_coef[i2+j*nc3+k] * f(par,IX(k));
                    }
                    out[IX(ix)] += ad2 * val * in[IX(j)];
                }
            }
        }

        // apply the operator near the right boundary if included
        for (int ix = nx-ixmax; ix < nc1; ix++){
            int i2=ix*nc2*nc3;
            for (int iz=izmin; iz<izmax; iz++){
                out[IX(nx-1-ix)] = add*out[IX(nx-1-ix)];
                for (int j = 0; j<nc2; j++){
                    val=0;
                    for (int k = 0; k<nc3; k++){
                        val += bnd_coef[i2+j*nc3+k] * f(par,IX(nx-1-k));
                    }
                    out[IX(nx-1-ix)] += ad2 * val * in[IX(nx-1-j)];
                }
            }
        }

        // apply the operator in the interior
        for (int ix=ixminb; ix<ixmaxb; ix++){
            #pragma omp simd
            for (int iz = izmin; iz<izmax; iz++){
                out[ IX(ix)] = add*out[ IX(ix)] + ad2 * ( (coef[0]* f(par, IX(ix))+coef[1]*( f(par,IX(ix-1))+ f(par,IX(ix+1)))+coef[2]*( f(par,IX(ix-2))+ f(par,IX(ix+2))))*in[ IX(ix)] 
                                            + (coef[3]* f(par,IX(ix-2))+coef[4]* f(par,IX(ix-1))+coef[5]* f(par, IX(ix))+coef[6]* f(par,IX(ix+1)))*in[IX(ix-1)]
                                            + (coef[3]* f(par,IX(ix+2))+coef[4]* f(par,IX(ix+1))+coef[5]* f(par, IX(ix))+coef[6]* f(par,IX(ix-1)))*in[IX(ix+1)]
                                            + (coef[7]* f(par,IX(ix-2))+coef[8]* f(par,IX(ix-1))+coef[9]* f(par, IX(ix)))*in[IX(ix-2)]
                                            + (coef[7]* f(par,IX(ix+2))+coef[8]* f(par,IX(ix+1))+coef[9]* f(par, IX(ix)))*in[IX(ix+2)]
                                            );
            }
        }
    }
}

template<expr f>
void Dyy_var(bool add, const data_t* in, data_t* __restrict out, int nx, int ny, int nz, data_t d, int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax, const data_t ** par, data_t a){

    data_t coef[10] = {-3.0/4, -5.0/6, -1.0/24, 1.0/6, 1.0/2, 1.0/2, 1.0/6, -1.0/8, 1.0/6, -1.0/8};
    data_t bnd_coef[384] = {920.0/289,-59.0/68,-7549318.0/34157643,-17440994.0/184112825,0,0,0,0,-1740.0/289,0,295314153.0/366719282,262545878.0/1218454145,0,0,0,0,1128.0/289,59.0/68,-19250923.0/26254840,-12192537.0/324892213,0,0,0,0,-308.0/289,0,42283069.0/254173229,-43013531.0/427521546,0,0,0,0,0,0,-18700293.0/1355757959,18700293.0/1355757959,0,0,0,0,0,0,-3.0/833,3.0/833,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                12.0/17,0,89562243.0/385991318,16914691.0/272440171,0,0,0,0,-59.0/68,0,-47979680.0/48852831,-18034913.0/120051851,0,0,0,0,2.0/17,0,299262497.0/373256703,16156647.0/200473586,0,0,0,0,3.0/68,0,-14723969.0/177744748,22633571.0/584543661,0,0,0,0,0,0,46802031.0/1628311862,-46802031.0/1628311862,0,0,0,0,0,0,441.0/181507,-441.0/181507,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                -96.0/731,59.0/172,-47632751.0/164317197,-5570723.0/375470930,0,0,0,0,118.0/731,0,54282598.0/49343777,23802793.0/215253532,0,0,0,0,-16.0/731,-59.0/172,-39119273.0/25083370,-35971870.0/61324629,-26254.0/557679,0,0,0,-6.0/731,0,360454121.0/368940022,17254963.0/80047776,1500708.0/7993399,0,0,0,0,0,-18024731.0/79673021,24178273.0/88099647,-26254.0/185893,0,0,0,0,0,-870707.0/620833782,960119.0/1147305747,13564.0/23980197,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                -36.0/833,0,54741803.0/948483020,-13602043.0/389676498,0,0,0,0,177.0/3332,0,-35820026.0/359121865,24921773.0/534548210,0,0,0,0,-6.0/833,0,813284372.0/948584067,30057666.0/158897885,1500708.0/9108757,0,0,0,-9.0/3332,0,-95056924.0/90903639,-23417695.0/47008088,-7476412.0/9108757,-2.0/49,0,0,0,0,23159719.0/99948527,110687545.0/265515812,4502124.0/9108757,8.0/49,0,0,0,0,-3671038.0/2687426923,-1063649.0/8893843,1473580.0/9108757,-6.0/49,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,-9437957.0/1931986386,9437957.0/1931986386,0,0,0,0,0,0,17289851.0/489388053,-17289851.0/489388053,0,0,0,0,0,0,-66355919.0/327412264,24178273.0/98343792,-564461.0/4461432,0,0,0,0,0,17638343.0/74566894,103749401.0/243793650,375177.0/743572,1.0/6,0,0,0,0,-19321801.0/295845927,-50677283.0/62943042,-280535.0/371786,-5.0/6,-1.0/24,0,0,0,5130905.0/5183662662,35039615.0/213452232,1118749.0/2230716,1.0/2,1.0/6,0,0,0,0,0,-1.0/8,1.0/6,-1.0/8,0,0,0,0,0,0,0,0,0,
                0,0,-1.0/784,1.0/784,0,0,0,0,0,0,8673.0/2904112,-8673.0/2904112,0,0,0,0,0,0,-403062.0/320810033,960119.0/1280713392,3391.0/6692148,0,0,0,0,0,-1920494.0/1377228165,-1063649.0/8712336,368395.0/2230716,-1.0/8,0,0,0,0,5130905.0/5183662662,35039615.0/213452232,1118749.0/2230716,1.0/2,1.0/6,0,0,0,-117380.0/2351569839,-3290636.0/80044587,-5580181.0/6692148,-3.0/4,-5.0/6,-1.0/24,0,0,0,0,1.0/6,1.0/2,1.0/2,1.0/6,0,0,0,0,0,-1.0/8,1.0/6,-1.0/8
    };

    int nc1=6, nc2=8, nc3=8;
    int iyminb=std::max(nc1,iymin);
    int iymaxb=std::min(ny - nc1, iymax);
    data_t  val=0;
    data_t ad2 = a/(d*d);
    
    // apply the operator near the front boundary if included
    #pragma omp parallel for private( val)
    for (int iy = iymin; iy < nc1; iy++){
        int i2=iy*nc2*nc3;
        for (int ix = ixmin; ix < ixmax; ix++){
            for (int iz=izmin; iz<izmax; iz++){
                out[IY(iy)] = add*out[IY(iy)];
                for (int j = 0; j<nc2; j++){
                    val=0;
                    for (int k = 0; k<nc3; k++){
                        val += bnd_coef[i2+j*nc3+k] * f(par,IY(k));
                    }
                    out[IY(iy)] += ad2 * val * in[IY(j)];
                }
            }
        }
    }

    // apply the operator near the back boundary if included
    #pragma omp parallel for private( val)
    for (int iy = ny-iymax; iy < nc1; iy++){
        int i2=iy*nc2*nc3;
        for (int ix = ixmin; ix < ixmax; ix++){
            for (int iz=izmin; iz<izmax; iz++){
                out[IY(ny-1-iy)] = add*out[IY(ny-1-iy)];
                for (int j = 0; j<nc2; j++){
                    val=0;
                    for (int k = 0; k<nc3; k++){
                        val += bnd_coef[i2+j*nc3+k] * f(par,IY(ny-1-k));
                    }
                    out[IY(ny-1-iy)] += ad2 * val * in[IY(ny-1-j)];
                }
            }
        }
    }

    // apply the operator in the interior
    #pragma omp parallel for private( val)
    for (int iy = iyminb; iy < iymaxb; iy++){
        for (int ix=ixmin; ix<ixmax; ix++){
            #pragma omp simd
            for (int iz = izmin; iz<izmax; iz++){
                out[ IY(iy)] = add*out[ IY(iy)] + ad2 * ( (coef[0]* f(par, IY(iy))+coef[1]*( f(par,IY(iy-1))+ f(par,IY(iy+1)))+coef[2]*( f(par,IY(iy-2))+ f(par,IY(iy+2))))*in[ IY(iy)] 
                                            + (coef[3]* f(par,IY(iy-2))+coef[4]* f(par,IY(iy-1))+coef[5]* f(par, IY(iy))+coef[6]* f(par,IY(iy+1)))*in[IY(iy-1)]
                                            + (coef[3]* f(par,IY(iy+2))+coef[4]* f(par,IY(iy+1))+coef[5]* f(par, IY(iy))+coef[6]* f(par,IY(iy-1)))*in[IY(iy+1)]
                                            + (coef[7]* f(par,IY(iy-2))+coef[8]* f(par,IY(iy-1))+coef[9]* f(par, IY(iy)))*in[IY(iy-2)]
                                            + (coef[7]* f(par,IY(iy+2))+coef[8]* f(par,IY(iy+1))+coef[9]* f(par, IY(iy)))*in[IY(iy+2)]
                                            );
            }
        }
    }
}

// Neumann SAT to impose free surface BC for elastic WE
// defined as template function to accomodate different expressions of parameters
template<expr f1, expr f2, expr f3, expr f4, expr f5>
void esat_neumann_top(bool add, const data_t** in, __restrict data_t* out, int nx, int ny, int nz, data_t dx, data_t dy, data_t dz, data_t dt, int ixmin, int ixmax, int iymin, int iymax, const data_t ** par){

    data_t coef[2] = {2.0/3,-1.0/12};
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3}; 
    data_t h0 = 17.0/48;
    int nc1=4, nc2=6;
    int ixminb=std::max(nc1,ixmin);
    int ixmaxb=std::min(nx - nc1, ixmax);
    int iyminb=std::max(nc1,iymin);
    int iymaxb=std::min(ny - nc1, iymax);
    data_t adh = 1.0/(dz*h0);

    data_t sumx, sumy, sumz, sumz0;

    // SAT = - Hz-1.(-f1.Dx.in1 -f2.Dy.in2 -f3.Dz.in3 + f4.Sz.in3 - f5.in4/dt)_0
    // Sz is the boundary derivative operator pointing outwards

    // top left front
    #pragma omp parallel for private(sumx,sumy,sumz,sumz0)
    for (int iy=iymin; iy<nc1; iy++){
        for (int ix=ixmin; ix<nc1; ix++){

            // (Dx.in1)_0 and (Dy.in2)_0
            int iz=0;
            sumx=0;
            sumy=0;
            for (int j=0; j<nc2; j++){
                sumx += bnd_coef[ix*nc2+j] * in[0][IX(j)];
                sumy += bnd_coef[iy*nc2+j] * in[1][IY(j)];
            }
            // (Sz.in3)_0
            sumz = 0;
            for (iz = 0; iz < 4; iz++){
                sumz += scoef[iz] * in[2][IZ(iz)];
            }
            // (Dz.in3)_0
            sumz0=0;
            for (iz = 0; iz < 6; iz++){
                sumz0 += bnd_coef[iz] * in[2][IZ(iz)];
            }
            iz=0;
            out[IZ(iz)] = add*out[IZ(iz)] - adh * (-f1(par,IZ(iz))*sumx/dx -f2(par,IZ(iz))*sumy/dy -f3(par,IZ(iz))*sumz0/dz + f4(par,IZ(iz))*sumz/dz - f5(par,IZ(iz))*in[3][IZ(iz)]/dt);
        }
    }

    // top left back
    #pragma omp parallel for private(sumx,sumy,sumz,sumz0)
    for (int iy=ny-iymax; iy<std::min(nc1,ny-iymin); iy++){
        for (int ix=ixmin; ix<nc1; ix++){

            // (Dx.in1)_0 and (Dy.in2)_0
            int iz=0;
            sumx=0;
            sumy=0;
            for (int j=0; j<nc2; j++){
                sumx += bnd_coef[ix*nc2+j] * in[0][IXYZ(j,ny-1-iy,iz)];
                sumy -= bnd_coef[iy*nc2+j] * in[1][IY(ny-1-j)];
            }

            // (Sz.in3)_0
            sumz = 0;
            for (iz = 0; iz < 4; iz++){
                sumz += scoef[iz] * in[2][IXYZ(ix,ny-1-iy,iz)];
            }
            // (Dz.in3)_0
            sumz0 = 0;
            for (iz = 0; iz < 6; iz++){
                sumz0 += bnd_coef[iz] * in[2][IXYZ(ix,ny-1-iy,iz)];
            }
            iz=0;
            out[IXYZ(ix,ny-1-iy,iz)] = add*out[IXYZ(ix,ny-1-iy,iz)] - adh * (-f1(par,IXYZ(ix,ny-1-iy,iz))*sumx/dx -f2(par,IXYZ(ix,ny-1-iy,iz))*sumy/dy -f3(par,IXYZ(ix,ny-1-iy,iz))*sumz0/dz + f4(par,IXYZ(ix,ny-1-iy,iz))*sumz/dz -f5(par,IXYZ(ix,ny-1-iy,iz))*in[3][IXYZ(ix,ny-1-iy,iz)]/dt);
        }
    }

    // top left middle
    #pragma omp parallel for private(sumx,sumy,sumz,sumz0)
    for (int iy=iyminb; iy<iymaxb; iy++){
        for (int ix=ixmin; ix<nc1; ix++){

            // (Dx.in1)_0 and (Dy.in2)_0
            int iz=0;
            sumx=0;
            for (int j=0; j<nc2; j++){
                sumx += bnd_coef[ix*nc2+j] * in[0][IX(j)];
            }
            sumy = coef[0] * (in[1][IY(iy+1)]-in[1][IY(iy-1)]) + coef[1] * (in[1][IY(iy+2)]-in[1][IY(iy-2)]);

            // (Sz.in3)_0
            sumz = 0;
            for (iz = 0; iz < 4; iz++){
                sumz += scoef[iz] * in[2][IZ(iz)];
            }
            // (Dz.in3)_0
            sumz0 = 0;
            for (iz = 0; iz < 6; iz++){
                sumz0 += bnd_coef[iz] * in[2][IZ(iz)];
            }
            iz=0;
            out[IZ(iz)] = add*out[IZ(iz)] - adh * (-f1(par,IZ(iz))*sumx/dx -f2(par,IZ(iz))*sumy/dy -f3(par,IZ(iz))*sumz0/dz + f4(par,IZ(iz))*sumz/dz - f5(par,IZ(iz))*in[3][IZ(iz)]/dt);
        }
    }

    // top right front
    #pragma omp parallel for private(sumx,sumy,sumz,sumz0)
    for (int iy=iymin; iy<nc1; iy++){
        for (int ix=nx-ixmax; ix<std::min(nc1,nx-ixmin); ix++){

            // (Dx.in1)_0 and (Dy.in2)_0
            int iz=0;
            sumx=0;
            sumy=0;
            for (int j=0; j<nc2; j++){
                sumx -= bnd_coef[ix*nc2+j] * in[0][IX(nx-1-j)];
                sumy += bnd_coef[iy*nc2+j] * in[1][IXYZ(nx-1-ix,j,iz)];
            }

            // (Sz.in3)_0
            sumz = 0;
            for (iz = 0; iz < 4; iz++){
                sumz += scoef[iz] * in[2][IXYZ(nx-1-ix,iy,iz)];
            }
            // (Dz.in3)_0
            sumz0 = 0;
            for (iz = 0; iz < 6; iz++){
                sumz0 += bnd_coef[iz] * in[2][IXYZ(nx-1-ix,iy,iz)];
            }
            iz=0;
            out[IXYZ(nx-1-ix,iy,iz)] = add*out[IXYZ(nx-1-ix,iy,iz)] - adh * (-f1(par,IXYZ(nx-1-ix,iy,iz))*sumx/dx -f2(par,IXYZ(nx-1-ix,iy,iz))*sumy/dy -f3(par,IXYZ(nx-1-ix,iy,iz))*sumz0/dz + f4(par,IXYZ(nx-1-ix,iy,iz))*sumz/dz - f5(par,IXYZ(nx-1-ix,iy,iz))*in[3][IXYZ(nx-1-ix,iy,iz)]/dt);
        }
    }

    // top right back
    #pragma omp parallel for private(sumx,sumy,sumz,sumz0)
    for (int iy=ny-iymax; iy<std::min(nc1,ny-iymin); iy++){
        for (int ix=nx-ixmax; ix<std::min(nc1,nx-ixmin); ix++){

            // (Dx.in1)_0 and (Dy.in2)_0
            int iz=0;
            sumx=0;
            sumy=0;
            for (int j=0; j<nc2; j++){
                sumx -= bnd_coef[ix*nc2+j] * in[0][IXYZ(nx-1-j,ny-1-iy,iz)];
                sumy -= bnd_coef[iy*nc2+j] * in[1][IXYZ(nx-1-ix,ny-1-j,iz)];
            }

            // (Sz.in3)_0
            sumz = 0;
            for (iz = 0; iz < 4; iz++){
                sumz += scoef[iz] * in[2][IXYZ(nx-1-ix,ny-1-iy,iz)];
            }
            // (Dz.in3)_0
            sumz0 = 0;
            for (iz = 0; iz < 6; iz++){
                sumz0 += bnd_coef[iz] * in[2][IXYZ(nx-1-ix,ny-1-iy,iz)];
            }
            iz=0;
            out[IXYZ(nx-1-ix,ny-1-iy,iz)] = add*out[IXYZ(nx-1-ix,ny-1-iy,iz)] - adh * (-f1(par,IXYZ(nx-1-ix,ny-1-iy,iz))*sumx/dx -f2(par,IXYZ(nx-1-ix,ny-1-iy,iz))*sumy/dy -f3(par,IXYZ(nx-1-ix,ny-1-iy,iz))*sumz0/dz + f4(par,IXYZ(nx-1-ix,ny-1-iy,iz))*sumz/dz - f5(par,IXYZ(nx-1-ix,ny-1-iy,iz))*in[3][IXYZ(nx-1-ix,ny-1-iy,iz)]/dt);
        }
    }

    // top right middle
    #pragma omp parallel for private(sumx,sumy,sumz,sumz0)
    for (int iy=iyminb; iy<iymaxb; iy++){
        for (int ix=nx-ixmax; ix<std::min(nc1,nx-ixmin); ix++){

            // (Dx.in1)_0 and (Dy.in2)_0
            int iz=0;
            sumx=0;
            for (int j=0; j<nc2; j++){
                sumx -= bnd_coef[ix*nc2+j] * in[0][IX(nx-1-j)];
            }
            sumy = coef[0] * (in[1][IXYZ(nx-1-ix,iy+1,iz)]-in[1][IXYZ(nx-1-ix,iy-1,iz)]) + coef[1] * (in[1][IXYZ(nx-1-ix,iy+2,iz)]-in[1][IXYZ(nx-1-ix,iy-2,iz)]);

            // (Sz.in3)_0
            sumz = 0;
            for (iz = 0; iz < 4; iz++){
                sumz += scoef[iz] * in[2][IXYZ(nx-1-ix,iy,iz)];
            }
            // (Dz.in3)_0
            sumz0 = 0;
            for (iz = 0; iz < 6; iz++){
                sumz0 += bnd_coef[iz] * in[2][IXYZ(nx-1-ix,iy,iz)];
            }
            iz=0;
            out[IXYZ(nx-1-ix,iy,iz)] = add*out[IXYZ(nx-1-ix,iy,iz)] - adh * (-f1(par,IXYZ(nx-1-ix,iy,iz))*sumx/dx -f2(par,IXYZ(nx-1-ix,iy,iz))*sumy/dy -f3(par,IXYZ(nx-1-ix,iy,iz))*sumz0/dz + f4(par,IXYZ(nx-1-ix,iy,iz))*sumz/dz - f5(par,IXYZ(nx-1-ix,iy,iz))*in[3][IXYZ(nx-1-ix,iy,iz)]/dt);
        }
    }

    // top middle front
    #pragma omp parallel for private(sumx,sumy,sumz,sumz0)
    for (int iy=iymin; iy<nc1; iy++){
        for (int ix=ixminb; ix<ixmaxb; ix++){

            // (Dx.in1)_0 and (Dy.in2)_0
            int iz=0;
            sumx = coef[0] * (in[0][IX(ix+1)]-in[0][IX(ix-1)]) + coef[1] * (in[0][IX(ix+2)]-in[0][IX(ix-2)]);
            sumy=0;
            for (int j=0; j<nc2; j++){
                sumy += bnd_coef[iy*nc2+j] * in[1][IY(j)];
            }

            // (Sz.in3)_0
            sumz = 0;
            for (iz = 0; iz < 4; iz++){
                sumz += scoef[iz] * in[2][IZ(iz)];
            }
            // (Dz.in3)_0
            sumz0 = 0;
            for (iz = 0; iz < 6; iz++){
                sumz0 += bnd_coef[iz] * in[2][IZ(iz)];
            }
            iz=0;
            out[IZ(iz)] = add*out[IZ(iz)] - adh * (-f1(par,IZ(iz))*sumx/dx -f2(par,IZ(iz))*sumy/dy -f3(par,IZ(iz))*sumz0/dz + f4(par,IZ(iz))*sumz/dz - f5(par,IZ(iz))*in[3][IZ(iz)]/dt);
        }
    }

    // top middle back
    #pragma omp parallel for private(sumx,sumy,sumz,sumz0)
    for (int iy=ny-iymax; iy<std::min(nc1,ny-iymin); iy++){
        for (int ix=ixminb; ix<ixmaxb; ix++){

            // (Dx.in1)_0 and (Dy.in2)_0
            int iz=0;
            sumx = coef[0] * (in[0][IXYZ(ix+1,ny-1-iy,iz)]-in[0][IXYZ(ix-1,ny-1-iy,iz)]) + coef[1] * (in[0][IXYZ(ix+2,ny-1-iy,iz)]-in[0][IXYZ(ix-2,ny-1-iy,iz)]);
            sumy=0;
            for (int j=0; j<nc2; j++){
                sumy -= bnd_coef[iy*nc2+j] * in[1][IY(ny-1-j)];
            }

            // (Sz.in3)_0
            sumz = 0;
            for (iz = 0; iz < 4; iz++){
                sumz += scoef[iz] * in[2][IXYZ(ix,ny-1-iy,iz)];
            }
            // (Dz.in3)_0
            sumz0 = 0;
            for (iz = 0; iz < 6; iz++){
                sumz0 += bnd_coef[iz] * in[2][IXYZ(ix,ny-1-iy,iz)];
            }
            iz=0;
            out[IXYZ(ix,ny-1-iy,iz)] = add*out[IXYZ(ix,ny-1-iy,iz)] - adh * (-f1(par,IXYZ(ix,ny-1-iy,iz))*sumx/dx -f2(par,IXYZ(ix,ny-1-iy,iz))*sumy/dy -f3(par,IXYZ(ix,ny-1-iy,iz))*sumz0/dz + f4(par,IXYZ(ix,ny-1-iy,iz))*sumz/dz - f5(par,IXYZ(ix,ny-1-iy,iz))*in[3][IXYZ(ix,ny-1-iy,iz)]/dt);
        }
    }

    // top middle
    #pragma omp parallel for private(sumx, sumy, sumz, sumz0)
    for (int iy=iyminb; iy<iymaxb; iy++){
        for (int ix=ixminb; ix<ixmaxb; ix++){

            // (Dx.in1)_0 and (Dy.in2)_0
            int iz=0;
            sumx = coef[0] * (in[0][IX(ix+1)]-in[0][IX(ix-1)]) + coef[1] * (in[0][IX(ix+2)]-in[0][IX(ix-2)]);
            sumy = coef[0] * (in[1][IY(iy+1)]-in[1][IY(iy-1)]) + coef[1] * (in[1][IY(iy+2)]-in[1][IY(iy-2)]);

            // (Sz.in3)_0
            sumz = 0;
            for (iz = 0; iz < 4; iz++){
                sumz += scoef[iz] * in[2][IZ(iz)];
            }
            // (Dz.in3)_0
            sumz0 = 0;
            for (iz = 0; iz < 6; iz++){
                sumz0 += bnd_coef[iz] * in[2][IZ(iz)];
            }
            iz=0;
            out[IZ(iz)] = add*out[IZ(iz)] - adh * (-f1(par,IZ(iz))*sumx/dx -f2(par,IZ(iz))*sumy/dy -f3(par,IZ(iz))*sumz0/dz + f4(par,IZ(iz))*sumz/dz - f5(par,IZ(iz))*in[3][IZ(iz)]/dt);
        }
    }
}

template<expr f1, expr f2, expr f3, expr f4, expr f5>
void esat_neumann_bottom(bool add, const data_t** in, __restrict data_t* out, int nx, int ny, int nz, data_t dx, data_t dy, data_t dz, data_t dt, int ixmin, int ixmax, int iymin, int iymax, const data_t ** par){

    data_t coef[2] = {2.0/3,-1.0/12};
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3}; 
    data_t h0 = 17.0/48;
    int nc1=4, nc2=6;
    int ixminb=std::max(nc1,ixmin);
    int ixmaxb=std::min(nx - nc1, ixmax);
    int iyminb=std::max(nc1,iymin);
    int iymaxb=std::min(ny - nc1, iymax);
    data_t adh = 1.0/(dz*h0);

    data_t sumx, sumy, sumz, sumz0;

    // SAT = - Hz-1.(f1.Dx.in1 + f2.Dy.in2 + f3.Dz.in3 + f4.Sz.in3 - f5.in4/dt)_0
    // Sz is the boundary derivative operator pointing outwards

    // bottom left front
    #pragma omp parallel for private(sumx,sumy,sumz,sumz0)
    for (int iy=iymin; iy<nc1; iy++){
        for (int ix=ixmin; ix<nc1; ix++){

            // (Dx.in1)_0 and (Dy.in2)_0
            int iz=nz-1;
            sumx=0;
            sumy=0;
            for (int j=0; j<nc2; j++){
                sumx += bnd_coef[ix*nc2+j] * in[0][IX(j)];
                sumy += bnd_coef[iy*nc2+j] * in[1][IY(j)];
            }

            // (Sz.in3)_0
            sumz = 0;
            for (iz = 0; iz < 4; iz++){
                sumz += scoef[iz] * in[2][IZ(nz-1-iz)];
            }
            // (Dz.in3)_0
            sumz0 = 0;
            for (iz = 0; iz < 6; iz++){
                sumz0 -= bnd_coef[iz] * in[2][IZ(nz-1-iz)];
            }
            iz=nz-1;
            out[IZ(iz)] = add*out[IZ(iz)] - adh * (f1(par,IZ(iz))*sumx/dx + f2(par,IZ(iz))*sumy/dy + f3(par,IZ(iz))*sumz0/dz + f4(par,IZ(iz))*sumz/dz -f5(par,IZ(iz))*in[3][IZ(iz)]/dt);
        }
    }

    // bottom left back
    #pragma omp parallel for private(sumx,sumy,sumz,sumz0)
    for (int iy=ny-iymax; iy<std::min(nc1,ny-iymin); iy++){
        for (int ix=ixmin; ix<nc1; ix++){

            // (Dx.in1)_0 and (Dy.in2)_0
            int iz=nz-1;
            sumx=0;
            sumy=0;
            for (int j=0; j<nc2; j++){
                sumx += bnd_coef[ix*nc2+j] * in[0][IXYZ(j,ny-1-iy,iz)];
                sumy -= bnd_coef[iy*nc2+j] * in[1][IY(ny-1-j)];
            }

            // (Sz.in3)_0
            sumz = 0;
            for (iz = 0; iz < 4; iz++){
                sumz += scoef[iz] * in[2][IXYZ(ix,ny-1-iy,nz-1-iz)];
            }
            // (Dz.in3)_0
            sumz0 = 0;
            for (iz = 0; iz < 6; iz++){
                sumz0 -= bnd_coef[iz] * in[2][IXYZ(ix,ny-1-iy,nz-1-iz)];
            }
            iz=nz-1;
            out[IXYZ(ix,ny-1-iy,iz)] = add*out[IXYZ(ix,ny-1-iy,iz)] - adh * (f1(par,IXYZ(ix,ny-1-iy,iz))*sumx/dx +f2(par,IXYZ(ix,ny-1-iy,iz))*sumy/dy +f3(par,IXYZ(ix,ny-1-iy,iz))*sumz0/dz + f4(par,IXYZ(ix,ny-1-iy,iz))*sumz/dz -f5(par,IXYZ(ix,ny-1-iy,iz))*in[3][IXYZ(ix,ny-1-iy,iz)]/dt);
        }
    }

    // bottom left middle
    #pragma omp parallel for private(sumx,sumy,sumz,sumz0)
    for (int iy=iyminb; iy<iymaxb; iy++){
        for (int ix=ixmin; ix<nc1; ix++){

            // (Dx.in1)_0 and (Dy.in2)_0
            int iz=nz-1;
            sumx=0;
            for (int j=0; j<nc2; j++){
                sumx += bnd_coef[ix*nc2+j] * in[0][IX(j)];
            }
            sumy = coef[0] * (in[1][IY(iy+1)]-in[1][IY(iy-1)]) + coef[1] * (in[1][IY(iy+2)]-in[1][IY(iy-2)]);

            // (Sz.in3)_0
            sumz = 0;
            for (iz = 0; iz < 4; iz++){
                sumz += scoef[iz] * in[2][IZ(nz-1-iz)];
            }
            // (Dz.in3)_0
            sumz0 = 0;
            for (iz = 0; iz < 6; iz++){
                sumz0 -= bnd_coef[iz] * in[2][IZ(nz-1-iz)];
            }
            iz=nz-1;
            out[IZ(iz)] = add*out[IZ(iz)] - adh * (f1(par,IZ(iz))*sumx/dx +f2(par,IZ(iz))*sumy/dy +f3(par,IZ(iz))*sumz0/dz + f4(par,IZ(iz))*sumz/dz -f5(par,IZ(iz))*in[3][IZ(iz)]/dt);
        }
    }

    // bottom right front
    #pragma omp parallel for private(sumx,sumy,sumz,sumz0)
    for (int iy=iymin; iy<nc1; iy++){
        for (int ix=nx-ixmax; ix<std::min(nc1,nx-ixmin); ix++){

            // (Dx.in1)_0 and (Dy.in2)_0
            int iz=nz-1;
            sumx=0;
            sumy=0;
            for (int j=0; j<nc2; j++){
                sumx -= bnd_coef[ix*nc2+j] * in[0][IX(nx-1-j)];
                sumy += bnd_coef[iy*nc2+j] * in[1][IXYZ(nx-1-ix,j,iz)];
            }

            // (Sz.in3)_0
            sumz = 0;
            for (iz = 0; iz < 4; iz++){
                sumz += scoef[iz] * in[2][IXYZ(nx-1-ix,iy,nz-1-iz)];
            }
            // (Dz.in3)_0
            sumz0 = 0;
            for (iz = 0; iz < 6; iz++){
                sumz0 -= bnd_coef[iz] * in[2][IXYZ(nx-1-ix,iy,nz-1-iz)];
            }
            iz=nz-1;
            out[IXYZ(nx-1-ix,iy,iz)] = add*out[IXYZ(nx-1-ix,iy,iz)] - adh * (f1(par,IXYZ(nx-1-ix,iy,iz))*sumx/dx +f2(par,IXYZ(nx-1-ix,iy,iz))*sumy/dy +f3(par,IXYZ(nx-1-ix,iy,iz))*sumz0/dz + f4(par,IXYZ(nx-1-ix,iy,iz))*sumz/dz -f5(par,IXYZ(nx-1-ix,iy,iz))*in[3][IXYZ(nx-1-ix,iy,iz)]/dt);
        }
    }

    // bottom right back
    #pragma omp parallel for private(sumx,sumy,sumz,sumz0)
    for (int iy=ny-iymax; iy<std::min(nc1,ny-iymin); iy++){
        for (int ix=nx-ixmax; ix<std::min(nc1,nx-ixmin); ix++){

            // (Dx.in1)_0 and (Dy.in2)_0
            int iz=nz-1;
            sumx=0;
            sumy=0;
            for (int j=0; j<nc2; j++){
                sumx -= bnd_coef[ix*nc2+j] * in[0][IXYZ(nx-1-j,ny-1-iy,iz)];
                sumy -= bnd_coef[iy*nc2+j] * in[1][IXYZ(nx-1-ix,ny-1-j,iz)];
            }

            // (Sz.in3)_0
            sumz = 0;
            for (iz = 0; iz < 4; iz++){
                sumz += scoef[iz] * in[2][IXYZ(nx-1-ix,ny-1-iy,nz-1-iz)];
            }
            // (Dz.in3)_0
            sumz0 = 0;
            for (iz = 0; iz < 6; iz++){
                sumz0 -= bnd_coef[iz] * in[2][IXYZ(nx-1-ix,ny-1-iy,nz-1-iz)];
            }
            iz=nz-1;
            out[IXYZ(nx-1-ix,ny-1-iy,iz)] = add*out[IXYZ(nx-1-ix,ny-1-iy,iz)] - adh * (f1(par,IXYZ(nx-1-ix,ny-1-iy,iz))*sumx/dx +f2(par,IXYZ(nx-1-ix,ny-1-iy,iz))*sumy/dy +f3(par,IXYZ(nx-1-ix,ny-1-iy,iz))*sumz0/dz + f4(par,IXYZ(nx-1-ix,ny-1-iy,iz))*sumz/dz -f5(par,IXYZ(nx-1-ix,ny-1-iy,iz))*in[3][IXYZ(nx-1-ix,ny-1-iy,iz)]/dt);
        }
    }

    // bottom right middle
    #pragma omp parallel for private(sumx,sumy,sumz,sumz0)
    for (int iy=iyminb; iy<iymaxb; iy++){
        for (int ix=nx-ixmax; ix<std::min(nc1,nx-ixmin); ix++){

            // (Dx.in1)_0 and (Dy.in2)_0
            int iz=nz-1;
            sumx=0;
            for (int j=0; j<nc2; j++){
                sumx -= bnd_coef[ix*nc2+j] * in[0][IX(nx-1-j)];
            }
            sumy = coef[0] * (in[1][IXYZ(nx-1-ix,iy+1,iz)]-in[1][IXYZ(nx-1-ix,iy-1,iz)]) + coef[1] * (in[1][IXYZ(nx-1-ix,iy+2,iz)]-in[1][IXYZ(nx-1-ix,iy-2,iz)]);

            // (Sz.in3)_0
            sumz = 0;
            for (iz = 0; iz < 4; iz++){
                sumz += scoef[iz] * in[2][IXYZ(nx-1-ix,iy,nz-1-iz)];
            }
            // (Dz.in3)_0
            sumz0 = 0;
            for (iz = 0; iz < 6; iz++){
                sumz0 -= bnd_coef[iz] * in[2][IXYZ(nx-1-ix,iy,nz-1-iz)];
            }
            iz=nz-1;
            out[IXYZ(nx-1-ix,iy,iz)] = add*out[IXYZ(nx-1-ix,iy,iz)] - adh * (f1(par,IXYZ(nx-1-ix,iy,iz))*sumx/dx +f2(par,IXYZ(nx-1-ix,iy,iz))*sumy/dy +f3(par,IXYZ(nx-1-ix,iy,iz))*sumz0/dz + f4(par,IXYZ(nx-1-ix,iy,iz))*sumz/dz -f5(par,IXYZ(nx-1-ix,iy,iz))*in[3][IXYZ(nx-1-ix,iy,iz)]/dt);
        }
    }

    // bottom middle front
    #pragma omp parallel for private(sumx,sumy,sumz,sumz0)
    for (int iy=iymin; iy<nc1; iy++){
        for (int ix=ixminb; ix<ixmaxb; ix++){

            // (Dx.in1)_0 and (Dy.in2)_0
            int iz=nz-1;
            sumx = coef[0] * (in[0][IX(ix+1)]-in[0][IX(ix-1)]) + coef[1] * (in[0][IX(ix+2)]-in[0][IX(ix-2)]);
            sumy=0;
            for (int j=0; j<nc2; j++){
                sumy += bnd_coef[iy*nc2+j] * in[1][IY(j)];
            }

            // (Sz.in3)_0
            sumz = 0;
            for (iz = 0; iz < 4; iz++){
                sumz += scoef[iz] * in[2][IZ(nz-1-iz)];
            }
            // (Dz.in3)_0
            sumz0 = 0;
            for (iz = 0; iz < 6; iz++){
                sumz0 -= bnd_coef[iz] * in[2][IZ(nz-1-iz)];
            }
            iz=nz-1;
            out[IZ(iz)] = add*out[IZ(iz)] - adh * (f1(par,IZ(iz))*sumx/dx +f2(par,IZ(iz))*sumy/dy +f3(par,IZ(iz))*sumz0/dz + f4(par,IZ(iz))*sumz/dz -f5(par,IZ(iz))*in[3][IZ(iz)]/dt);
        }
    }

    // bottom middle back
    #pragma omp parallel for private(sumx,sumy,sumz,sumz0)
    for (int iy=ny-iymax; iy<std::min(nc1,ny-iymin); iy++){
        for (int ix=ixminb; ix<ixmaxb; ix++){

            // (Dx.in1)_0 and (Dy.in2)_0
            int iz=nz-1;
            sumx = coef[0] * (in[0][IXYZ(ix+1,ny-1-iy,iz)]-in[0][IXYZ(ix-1,ny-1-iy,iz)]) + coef[1] * (in[0][IXYZ(ix+2,ny-1-iy,iz)]-in[0][IXYZ(ix-2,ny-1-iy,iz)]);
            sumy=0;
            for (int j=0; j<nc2; j++){
                sumy -= bnd_coef[iy*nc2+j] * in[1][IY(ny-1-j)];
            }

            // (Sz.in3)_0
            sumz = 0;
            for (iz = 0; iz < 4; iz++){
                sumz += scoef[iz] * in[2][IXYZ(ix,ny-1-iy,nz-1-iz)];
            }
            // (Dz.in3)_0
            sumz0 = 0;
            for (iz = 0; iz < 6; iz++){
                sumz0 -= bnd_coef[iz] * in[2][IXYZ(ix,ny-1-iy,nz-1-iz)];
            }
            iz=nz-1;
            out[IXYZ(ix,ny-1-iy,iz)] = add*out[IXYZ(ix,ny-1-iy,iz)] - adh * (f1(par,IXYZ(ix,ny-1-iy,iz))*sumx/dx +f2(par,IXYZ(ix,ny-1-iy,iz))*sumy/dy +f3(par,IXYZ(ix,ny-1-iy,iz))*sumz0/dz + f4(par,IXYZ(ix,ny-1-iy,iz))*sumz/dz -f5(par,IXYZ(ix,ny-1-iy,iz))*in[3][IXYZ(ix,ny-1-iy,iz)]/dt);
        }
    }

    // bottom middle
    #pragma omp parallel for private(sumx, sumy, sumz,sumz0)
    for (int iy=iyminb; iy<iymaxb; iy++){
        for (int ix=ixminb; ix<ixmaxb; ix++){

            // (Dx.in1)_0 and (Dy.in2)_0
            int iz=nz-1;
            sumx = coef[0] * (in[0][IX(ix+1)]-in[0][IX(ix-1)]) + coef[1] * (in[0][IX(ix+2)]-in[0][IX(ix-2)]);
            sumy = coef[0] * (in[1][IY(iy+1)]-in[1][IY(iy-1)]) + coef[1] * (in[1][IY(iy+2)]-in[1][IY(iy-2)]);

            // (Sz.in3)_0
            sumz = 0;
            for (iz = 0; iz < 4; iz++){
                sumz += scoef[iz] * in[2][IZ(nz-1-iz)];
            }
            // (Dz.in3)_0
            sumz0 = 0;
            for (iz = 0; iz < 6; iz++){
                sumz0 -= bnd_coef[iz] * in[2][IZ(nz-1-iz)];
            }
            iz=nz-1;
            out[IZ(iz)] = add*out[IZ(iz)] - adh * (f1(par,IZ(iz))*sumx/dx +f2(par,IZ(iz))*sumy/dy +f3(par,IZ(iz))*sumz0/dz + f4(par,IZ(iz))*sumz/dz -f5(par,IZ(iz))*in[3][IZ(iz)]/dt);
        }
    }
}

template<expr f1, expr f2, expr f3, expr f4, expr f5>
void esat_neumann_left(bool add, const data_t** in, __restrict data_t* out, int nx, int ny, int nz, data_t dx, data_t dy, data_t dz, data_t dt, int iymin, int iymax, int izmin, int izmax, const data_t ** par){

    data_t coef[2] = {2.0/3,-1.0/12};
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3}; 
    data_t h0 = 17.0/48;
    int nc1=4, nc2=6;
    int iyminb=std::max(nc1,iymin);
    int iymaxb=std::min(ny - nc1, iymax);
    int izminb=std::max(nc1,izmin);
    int izmaxb=std::min(nz - nc1, izmax);
    data_t adh = 1.0/(dx*h0);

    data_t sumx, sumy, sumz, sumx0;

    // SAT = - Hx-1.(-f1.Dy.in1 -f2.Dz.in2 -f3.Dx.in3 + f4.Sx.in3 - f5.in4/dt)_0
    // Sx is the boundary derivative operator pointing outwards

    // left top front
    #pragma omp parallel for private(sumx,sumy,sumz,sumx0)
    for (int iy=iymin; iy<nc1; iy++){
        for (int iz=izmin; iz<nc1; iz++){

            // (Dy.in1)_0 and (Dz.in2)_0
            int ix=0;
            sumy=0;
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumy += bnd_coef[iy*nc2+j] * in[0][IY(j)];
                sumz += bnd_coef[iz*nc2+j] * in[1][IZ(j)];
            }

            // (Sx.in3)_0
            sumx = 0;
            for (ix = 0; ix < 4; ix++){
                sumx += scoef[ix] * in[2][IX(ix)];
            }
            // (Dx.in3)_0
            sumx0 = 0;
            for (ix = 0; ix < 6; ix++){
                sumx0 += bnd_coef[ix] * in[2][IX(ix)];
            }
            ix=0;
            out[IX(ix)] = add*out[IX(ix)] - adh * (-f1(par,IX(ix))*sumy/dy -f2(par,IX(ix))*sumz/dz -f3(par,IX(ix))*sumx0/dx + f4(par,IX(ix))*sumx/dx -f5(par,IX(ix))*in[3][IX(ix)]/dt);
        }
    }

    // left top back
    #pragma omp parallel for private(sumx,sumy,sumz,sumx0)
    for (int iy=ny-iymax; iy<std::min(nc1,ny-iymin); iy++){
        for (int iz=izmin; iz<nc1; iz++){

            // (Dy.in1)_0 and (Dz.in2)_0
            int ix=0;
            sumy=0;
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumy -= bnd_coef[iy*nc2+j] * in[0][IY(ny-1-j)];
                sumz += bnd_coef[iz*nc2+j] * in[1][IXYZ(ix,ny-1-iy,j)];
            }

            // (Sx.in3)_0
            sumx = 0;
            for (ix = 0; ix < 4; ix++){
                sumx += scoef[ix] * in[2][IXYZ(ix,ny-1-iy,iz)];
            }
            // (Dx.in3)_0
            sumx0 = 0;
            for (ix = 0; ix < 6; ix++){
                sumx0 += bnd_coef[ix] * in[2][IXYZ(ix,ny-1-iy,iz)];
            }
            ix=0;
            out[IXYZ(ix,ny-1-iy,iz)] = add*out[IXYZ(ix,ny-1-iy,iz)] - adh * (-f1(par,IXYZ(ix,ny-1-iy,iz))*sumy/dy -f2(par,IXYZ(ix,ny-1-iy,iz))*sumz/dz -f3(par,IXYZ(ix,ny-1-iy,iz))*sumx0/dx + f4(par,IXYZ(ix,ny-1-iy,iz))*sumx/dx -f5(par,IXYZ(ix,ny-1-iy,iz))*in[3][IXYZ(ix,ny-1-iy,iz)]/dt);
        }
    }

    // left top middle
    #pragma omp parallel for private(sumx,sumy,sumz,sumx0)
    for (int iy=iyminb; iy<iymaxb; iy++){
        for (int iz=izmin; iz<nc1; iz++){

            // (Dy.in1)_0 and (Dz.in2)_0
            int ix=0;
            sumy = coef[0] * (in[0][IY(iy+1)]-in[0][IY(iy-1)]) + coef[1] * (in[0][IY(iy+2)]-in[0][IY(iy-2)]);
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumz += bnd_coef[iz*nc2+j] * in[1][IZ(j)];
            }

            // (Sx.in3)_0
            sumx = 0;
            for (ix = 0; ix < 4; ix++){
                sumx += scoef[ix] * in[2][IX(ix)];
            }
            // (Dx.in3)_0
            sumx0 = 0;
            for (ix = 0; ix < 6; ix++){
                sumx0 += bnd_coef[ix] * in[2][IX(ix)];
            }
            ix=0;
            out[IX(ix)] = add*out[IX(ix)] - adh * (-f1(par,IX(ix))*sumy/dy -f2(par,IX(ix))*sumz/dz -f3(par,IX(ix))*sumx0/dx + f4(par,IX(ix))*sumx/dx -f5(par,IX(ix))*in[3][IX(ix)]/dt);;
        }
    }

    // left bottom front
    #pragma omp parallel for private(sumx,sumy,sumz,sumx0)
    for (int iy=iymin; iy<nc1; iy++){
        for (int iz=nz-izmax; iz<std::min(nc1,nz-izmin); iz++){

            // (Dy.in1)_0 and (Dz.in2)_0
            int ix=0;
            sumy=0;
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumy += bnd_coef[iy*nc2+j] * in[0][IXYZ(ix,j,nz-1-iz)];
                sumz -= bnd_coef[iz*nc2+j] * in[1][IZ(nz-1-j)];
            }

            // (Sx.in3)_0
            sumx = 0;
            for (ix = 0; ix < 4; ix++){
                sumx += scoef[ix] * in[2][IXYZ(ix,iy,nz-1-iz)];
            }
            // (Dx.in3)_0
            sumx0 = 0;
            for (ix = 0; ix < 6; ix++){
                sumx0 += bnd_coef[ix] * in[2][IXYZ(ix,iy,nz-1-iz)];
            }
            ix=0;
            out[IXYZ(ix,iy,nz-1-iz)] = add*out[IXYZ(ix,iy,nz-1-iz)] - adh * (-f1(par,IXYZ(ix,iy,nz-1-iz))*sumy/dy -f2(par,IXYZ(ix,iy,nz-1-iz))*sumz/dz -f3(par,IXYZ(ix,iy,nz-1-iz))*sumx0/dx + f4(par,IXYZ(ix,iy,nz-1-iz))*sumx/dx -f5(par,IXYZ(ix,iy,nz-1-iz))*in[3][IXYZ(ix,iy,nz-1-iz)]/dt);;
        }
    }

    // left bottom back
    #pragma omp parallel for private(sumx,sumy,sumz,sumx0)
    for (int iy=ny-iymax; iy<std::min(nc1,ny-iymin); iy++){
        for (int iz=nz-izmax; iz<std::min(nc1,nz-izmin); iz++){

            // (Dy.in1)_0 and (Dz.in2)_0
            int ix=0;
            sumy=0;
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumy -= bnd_coef[iy*nc2+j] * in[0][IXYZ(ix,ny-1-j,nz-1-iz)];
                sumz -= bnd_coef[iz*nc2+j] * in[1][IXYZ(ix,ny-1-iy,nz-1-j)];
            }

            // (Sx.in3)_0
            sumx = 0;
            for (ix = 0; ix < 4; ix++){
                sumx += scoef[ix] * in[2][IXYZ(ix,ny-1-iy,nz-1-iz)];
            }
            // (Dx.in3)_0
            sumx0 = 0;
            for (ix = 0; ix < 6; ix++){
                sumx0 += bnd_coef[ix] * in[2][IXYZ(ix,ny-1-iy,nz-1-iz)];
            }
            ix=0;
            out[IXYZ(ix,ny-1-iy,nz-1-iz)] = add*out[IXYZ(ix,ny-1-iy,nz-1-iz)] - adh * (-f1(par,IXYZ(ix,ny-1-iy,nz-1-iz))*sumy/dy -f2(par,IXYZ(ix,ny-1-iy,nz-1-iz))*sumz/dz -f3(par,IXYZ(ix,ny-1-iy,nz-1-iz))*sumx0/dx + f4(par,IXYZ(ix,ny-1-iy,nz-1-iz))*sumx/dx -f5(par,IXYZ(ix,ny-1-iy,nz-1-iz))*in[3][IXYZ(ix,ny-1-iy,nz-1-iz)]/dt);;
        }
    }

    // left bottom middle
    #pragma omp parallel for private(sumx,sumy,sumz,sumx0)
    for (int iy=iyminb; iy<iymaxb; iy++){
        for (int iz=nz-izmax; iz<std::min(nc1,nz-izmin); iz++){

            // (Dy.in1)_0 and (Dz.in2)_0
            int ix=0;
            sumy = coef[0] * (in[0][IXYZ(ix,iy+1,nz-1-iz)]-in[0][IXYZ(ix,iy-1,nz-1-iz)]) + coef[1] * (in[0][IXYZ(ix,iy+2,nz-1-iz)]-in[0][IXYZ(ix,iy-2,nz-1-iz)]);
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumz -= bnd_coef[iz*nc2+j] * in[1][IZ(nz-1-j)];
            }

            // (Sx.in3)_0
            sumx = 0;
            for (ix = 0; ix < 4; ix++){
                sumx += scoef[ix] * in[2][IXYZ(ix,iy,nz-1-iz)];
            }
            // (Dx.in3)_0
            sumx0 = 0;
            for (ix = 0; ix < 6; ix++){
                sumx0 += bnd_coef[ix] * in[2][IXYZ(ix,iy,nz-1-iz)];
            }
            ix=0;
            out[IXYZ(ix,iy,nz-1-iz)] = add*out[IXYZ(ix,iy,nz-1-iz)] - adh * (-f1(par,IXYZ(ix,iy,nz-1-iz))*sumy/dy -f2(par,IXYZ(ix,iy,nz-1-iz))*sumz/dz -f3(par,IXYZ(ix,iy,nz-1-iz))*sumx0/dx + f4(par,IXYZ(ix,iy,nz-1-iz))*sumx/dx -f5(par,IXYZ(ix,iy,nz-1-iz))*in[3][IXYZ(ix,iy,nz-1-iz)]/dt);;
        }
    }

    // left middle front
    #pragma omp parallel for private(sumx,sumy,sumz,sumx0)
    for (int iy=iymin; iy<nc1; iy++){
        for (int iz=izminb; iz<izmaxb; iz++){

            // (Dy.in1)_0 and (Dz.in2)_0
            int ix=0;
            sumy=0;
            for (int j=0; j<nc2; j++){
                sumy += bnd_coef[iy*nc2+j] * in[0][IY(j)];
            }
            sumz = coef[0] * (in[1][IZ(iz+1)]-in[1][IZ(iz-1)]) + coef[1] * (in[1][IZ(iz+2)]-in[1][IZ(iz-2)]);

            // (Sx.in3)_0
            sumx = 0;
            for (ix = 0; ix < 4; ix++){
                sumx += scoef[ix] * in[2][IX(ix)];
            }
            // (Dx.in3)_0
            sumx0 = 0;
            for (ix = 0; ix < 6; ix++){
                sumx0 += bnd_coef[ix] * in[2][IX(ix)];
            }
            ix=0;
            out[IX(ix)] = add*out[IX(ix)] - adh * (-f1(par,IX(ix))*sumy/dy -f2(par,IX(ix))*sumz/dz -f3(par,IX(ix))*sumx0/dx + f4(par,IX(ix))*sumx/dx -f5(par,IX(ix))*in[3][IX(ix)]/dt);;
        }
    }

    // left middle back
    #pragma omp parallel for private(sumx,sumy,sumz,sumx0)
    for (int iy=ny-iymax; iy<std::min(nc1,ny-iymin); iy++){
        for (int iz=izminb; iz<izmaxb; iz++){

            // (Dy.in1)_0 and (Dz.in2)_0
            int ix=0;
            sumy=0;
            for (int j=0; j<nc2; j++){
                sumy -= bnd_coef[iy*nc2+j] * in[0][IY(ny-1-j)];
            }
            sumz = coef[0] * (in[1][IXYZ(ix,ny-1-iy,iz+1)]-in[1][IXYZ(ix,ny-1-iy,iz-1)]) + coef[1] * (in[1][IXYZ(ix,ny-1-iy,iz+2)]-in[1][IXYZ(ix,ny-1-iy,iz-2)]);

            // (Sx.in3)_0
            sumx = 0;
            for (ix = 0; ix < 4; ix++){
                sumx += scoef[ix] * in[2][IXYZ(ix,ny-1-iy,iz)];
            }
            // (Dx.in3)_0
            sumx0 = 0;
            for (ix = 0; ix < 6; ix++){
                sumx0 += bnd_coef[ix] * in[2][IXYZ(ix,ny-1-iy,iz)];
            }
            ix=0;
            out[IXYZ(ix,ny-1-iy,iz)] = add*out[IXYZ(ix,ny-1-iy,iz)] - adh * (-f1(par,IXYZ(ix,ny-1-iy,iz))*sumy/dy -f2(par,IXYZ(ix,ny-1-iy,iz))*sumz/dz -f3(par,IXYZ(ix,ny-1-iy,iz))*sumx0/dx + f4(par,IXYZ(ix,ny-1-iy,iz))*sumx/dx -f5(par,IXYZ(ix,ny-1-iy,iz))*in[3][IXYZ(ix,ny-1-iy,iz)]/dt);;
        }
    }

    // left middle
    #pragma omp parallel for private(sumx, sumy, sumz)
    for (int iy=iyminb; iy<iymaxb; iy++){
        for (int iz=izminb; iz<izmaxb; iz++){

            // (Dy.in1)_0 and (Dz.in2)_0
            int ix=0;
            sumy = coef[0] * (in[0][IY(iy+1)]-in[0][IY(iy-1)]) + coef[1] * (in[0][IY(iy+2)]-in[0][IY(iy-2)]);
            sumz = coef[0] * (in[1][IZ(iz+1)]-in[1][IZ(iz-1)]) + coef[1] * (in[1][IZ(iz+2)]-in[1][IZ(iz-2)]);

            // (Sx.in3)_0
            sumx = 0;
            for (ix = 0; ix < 4; ix++){
                sumx += scoef[ix] * in[2][IX(ix)];
            }
            // (Dx.in3)_0
            sumx0 = 0;
            for (ix = 0; ix < 6; ix++){
                sumx0 += bnd_coef[ix] * in[2][IX(ix)];
            }
            ix=0;
            out[IX(ix)] = add*out[IX(ix)] - adh * (-f1(par,IX(ix))*sumy/dy -f2(par,IX(ix))*sumz/dz -f3(par,IX(ix))*sumx0/dx + f4(par,IX(ix))*sumx/dx -f5(par,IX(ix))*in[3][IX(ix)]/dt);;
        }
    }
}

template<expr f1, expr f2, expr f3, expr f4, expr f5>
void esat_neumann_right(bool add, const data_t** in, __restrict data_t* out, int nx, int ny, int nz, data_t dx, data_t dy, data_t dz, data_t dt, int iymin, int iymax, int izmin, int izmax, const data_t ** par){

    data_t coef[2] = {2.0/3,-1.0/12};
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3}; 
    data_t h0 = 17.0/48;
    int nc1=4, nc2=6;
    int iyminb=std::max(nc1,iymin);
    int iymaxb=std::min(ny - nc1, iymax);
    int izminb=std::max(nc1,izmin);
    int izmaxb=std::min(nz - nc1, izmax);
    data_t adh = 1.0/(dx*h0);

    data_t sumx, sumy, sumz, sumx0;

    // SAT = - Hx-1.(f1.Dy.in1 +f2.Dz.in2 + f3.Dx.in3 + f4.Sx.in3 - f5.in4/dt)_0
    // Sx is the boundary derivative operator pointing outwards

    // right top front
    #pragma omp parallel for private(sumx,sumy,sumz,sumx0)
    for (int iy=iymin; iy<nc1; iy++){
        for (int iz=izmin; iz<nc1; iz++){

            // (Dy.in1)_0 and (Dz.in2)_0
            int ix=nx-1;
            sumy=0;
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumy += bnd_coef[iy*nc2+j] * in[0][IY(j)];
                sumz += bnd_coef[iz*nc2+j] * in[1][IZ(j)];
            }

            // (Sx.in3)_0
            sumx = 0;
            for (ix = 0; ix < 4; ix++){
                sumx += scoef[ix] * in[2][IX(nx-1-ix)];
            }
            // (Dx.in3)_0
            sumx0 = 0;
            for (ix = 0; ix < 6; ix++){
                sumx0 -= bnd_coef[ix] * in[2][IX(nx-1-ix)];
            }
            ix=nx-1;
            out[IX(ix)] = add*out[IX(ix)] - adh * (f1(par,IX(ix))*sumy/dy +f2(par,IX(ix))*sumz/dz +f3(par,IX(ix))*sumx0/dx + f4(par,IX(ix))*sumx/dx -f5(par,IX(ix))*in[3][IX(ix)]/dt);
        }
    }

    // right top back
    #pragma omp parallel for private(sumx,sumy,sumz,sumx0)
    for (int iy=ny-iymax; iy<std::min(nc1,ny-iymin); iy++){
        for (int iz=izmin; iz<nc1; iz++){

            // (Dy.in1)_0 and (Dz.in2)_0
            int ix=nx-1;
            sumy=0;
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumy -= bnd_coef[iy*nc2+j] * in[0][IY(ny-1-j)];
                sumz += bnd_coef[iz*nc2+j] * in[1][IXYZ(ix,ny-1-iy,j)];
            }

            // (Sx.in3)_0
            sumx = 0;
            for (ix = 0; ix < 4; ix++){
                sumx += scoef[ix] * in[2][IXYZ(nx-1-ix,ny-1-iy,iz)];
            }
            // (Dx.in3)_0
            sumx0 = 0;
            for (ix = 0; ix < 6; ix++){
                sumx0 -= bnd_coef[ix] * in[2][IXYZ(nx-1-ix,ny-1-iy,iz)];
            }
            ix=nx-1;
            out[IXYZ(ix,ny-1-iy,iz)] = add*out[IXYZ(ix,ny-1-iy,iz)] - adh * (f1(par,IXYZ(ix,ny-1-iy,iz))*sumy/dy +f2(par,IXYZ(ix,ny-1-iy,iz))*sumz/dz +f3(par,IXYZ(ix,ny-1-iy,iz))*sumx0/dx + f4(par,IXYZ(ix,ny-1-iy,iz))*sumx/dx -f5(par,IXYZ(ix,ny-1-iy,iz))*in[3][IXYZ(ix,ny-1-iy,iz)]/dt);
        }
    }

    // right top middle
    #pragma omp parallel for private(sumx,sumy,sumz,sumx0)
    for (int iy=iyminb; iy<iymaxb; iy++){
        for (int iz=izmin; iz<nc1; iz++){

            // (Dy.in1)_0 and (Dz.in2)_0
            int ix=nx-1;
            sumy = coef[0] * (in[0][IY(iy+1)]-in[0][IY(iy-1)]) + coef[1] * (in[0][IY(iy+2)]-in[0][IY(iy-2)]);
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumz += bnd_coef[iz*nc2+j] * in[1][IZ(j)];
            }

            // (Sx.in3)_0
            sumx = 0;
            for (ix = 0; ix < 4; ix++){
                sumx += scoef[ix] * in[2][IX(nx-1-ix)];
            }
            // (Dx.in3)_0
            sumx0 = 0;
            for (ix = 0; ix < 6; ix++){
                sumx0 -= bnd_coef[ix] * in[2][IX(nx-1-ix)];
            }
            ix=nx-1;
            out[IX(ix)] = add*out[IX(ix)] - adh * (f1(par,IX(ix))*sumy/dy +f2(par,IX(ix))*sumz/dz +f3(par,IX(ix))*sumx0/dx + f4(par,IX(ix))*sumx/dx -f5(par,IX(ix))*in[3][IX(ix)]/dt);
        }
    }

    // right bottom front
    #pragma omp parallel for private(sumx,sumy,sumz,sumx0)
    for (int iy=iymin; iy<nc1; iy++){
        for (int iz=nz-izmax; iz<std::min(nc1,nz-izmin); iz++){

            // (Dy.in1)_0 and (Dz.in2)_0
            int ix=nx-1;
            sumy=0;
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumy += bnd_coef[iy*nc2+j] * in[0][IXYZ(ix,j,nz-1-iz)];
                sumz -= bnd_coef[iz*nc2+j] * in[1][IZ(nz-1-j)];
            }

            // (Sx.in3)_0
            sumx = 0;
            for (ix = 0; ix < 4; ix++){
                sumx += scoef[ix] * in[2][IXYZ(nx-1-ix,iy,nz-1-iz)];
            }
            // (Dx.in3)_0
            sumx0 = 0;
            for (ix = 0; ix < 6; ix++){
                sumx0 -= bnd_coef[ix] * in[2][IXYZ(nx-1-ix,iy,nz-1-iz)];
            }
            ix=nx-1;
            out[IXYZ(ix,iy,nz-1-iz)] = add*out[IXYZ(ix,iy,nz-1-iz)] - adh * (f1(par,IXYZ(ix,iy,nz-1-iz))*sumy/dy +f2(par,IXYZ(ix,iy,nz-1-iz))*sumz/dz +f3(par,IXYZ(ix,iy,nz-1-iz))*sumx0/dx + f4(par,IXYZ(ix,iy,nz-1-iz))*sumx/dx -f5(par,IXYZ(ix,iy,nz-1-iz))*in[3][IXYZ(ix,iy,nz-1-iz)]/dt);
        }
    }

    // right bottom back
    #pragma omp parallel for private(sumx,sumy,sumz,sumx0)
    for (int iy=ny-iymax; iy<std::min(nc1,ny-iymin); iy++){
        for (int iz=nz-izmax; iz<std::min(nc1,nz-izmin); iz++){

            // (Dy.in1)_0 and (Dz.in2)_0
            int ix=nx-1;
            sumy=0;
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumy -= bnd_coef[iy*nc2+j] * in[0][IXYZ(ix,ny-1-j,nz-1-iz)];
                sumz -= bnd_coef[iz*nc2+j] * in[1][IXYZ(ix,ny-1-iy,nz-1-j)];
            }

            // (Sx.in3)_0
            sumx = 0;
            for (ix = 0; ix < 4; ix++){
                sumx += scoef[ix] * in[2][IXYZ(nx-1-ix,ny-1-iy,nz-1-iz)];
            }
            // (Dx.in3)_0
            sumx0 = 0;
            for (ix = 0; ix < 6; ix++){
                sumx0 -= bnd_coef[ix] * in[2][IXYZ(nx-1-ix,ny-1-iy,nz-1-iz)];
            }
            ix=nx-1;
            out[IXYZ(ix,ny-1-iy,nz-1-iz)] = add*out[IXYZ(ix,ny-1-iy,nz-1-iz)] - adh * (f1(par,IXYZ(ix,ny-1-iy,nz-1-iz))*sumy/dy +f2(par,IXYZ(ix,ny-1-iy,nz-1-iz))*sumz/dz +f3(par,IXYZ(ix,ny-1-iy,nz-1-iz))*sumx0/dx + f4(par,IXYZ(ix,ny-1-iy,nz-1-iz))*sumx/dx -f5(par,IXYZ(ix,ny-1-iy,nz-1-iz))*in[3][IXYZ(ix,ny-1-iy,nz-1-iz)]/dt);
        }
    }

    // right bottom middle
    #pragma omp parallel for private(sumx,sumy,sumz,sumx0)
    for (int iy=iyminb; iy<iymaxb; iy++){
        for (int iz=nz-izmax; iz<std::min(nc1,nz-izmin); iz++){

            // (Dy.in1)_0 and (Dz.in2)_0
            int ix=nx-1;
            sumy = coef[0] * (in[0][IXYZ(ix,iy+1,nz-1-iz)]-in[0][IXYZ(ix,iy-1,nz-1-iz)]) + coef[1] * (in[0][IXYZ(ix,iy+2,nz-1-iz)]-in[0][IXYZ(ix,iy-2,nz-1-iz)]);
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumz -= bnd_coef[iz*nc2+j] * in[1][IZ(nz-1-j)];
            }

            // (Sx.in3)_0
            sumx = 0;
            for (ix = 0; ix < 4; ix++){
                sumx += scoef[ix] * in[2][IXYZ(nx-1-ix,iy,nz-1-iz)];
            }
            // (Dx.in3)_0
            sumx0 = 0;
            for (ix = 0; ix < 6; ix++){
                sumx0 -= bnd_coef[ix] * in[2][IXYZ(nx-1-ix,iy,nz-1-iz)];
            }
            ix=nx-1;
            out[IXYZ(ix,iy,nz-1-iz)] = add*out[IXYZ(ix,iy,nz-1-iz)] - adh * (f1(par,IXYZ(ix,iy,nz-1-iz))*sumy/dy +f2(par,IXYZ(ix,iy,nz-1-iz))*sumz/dz +f3(par,IXYZ(ix,iy,nz-1-iz))*sumx0/dx + f4(par,IXYZ(ix,iy,nz-1-iz))*sumx/dx -f5(par,IXYZ(ix,iy,nz-1-iz))*in[3][IXYZ(ix,iy,nz-1-iz)]/dt);
        }
    }

    // right middle front
    #pragma omp parallel for private(sumx,sumy,sumz,sumx0)
    for (int iy=iymin; iy<nc1; iy++){
        for (int iz=izminb; iz<izmaxb; iz++){

            // (Dy.in1)_0 and (Dz.in2)_0
            int ix=nx-1;
            sumy=0;
            for (int j=0; j<nc2; j++){
                sumy += bnd_coef[iy*nc2+j] * in[0][IY(j)];
            }
            sumz = coef[0] * (in[1][IZ(iz+1)]-in[1][IZ(iz-1)]) + coef[1] * (in[1][IZ(iz+2)]-in[1][IZ(iz-2)]);

            // (Sx.in3)_0
            sumx = 0;
            for (ix = 0; ix < 4; ix++){
                sumx += scoef[ix] * in[2][IX(nx-1-ix)];
            }
            // (Dx.in3)_0
            sumx0 = 0;
            for (ix = 0; ix < 6; ix++){
                sumx0 -= bnd_coef[ix] * in[2][IX(nx-1-ix)];
            }
            ix=nx-1;
            out[IX(ix)] = add*out[IX(ix)] - adh * (f1(par,IX(ix))*sumy/dy +f2(par,IX(ix))*sumz/dz +f3(par,IX(ix))*sumx0/dx + f4(par,IX(ix))*sumx/dx -f5(par,IX(ix))*in[3][IX(ix)]/dt);
        }
    }

    // right middle back
    #pragma omp parallel for private(sumx,sumy,sumz,sumx0)
    for (int iy=ny-iymax; iy<std::min(nc1,ny-iymin); iy++){
        for (int iz=izminb; iz<izmaxb; iz++){

            // (Dy.in1)_0 and (Dz.in2)_0
            int ix=nx-1;
            sumy=0;
            for (int j=0; j<nc2; j++){
                sumy -= bnd_coef[iy*nc2+j] * in[0][IY(ny-1-j)];
            }
            sumz = coef[0] * (in[1][IXYZ(ix,ny-1-iy,iz+1)]-in[1][IXYZ(ix,ny-1-iy,iz-1)]) + coef[1] * (in[1][IXYZ(ix,ny-1-iy,iz+2)]-in[1][IXYZ(ix,ny-1-iy,iz-2)]);

            // (Sx.in3)_0
            sumx = 0;
            for (ix = 0; ix < 4; ix++){
                sumx += scoef[ix] * in[2][IXYZ(nx-1-ix,ny-1-iy,iz)];
            }
            // (Dx.in3)_0
            sumx0 = 0;
            for (ix = 0; ix < 6; ix++){
                sumx0 -= bnd_coef[ix] * in[2][IXYZ(nx-1-ix,ny-1-iy,iz)];
            }
            ix=nx-1;
            out[IXYZ(ix,ny-1-iy,iz)] = add*out[IXYZ(ix,ny-1-iy,iz)] - adh * (f1(par,IXYZ(ix,ny-1-iy,iz))*sumy/dy +f2(par,IXYZ(ix,ny-1-iy,iz))*sumz/dz +f3(par,IXYZ(ix,ny-1-iy,iz))*sumx0/dx + f4(par,IXYZ(ix,ny-1-iy,iz))*sumx/dx -f5(par,IXYZ(ix,ny-1-iy,iz))*in[3][IXYZ(ix,ny-1-iy,iz)]/dt);
        }
    }

    // right middle
    #pragma omp parallel for private(sumx, sumy, sumz)
    for (int iy=iyminb; iy<iymaxb; iy++){
        for (int iz=izminb; iz<izmaxb; iz++){

            // (Dy.in1)_0 and (Dz.in2)_0
            int ix=nx-1;
            sumy = coef[0] * (in[0][IY(iy+1)]-in[0][IY(iy-1)]) + coef[1] * (in[0][IY(iy+2)]-in[0][IY(iy-2)]);
            sumz = coef[0] * (in[1][IZ(iz+1)]-in[1][IZ(iz-1)]) + coef[1] * (in[1][IZ(iz+2)]-in[1][IZ(iz-2)]);

            // (Sx.in3)_0
            sumx = 0;
            for (ix = 0; ix < 4; ix++){
                sumx += scoef[ix] * in[2][IX(nx-1-ix)];
            }
            // (Dx.in3)_0
            sumx0 = 0;
            for (ix = 0; ix < 6; ix++){
                sumx0 -= bnd_coef[ix] * in[2][IX(nx-1-ix)];
            }
            ix=nx-1;
            out[IX(ix)] = add*out[IX(ix)] - adh * (f1(par,IX(ix))*sumy/dy +f2(par,IX(ix))*sumz/dz +f3(par,IX(ix))*sumx0/dx + f4(par,IX(ix))*sumx/dx -f5(par,IX(ix))*in[3][IX(ix)]/dt);
        }
    }
}

template<expr f1, expr f2, expr f3, expr f4, expr f5>
void esat_neumann_front(bool add, const data_t** in, __restrict data_t* out, int nx, int ny, int nz, data_t dx, data_t dy, data_t dz, data_t dt, int ixmin, int ixmax, int izmin, int izmax, const data_t ** par){

    data_t coef[2] = {2.0/3,-1.0/12};
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3}; 
    data_t h0 = 17.0/48;
    int nc1=4, nc2=6;
    int ixminb=std::max(nc1,ixmin);
    int ixmaxb=std::min(nx - nc1, ixmax);
    int izminb=std::max(nc1,izmin);
    int izmaxb=std::min(nz - nc1, izmax);
    data_t adh = 1.0/(dy*h0);

    data_t sumx, sumy, sumz, sumy0;

    // SAT = - Hy-1.(-f1.Dx.in1 -f2.Dz.in2 -f3.Dy.in3 + f4.Sy.in3 - f5.in4/dt)_0
    // Sy is the boundary derivative operator pointing outwards

    // front top left
    #pragma omp parallel for private(sumx,sumy,sumz,sumy0)
    for (int ix=ixmin; ix<nc1; ix++){
        for (int iz=izmin; iz<nc1; iz++){

            // (Dx.in1)_0 and (Dz.in2)_0
            int iy=0;
            sumx=0;
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumx += bnd_coef[ix*nc2+j] * in[0][IX(j)];
                sumz += bnd_coef[iz*nc2+j] * in[1][IZ(j)];
            }

            // (Sy.in3)_0
            sumy = 0;
            for (iy = 0; iy < 4; iy++){
                sumy += scoef[iy] * in[2][IY(iy)];
            }
            // (Dy.in3)_0
            sumy0 = 0;
            for (iy = 0; iy < 6; iy++){
                sumy0 += bnd_coef[iy] * in[2][IY(iy)];
            }
            iy=0;
            out[IY(iy)] = add*out[IY(iy)] - adh * (-f1(par,IY(iy))*sumx/dx -f2(par,IY(iy))*sumz/dz -f3(par,IY(iy))*sumy0/dy + f4(par,IY(iy))*sumy/dy -f5(par,IY(iy))*in[3][IY(iy)]/dt);
        }
    }

    // front bottom left
    #pragma omp parallel for private(sumx,sumy,sumz,sumy0)
    for (int ix=ixmin; ix<nc1; ix++){
        for (int iz=nz-izmax; iz<std::min(nc1,nz-izmin); iz++){

            // (Dx.in1)_0 and (Dz.in2)_0
            int iy=0;
            sumx=0;
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumx += bnd_coef[ix*nc2+j] * in[0][IXYZ(j,iy,nz-1-iz)];
                sumz -= bnd_coef[iz*nc2+j] * in[1][IZ(nz-1-j)];
            }

            // (Sy.in3)_0
            sumy = 0;
            for (iy = 0; iy < 4; iy++){
                sumy += scoef[iy] * in[2][IXYZ(ix,iy,nz-1-iz)];
            }
            // (Dy.in3)_0
            sumy0 = 0;
            for (iy = 0; iy < 6; iy++){
                sumy0 += bnd_coef[iy] * in[2][IXYZ(ix,iy,nz-1-iz)];
            }
            iy=0;
            out[IXYZ(ix,iy,nz-1-iz)] = add*out[IXYZ(ix,iy,nz-1-iz)] - adh * (-f1(par,IXYZ(ix,iy,nz-1-iz))*sumx/dx -f2(par,IXYZ(ix,iy,nz-1-iz))*sumz/dz -f3(par,IXYZ(ix,iy,nz-1-iz))*sumy0/dy + f4(par,IXYZ(ix,iy,nz-1-iz))*sumy/dy -f5(par,IXYZ(ix,iy,nz-1-iz))*in[3][IXYZ(ix,iy,nz-1-iz)]/dt);
        }
    }

    // front middle left
    #pragma omp parallel for private(sumx,sumy,sumz,sumy0)
    for (int ix=ixmin; ix<nc1; ix++){
        for (int iz=izminb; iz<izmaxb; iz++){
        
            // (Dx.in1)_0 and (Dz.in2)_0
            int iy=0;
            sumx=0;
            for (int j=0; j<nc2; j++){
                sumx += bnd_coef[ix*nc2+j] * in[0][IX(j)];
            }
            sumz = coef[0] * (in[1][IZ(iz+1)]-in[1][IZ(iz-1)]) + coef[1] * (in[1][IZ(iz+2)]-in[1][IZ(iz-2)]);

            // (Sy.in3)_0
            sumy = 0;
            for (iy = 0; iy < 4; iy++){
                sumy += scoef[iy] * in[2][IY(iy)];
            }
            // (Dy.in3)_0
            sumy0 = 0;
            for (iy = 0; iy < 6; iy++){
                sumy0 += bnd_coef[iy] * in[2][IY(iy)];
            }
            iy=0;
            out[IY(iy)] = add*out[IY(iy)] - adh * (-f1(par,IY(iy))*sumx/dx -f2(par,IY(iy))*sumz/dz -f3(par,IY(iy))*sumy0/dy + f4(par,IY(iy))*sumy/dy -f5(par,IY(iy))*in[3][IY(iy)]/dt);
        }
    }

    // front top right
    #pragma omp parallel for private(sumx,sumy,sumz,sumy0)
    for (int ix=nx-ixmax; ix<std::min(nc1,nx-ixmin); ix++){
        for (int iz=izmin; iz<nc1; iz++){

            // (Dx.in1)_0 and (Dz.in2)_0
            int iy=0;
            sumx=0;
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumx -= bnd_coef[ix*nc2+j] * in[0][IX(nx-1-j)];
                sumz += bnd_coef[iz*nc2+j] * in[1][IXYZ(nx-1-ix,iy,j)];
            }

            // (Sy.in3)_0
            sumy = 0;
            for (iy = 0; iy < 4; iy++){
                sumy += scoef[iy] * in[2][IXYZ(nx-1-ix,iy,iz)];
            }
            // (Dy.in3)_0
            sumy0 = 0;
            for (iy = 0; iy < 6; iy++){
                sumy0 += bnd_coef[iy] * in[2][IXYZ(nx-1-ix,iy,iz)];
            }
            iy=0;
            out[IXYZ(nx-1-ix,iy,iz)] = add*out[IXYZ(nx-1-ix,iy,iz)] - adh * (-f1(par,IXYZ(nx-1-ix,iy,iz))*sumx/dx -f2(par,IXYZ(nx-1-ix,iy,iz))*sumz/dz -f3(par,IXYZ(nx-1-ix,iy,iz))*sumy0/dy + f4(par,IXYZ(nx-1-ix,iy,iz))*sumy/dy -f5(par,IXYZ(nx-1-ix,iy,iz))*in[3][IXYZ(nx-1-ix,iy,iz)]/dt);
        }
    }

    // front bottom right
    #pragma omp parallel for private(sumx,sumy,sumz,sumy0)
    for (int ix=nx-ixmax; ix<std::min(nc1,nx-ixmin); ix++){
        for (int iz=nz-izmax; iz<std::min(nc1,nz-izmin); iz++){

            // (Dx.in1)_0 and (Dz.in2)_0
            int iy=0;
            sumx=0;
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumx -= bnd_coef[ix*nc2+j] * in[0][IXYZ(nx-1-j,iy,nz-1-iz)];
                sumz -= bnd_coef[iz*nc2+j] * in[1][IXYZ(nx-1-ix,iy,nz-1-j)];
            }

            // (Sy.in3)_0
            sumy = 0;
            for (iy = 0; iy < 4; iy++){
                sumy += scoef[iy] * in[2][IXYZ(nx-1-ix,iy,nz-1-iz)];
            }
            // (Dy.in3)_0
            sumy0 = 0;
            for (iy = 0; iy < 6; iy++){
                sumy0 += bnd_coef[iy] * in[2][IXYZ(nx-1-ix,iy,nz-1-iz)];
            }
            iy=0;
            out[IXYZ(nx-1-ix,iy,nz-1-iz)] = add*out[IXYZ(nx-1-ix,iy,nz-1-iz)] - adh * (-f1(par,IXYZ(nx-1-ix,iy,nz-1-iz))*sumx/dx -f2(par,IXYZ(nx-1-ix,iy,nz-1-iz))*sumz/dz -f3(par,IXYZ(nx-1-ix,iy,nz-1-iz))*sumy0/dy + f4(par,IXYZ(nx-1-ix,iy,nz-1-iz))*sumy/dy -f5(par,IXYZ(nx-1-ix,iy,nz-1-iz))*in[3][IXYZ(nx-1-ix,iy,nz-1-iz)]/dt);
        }
    }

    // front middle right
    #pragma omp parallel for private(sumx,sumy,sumz,sumy0)
    for (int ix=nx-ixmax; ix<std::min(nc1,nx-ixmin); ix++){
        for (int iz=izminb; iz<izmaxb; iz++){

            // (Dx.in1)_0 and (Dz.in2)_0
            int iy=0;
            sumx=0;
            for (int j=0; j<nc2; j++){
                sumx -= bnd_coef[ix*nc2+j] * in[0][IX(nx-1-j)];
            }
            sumz = coef[0] * (in[1][IXYZ(nx-1-ix,iy,iz+1)]-in[1][IXYZ(nx-1-ix,iy,iz-1)]) + coef[1] * (in[1][IXYZ(nx-1-ix,iy,iz+2)]-in[1][IXYZ(nx-1-ix,iy,iz-2)]);

            // (Sy.in3)_0
            sumy = 0;
            for (iy = 0; iy < 4; iy++){
                sumy += scoef[iy] * in[2][IXYZ(nx-1-ix,iy,iz)];
            }
            // (Dy.in3)_0
            sumy0 = 0;
            for (iy = 0; iy < 6; iy++){
                sumy0 += bnd_coef[iy] * in[2][IXYZ(nx-1-ix,iy,iz)];
            }
            iy=0;
            out[IXYZ(nx-1-ix,iy,iz)] = add*out[IXYZ(nx-1-ix,iy,iz)] - adh * (-f1(par,IXYZ(nx-1-ix,iy,iz))*sumx/dx -f2(par,IXYZ(nx-1-ix,iy,iz))*sumz/dz -f3(par,IXYZ(nx-1-ix,iy,iz))*sumy0/dy + f4(par,IXYZ(nx-1-ix,iy,iz))*sumy/dy -f5(par,IXYZ(nx-1-ix,iy,iz))*in[3][IXYZ(nx-1-ix,iy,iz)]/dt);
        }
    }

    // front top middle
    #pragma omp parallel for private(sumx,sumy,sumz,sumy0)
    for (int ix=ixminb; ix<ixmaxb; ix++){
        for (int iz=izmin; iz<nc1; iz++){

            // (Dx.in1)_0 and (Dz.in2)_0
            int iy=0;
            sumx = coef[0] * (in[0][IX(ix+1)]-in[0][IX(ix-1)]) + coef[1] * (in[0][IX(ix+2)]-in[0][IX(ix-2)]);
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumz += bnd_coef[iz*nc2+j] * in[1][IZ(j)];
            }

            // (Sy.in3)_0
            sumy = 0;
            for (iy = 0; iy < 4; iy++){
                sumy += scoef[iy] * in[2][IY(iy)];
            }
            // (Dy.in3)_0
            sumy0 = 0;
            for (iy = 0; iy < 6; iy++){
                sumy0 += bnd_coef[iy] * in[2][IY(iy)];
            }
            iy=0;
            out[IY(iy)] = add*out[IY(iy)] - adh * (-f1(par,IY(iy))*sumx/dx -f2(par,IY(iy))*sumz/dz -f3(par,IY(iy))*sumy0/dy + f4(par,IY(iy))*sumy/dy -f5(par,IY(iy))*in[3][IY(iy)]/dt);
        }
    }

    // front bottom middle
    #pragma omp parallel for private(sumx,sumy,sumz,sumy0)
    for (int ix=ixminb; ix<ixmaxb; ix++){
        for (int iz=nz-izmax; iz<std::min(nc1,nz-izmin); iz++){

            // (Dx.in1)_0 and (Dz.in2)_0
            int iy=0;
            sumx = coef[0] * (in[0][IXYZ(ix+1,iy,nz-1-iz)]-in[0][IXYZ(ix-1,iy,nz-1-iz)]) + coef[1] * (in[0][IXYZ(ix+2,iy,nz-1-iz)]-in[0][IXYZ(ix-2,iy,nz-1-iz)]);
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumz -= bnd_coef[iz*nc2+j] * in[1][IZ(nz-1-j)];
            }

            // (Sy.in3)_0
            sumy = 0;
            for (iy = 0; iy < 4; iy++){
                sumy += scoef[iy] * in[2][IXYZ(ix,iy,nz-1-iz)];
            }
            // (Dy.in3)_0
            sumy0 = 0;
            for (iy = 0; iy < 6; iy++){
                sumy0 += bnd_coef[iy] * in[2][IXYZ(ix,iy,nz-1-iz)];
            }
            iy=0;
            out[IXYZ(ix,iy,nz-1-iz)] = add*out[IXYZ(ix,iy,nz-1-iz)] - adh * (-f1(par,IXYZ(ix,iy,nz-1-iz))*sumx/dx -f2(par,IXYZ(ix,iy,nz-1-iz))*sumz/dz -f3(par,IXYZ(ix,iy,nz-1-iz))*sumy0/dy + f4(par,IXYZ(ix,iy,nz-1-iz))*sumy/dy -f5(par,IXYZ(ix,iy,nz-1-iz))*in[3][IXYZ(ix,iy,nz-1-iz)]/dt);
        }
    }

    // front middle
    #pragma omp parallel for private(sumx,sumy,sumz,sumy0)
    for (int ix=ixminb; ix<ixmaxb; ix++){
        for (int iz=izminb; iz<izmaxb; iz++){

            // (Dx.in1)_0 and (Dz.in2)_0
            int iy=0;
            sumx = coef[0] * (in[0][IX(ix+1)]-in[0][IX(ix-1)]) + coef[1] * (in[0][IX(ix+2)]-in[0][IX(ix-2)]);
            sumz = coef[0] * (in[1][IZ(iz+1)]-in[1][IZ(iz-1)]) + coef[1] * (in[1][IZ(iz+2)]-in[1][IZ(iz-2)]);

            // (Sy.in3)_0
            sumy = 0;
            for (iy = 0; iy < 4; iy++){
                sumy += scoef[iy] * in[2][IY(iy)];
            }
            // (Dy.in3)_0
            sumy0 = 0;
            for (iy = 0; iy < 6; iy++){
                sumy0 += bnd_coef[iy] * in[2][IY(iy)];
            }
            iy=0;
            out[IY(iy)] = add*out[IY(iy)] - adh * (-f1(par,IY(iy))*sumx/dx -f2(par,IY(iy))*sumz/dz -f3(par,IY(iy))*sumy0/dy + f4(par,IY(iy))*sumy/dy -f5(par,IY(iy))*in[3][IY(iy)]/dt);
        }
    }
}

template<expr f1, expr f2, expr f3, expr f4, expr f5>
void esat_neumann_back(bool add, const data_t** in, __restrict data_t* out, int nx, int ny, int nz, data_t dx, data_t dy, data_t dz, data_t dt, int ixmin, int ixmax, int izmin, int izmax, const data_t ** par){

    data_t coef[2] = {2.0/3,-1.0/12};
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3}; 
    data_t h0 = 17.0/48;
    int nc1=4, nc2=6;
    int ixminb=std::max(nc1,ixmin);
    int ixmaxb=std::min(nx - nc1, ixmax);
    int izminb=std::max(nc1,izmin);
    int izmaxb=std::min(ny - nc1, izmax);
    data_t adh = 1.0/(dy*h0);

    data_t sumx, sumy, sumz, sumy0;

    // SAT = - Hy-1.(f1.Dx.in1 + f2.Dz.in2 + f3.Dy.in3 + f4.Sy.in3 - f5.in4/dt)_0
    // Sy is the boundary derivative operator pointing outwards

    // back top left
    #pragma omp parallel for private(sumx,sumy,sumz,sumy0)
    for (int ix=ixmin; ix<nc1; ix++){
        for (int iz=izmin; iz<nc1; iz++){

            // (Dx.in1)_0 and (Dz.in2)_0
            int iy=ny-1;
            sumx=0;
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumx += bnd_coef[ix*nc2+j] * in[0][IX(j)];
                sumz += bnd_coef[iz*nc2+j] * in[1][IZ(j)];
            }

            // (Sy.in3)_0
            sumy = 0;
            for (iy = 0; iy < 4; iy++){
                sumy += scoef[iy] * in[2][IY(ny-1-iy)];
            }
            // (Dy.in3)_0
            sumy0 = 0;
            for (iy = 0; iy < 6; iy++){
                sumy0 -= bnd_coef[iy] * in[2][IY(ny-1-iy)];
            }
            iy=ny-1;
            out[IY(iy)] = add*out[IY(iy)] - adh * (f1(par,IY(iy))*sumx/dx +f2(par,IY(iy))*sumz/dz +f3(par,IY(iy))*sumy0/dy + f4(par,IY(iy))*sumy/dy -f5(par,IY(iy))*in[3][IY(iy)]/dt);
        }
    }

    // back bottom left
    #pragma omp parallel for private(sumx,sumy,sumz,sumy0)
    for (int ix=ixmin; ix<nc1; ix++){
        for (int iz=nz-izmax; iz<std::min(nc1,nz-izmin); iz++){

            // (Dx.in1)_0 and (Dz.in2)_0
            int iy=ny-1;
            sumx=0;
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumx += bnd_coef[ix*nc2+j] * in[0][IXYZ(j,iy,nz-1-iz)];
                sumz -= bnd_coef[iz*nc2+j] * in[1][IZ(nz-1-j)];
            }

            // (Sy.in3)_0
            sumy = 0;
            for (iy = 0; iy < 4; iy++){
                sumy += scoef[iy] * in[2][IXYZ(ix,ny-1-iy,nz-1-iz)];
            }
            // (Dy.in3)_0
            sumy0 = 0;
            for (iy = 0; iy < 6; iy++){
                sumy0 -= bnd_coef[iy] * in[2][IXYZ(ix,ny-1-iy,nz-1-iz)];
            }
            iy=ny-1;
            out[IXYZ(ix,iy,nz-1-iz)] = add*out[IXYZ(ix,iy,nz-1-iz)] - adh * (f1(par,IXYZ(ix,iy,nz-1-iz))*sumx/dx +f2(par,IXYZ(ix,iy,nz-1-iz))*sumz/dz +f3(par,IXYZ(ix,iy,nz-1-iz))*sumy0/dy + f4(par,IXYZ(ix,iy,nz-1-iz))*sumy/dy -f5(par,IXYZ(ix,iy,nz-1-iz))*in[3][IXYZ(ix,iy,nz-1-iz)]/dt);
        }
    }

    // back middle left
    #pragma omp parallel for private(sumx,sumy,sumz,sumy0)
    for (int ix=ixmin; ix<nc1; ix++){
        for (int iz=izminb; iz<izmaxb; iz++){
        
            // (Dx.in1)_0 and (Dz.in2)_0
            int iy=ny-1;
            sumx=0;
            for (int j=0; j<nc2; j++){
                sumx += bnd_coef[ix*nc2+j] * in[0][IX(j)];
            }
            sumz = coef[0] * (in[1][IZ(iz+1)]-in[1][IZ(iz-1)]) + coef[1] * (in[1][IZ(iz+2)]-in[1][IZ(iz-2)]);

            // (Sy.in3)_0
            sumy = 0;
            for (iy = 0; iy < 4; iy++){
                sumy += scoef[iy] * in[2][IY(ny-1-iy)];
            }
            // (Dy.in3)_0
            sumy0 = 0;
            for (iy = 0; iy < 6; iy++){
                sumy0 -= bnd_coef[iy] * in[2][IY(ny-1-iy)];
            }
            iy=ny-1;
            out[IY(iy)] = add*out[IY(iy)] - adh * (f1(par,IY(iy))*sumx/dx +f2(par,IY(iy))*sumz/dz +f3(par,IY(iy))*sumy0/dy + f4(par,IY(iy))*sumy/dy -f5(par,IY(iy))*in[3][IY(iy)]/dt);
        }
    }

    // back top right
    #pragma omp parallel for private(sumx,sumy,sumz,sumy0)
    for (int ix=nx-ixmax; ix<std::min(nc1,nx-ixmin); ix++){
        for (int iz=izmin; iz<nc1; iz++){

            // (Dx.in1)_0 and (Dz.in2)_0
            int iy=ny-1;
            sumx=0;
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumx -= bnd_coef[ix*nc2+j] * in[0][IX(nx-1-j)];
                sumz += bnd_coef[iz*nc2+j] * in[1][IXYZ(nx-1-ix,iy,j)];
            }

            // (Sy.in3)_0
            sumy = 0;
            for (iy = 0; iy < 4; iy++){
                sumy += scoef[iy] * in[2][IXYZ(nx-1-ix,ny-1-iy,iz)];
            }
            // (Dy.in3)_0
            sumy0 = 0;
            for (iy = 0; iy < 6; iy++){
                sumy0 -= bnd_coef[iy] * in[2][IXYZ(nx-1-ix,ny-1-iy,iz)];
            }
            iy=ny-1;
            out[IXYZ(nx-1-ix,iy,iz)] = add*out[IXYZ(nx-1-ix,iy,iz)] - adh * (f1(par,IXYZ(nx-1-ix,iy,iz))*sumx/dx +f2(par,IXYZ(nx-1-ix,iy,iz))*sumz/dz +f3(par,IXYZ(nx-1-ix,iy,iz))*sumy0/dy + f4(par,IXYZ(nx-1-ix,iy,iz))*sumy/dy -f5(par,IXYZ(nx-1-ix,iy,iz))*in[3][IXYZ(nx-1-ix,iy,iz)]/dt);
        }
    }

    // back bottom right
    #pragma omp parallel for private(sumx,sumy,sumz,sumy0)
    for (int ix=nx-ixmax; ix<std::min(nc1,nx-ixmin); ix++){
        for (int iz=nz-izmax; iz<std::min(nc1,nz-izmin); iz++){

            // (Dx.in1)_0 and (Dz.in2)_0
            int iy=ny-1;
            sumx=0;
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumx -= bnd_coef[ix*nc2+j] * in[0][IXYZ(nx-1-j,iy,nz-1-iz)];
                sumz -= bnd_coef[iz*nc2+j] * in[1][IXYZ(nx-1-ix,iy,nz-1-j)];
            }

            // (Sy.in3)_0
            sumy = 0;
            for (iy = 0; iy < 4; iy++){
                sumy += scoef[iy] * in[2][IXYZ(nx-1-ix,ny-1-iy,nz-1-iz)];
            }
            // (Dy.in3)_0
            sumy0 = 0;
            for (iy = 0; iy < 6; iy++){
                sumy0 -= bnd_coef[iy] * in[2][IXYZ(nx-1-ix,ny-1-iy,nz-1-iz)];
            }
            iy=ny-1;
            out[IXYZ(nx-1-ix,iy,nz-1-iz)] = add*out[IXYZ(nx-1-ix,iy,nz-1-iz)] - adh * (f1(par,IXYZ(nx-1-ix,iy,nz-1-iz))*sumx/dx +f2(par,IXYZ(nx-1-ix,iy,nz-1-iz))*sumz/dz +f3(par,IXYZ(nx-1-ix,iy,nz-1-iz))*sumy0/dy + f4(par,IXYZ(nx-1-ix,iy,nz-1-iz))*sumy/dy -f5(par,IXYZ(nx-1-ix,iy,nz-1-iz))*in[3][IXYZ(nx-1-ix,iy,nz-1-iz)]/dt);
        }
    }

    // back middle right
    #pragma omp parallel for private(sumx,sumy,sumz,sumy0)
    for (int ix=nx-ixmax; ix<std::min(nc1,nx-ixmin); ix++){
        for (int iz=izminb; iz<izmaxb; iz++){

            // (Dx.in1)_0 and (Dz.in2)_0
            int iy=ny-1;
            sumx=0;
            for (int j=0; j<nc2; j++){
                sumx -= bnd_coef[ix*nc2+j] * in[0][IX(nx-1-j)];
            }
            sumz = coef[0] * (in[1][IXYZ(nx-1-ix,iy,iz+1)]-in[1][IXYZ(nx-1-ix,iy,iz-1)]) + coef[1] * (in[1][IXYZ(nx-1-ix,iy,iz+2)]-in[1][IXYZ(nx-1-ix,iy,iz-2)]);

            // (Sy.in3)_0
            sumy = 0;
            for (iy = 0; iy < 4; iy++){
                sumy += scoef[iy] * in[2][IXYZ(nx-1-ix,ny-1-iy,iz)];
            }
            // (Dy.in3)_0
            sumy0 = 0;
            for (iy = 0; iy < 6; iy++){
                sumy0 -= bnd_coef[iy] * in[2][IXYZ(nx-1-ix,ny-1-iy,iz)];
            }
            iy=ny-1;
            out[IXYZ(nx-1-ix,iy,iz)] = add*out[IXYZ(nx-1-ix,iy,iz)] - adh * (f1(par,IXYZ(nx-1-ix,iy,iz))*sumx/dx +f2(par,IXYZ(nx-1-ix,iy,iz))*sumz/dz +f3(par,IXYZ(nx-1-ix,iy,iz))*sumy0/dy + f4(par,IXYZ(nx-1-ix,iy,iz))*sumy/dy -f5(par,IXYZ(nx-1-ix,iy,iz))*in[3][IXYZ(nx-1-ix,iy,iz)]/dt);
        }
    }

    // back top middle
    #pragma omp parallel for private(sumx,sumy,sumz,sumy0)
    for (int ix=ixminb; ix<ixmaxb; ix++){
        for (int iz=izmin; iz<nc1; iz++){

            // (Dx.in1)_0 and (Dz.in2)_0
            int iy=ny-1;
            sumx = coef[0] * (in[0][IX(ix+1)]-in[0][IX(ix-1)]) + coef[1] * (in[0][IX(ix+2)]-in[0][IX(ix-2)]);
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumz += bnd_coef[iz*nc2+j] * in[1][IZ(j)];
            }

            // (Sy.in3)_0
            sumy = 0;
            for (iy = 0; iy < 4; iy++){
                sumy += scoef[iy] * in[2][IY(ny-1-iy)];
            }
            // (Dy.in3)_0
            sumy0 = 0;
            for (iy = 0; iy < 6; iy++){
                sumy0 -= bnd_coef[iy] * in[2][IY(ny-1-iy)];
            }
            iy=ny-1;
            out[IY(iy)] = add*out[IY(iy)] - adh * (f1(par,IY(iy))*sumx/dx +f2(par,IY(iy))*sumz/dz +f3(par,IY(iy))*sumy0/dy + f4(par,IY(iy))*sumy/dy -f5(par,IY(iy))*in[3][IY(iy)]/dt);
        }
    }

    // back bottom middle
    #pragma omp parallel for private(sumx,sumy,sumz,sumy0)
    for (int ix=ixminb; ix<ixmaxb; ix++){
        for (int iz=nz-izmax; iz<std::min(nc1,nz-izmin); iz++){

            // (Dx.in1)_0 and (Dz.in2)_0
            int iy=ny-1;
            sumx = coef[0] * (in[0][IXYZ(ix+1,iy,nz-1-iz)]-in[0][IXYZ(ix-1,iy,nz-1-iz)]) + coef[1] * (in[0][IXYZ(ix+2,iy,nz-1-iz)]-in[0][IXYZ(ix-2,iy,nz-1-iz)]);
            sumz=0;
            for (int j=0; j<nc2; j++){
                sumz -= bnd_coef[iz*nc2+j] * in[1][IZ(nz-1-j)];
            }

            // (Sy.in3)_0
            sumy = 0;
            for (iy = 0; iy < 4; iy++){
                sumy += scoef[iy] * in[2][IXYZ(ix,ny-1-iy,nz-1-iz)];
            }
            // (Dy.in3)_0
            sumy0 = 0;
            for (iy = 0; iy < 6; iy++){
                sumy0 -= bnd_coef[iy] * in[2][IXYZ(ix,ny-1-iy,nz-1-iz)];
            }
            iy=ny-1;
            out[IXYZ(ix,iy,nz-1-iz)] = add*out[IXYZ(ix,iy,nz-1-iz)] - adh * (f1(par,IXYZ(ix,iy,nz-1-iz))*sumx/dx +f2(par,IXYZ(ix,iy,nz-1-iz))*sumz/dz +f3(par,IXYZ(ix,iy,nz-1-iz))*sumy0/dy + f4(par,IXYZ(ix,iy,nz-1-iz))*sumy/dy -f5(par,IXYZ(ix,iy,nz-1-iz))*in[3][IXYZ(ix,iy,nz-1-iz)]/dt);
        }
    }

    // back middle
    #pragma omp parallel for private(sumx,sumy,sumz,sumy0)
    for (int ix=ixminb; ix<ixmaxb; ix++){
        for (int iz=izminb; iz<izmaxb; iz++){

            // (Dx.in1)_0 and (Dz.in2)_0
            int iy=ny-1;
            sumx = coef[0] * (in[0][IX(ix+1)]-in[0][IX(ix-1)]) + coef[1] * (in[0][IX(ix+2)]-in[0][IX(ix-2)]);
            sumz = coef[0] * (in[1][IZ(iz+1)]-in[1][IZ(iz-1)]) + coef[1] * (in[1][IZ(iz+2)]-in[1][IZ(iz-2)]);

            // (Sy.in3)_0
            sumy = 0;
            for (iy = 0; iy < 4; iy++){
                sumy += scoef[iy] * in[2][IY(ny-1-iy)];
            }
            // (Dy.in3)_0
            sumy0 = 0;
            for (iy = 0; iy < 6; iy++){
                sumy0 -= bnd_coef[iy] * in[2][IY(ny-1-iy)];
            }
            iy=ny-1;
            out[IY(iy)] = add*out[IY(iy)] - adh * (f1(par,IY(iy))*sumx/dx +f2(par,IY(iy))*sumz/dz +f3(par,IY(iy))*sumy0/dy + f4(par,IY(iy))*sumy/dy -f5(par,IY(iy))*in[3][IY(iy)]/dt);
        }
    }
}

#undef IZ
#undef IX
#undef IY
#undef IXYZ