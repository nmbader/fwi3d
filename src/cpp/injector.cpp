#include "injector.hpp"

void delta_m3::findIndicesWeights(const hypercube<data_t> &range, const std::vector<std::vector<data_t> > &loc){

    data_t dx = range.getAxis(2).d;
    data_t dy = range.getAxis(3).d;
    data_t dz = range.getAxis(1).d;
    int nx = range.getAxis(2).n;
    int ny = range.getAxis(3).n;
    int nz = range.getAxis(1).n;
    data_t ox = range.getAxis(2).o;
    data_t oy = range.getAxis(3).o;
    data_t oz = range.getAxis(1).o;

    int (* xind) [3] = (int (*) [3]) _xind.data();
    int (* yind) [3] = (int (*) [3]) _yind.data();
    int (* zind) [3] = (int (*) [3]) _zind.data();
    data_t (* xw) [2][3] = (data_t (*) [2][3]) _xw.data();
    data_t (* yw) [2][3] = (data_t (*) [2][3]) _yw.data();
    data_t (* zw) [2][3] = (data_t (*) [2][3]) _zw.data();

    data_t alpha_x=0, alpha_y=0, alpha_z=0;

    data_t Hcoef[4] = {17.0/48, 59.0/48, 43.0/48, 49.0/48};
    int npts=loc.size();

    for (int itr = 0; itr < npts; itr++){

        // compute the index of the pivot
        xind[itr][1] = std::round((loc[itr][0] - ox)/dx);
        yind[itr][1] = std::round((loc[itr][1] - oy)/dy);
        zind[itr][1] = std::round((loc[itr][2] - oz)/dz);

        // check if the pivot falls at or near the boundary
        // if so, shift it forward or backward one or two samples to make room for the other points xi-1 xi+1
        if (xind[itr][1] == 0) xind[itr][1] = 1;
        if (xind[itr][1] == nx - 1) xind[itr][1] = nx - 2;
        if (yind[itr][1] == 0) yind[itr][1] = 1;
        if (yind[itr][1] == ny - 1) yind[itr][1] = ny - 2;
        if (zind[itr][1] == 0) zind[itr][1] = 1;
        if (zind[itr][1] == nz - 1) zind[itr][1] = nz - 2;

        xind[itr][0] = xind[itr][1] - 1;
        xind[itr][2] = xind[itr][1] + 1;
        yind[itr][0] = yind[itr][1] - 1;
        yind[itr][2] = yind[itr][1] + 1;
        zind[itr][0] = zind[itr][1] - 1;
        zind[itr][2] = zind[itr][1] + 1;

        // compute alpha based on the pivot: x = xi + alpha * h ; xi is the pivot
        alpha_x = (loc[itr][0] - dx * xind[itr][1] - ox) / dx;
        alpha_y = (loc[itr][1] - dy * yind[itr][1] - oy) / dy;
        alpha_z = (loc[itr][2] - dz * zind[itr][1] - oz) / dz;

        // compute the weights based on the formulas given in the header file
        // weights for the injection phase
        xw[itr][0][0] = (pow(alpha_x,2) - alpha_x) / (2 * dx * getHcoef(Hcoef, 4, nx, xind[itr][0]));
        xw[itr][0][1] = ( -pow(alpha_x,2) + 1) / (dx * getHcoef(Hcoef, 4, nx, xind[itr][1]));
        xw[itr][0][2] = (pow(alpha_x,2) + alpha_x) / (2 * dx * getHcoef(Hcoef, 4, nx, xind[itr][2]));

        yw[itr][0][0] = (pow(alpha_y,2) - alpha_y) / (2 * dy * getHcoef(Hcoef, 4, ny, yind[itr][0]));
        yw[itr][0][1] = ( -pow(alpha_y,2) + 1) / (dy * getHcoef(Hcoef, 4, ny, yind[itr][1]));
        yw[itr][0][2] = (pow(alpha_y,2) + alpha_y) / (2 * dy * getHcoef(Hcoef, 4, ny, yind[itr][2]));

        zw[itr][0][0] = (pow(alpha_z,2) - alpha_z) / (2 * dz * getHcoef(Hcoef, 4, nz, zind[itr][0]));
        zw[itr][0][1] = (-pow(alpha_z,2) + 1) / (dz * getHcoef(Hcoef, 4, nz, zind[itr][1]));
        zw[itr][0][2] = (pow(alpha_z,2) + alpha_z) / (2 * dz * getHcoef(Hcoef, 4, nz, zind[itr][2]));

        // weights for the extraction phase
        xw[itr][1][0] = (pow(alpha_x,2) - alpha_x) / 2;
        xw[itr][1][1] = ( -pow(alpha_x,2) + 1);
        xw[itr][1][2] = (pow(alpha_x,2) + alpha_x) / 2;

        yw[itr][1][0] = (pow(alpha_y,2) - alpha_y) / 2;
        yw[itr][1][1] = ( -pow(alpha_y,2) + 1);
        yw[itr][1][2] = (pow(alpha_y,2) + alpha_y) / 2;

        zw[itr][1][0] = (pow(alpha_z,2) - alpha_z) / 2;
        zw[itr][1][1] = ( -pow(alpha_z,2) + 1);
        zw[itr][1][2] = (pow(alpha_z,2) + alpha_z) / 2;
    }
}

void delta_m3::inject(bool add, const data_t ** in, data_t * out, int nx, int ny, int nz, int nt, int npts, int it, int itr_min, int itr_max, const int * xind, const int * yind, const int * zind, const data_t * xw, const data_t * yw, const data_t * zw) const {

    const int (* p_xind) [3] = (const int (*) [3]) xind;
    const int (* p_yind) [3] = (const int (*) [3]) yind;
    const int (* p_zind) [3] = (const int (*) [3]) zind;
    const data_t (* p_xw) [2][3] = (const data_t (*) [2][3]) xw;
    const data_t (* p_yw) [2][3] = (const data_t (*) [2][3]) yw;
    const data_t (* p_zw) [2][3] = (const data_t (*) [2][3]) zw;
    const data_t (* p_in) [nt] = (const data_t (*) [nt]) in[0];
    data_t (* p_out) [nx][nz] = (data_t (*) [nx][nz]) out;

    if (!add) {
        for (int itr = itr_min; itr < itr_max; itr++){
            for (int iy=0; iy<3; iy++){
                for (int ix=0; ix<3; ix++){
                    for (int iz=0; iz<3; iz++){
                        p_out[p_yind[itr][iy]][p_xind[itr][ix]][p_zind[itr][iz]] = 0;
                    }
                }
            } 
        }
    }

    for (int itr = itr_min; itr < itr_max; itr++){
        for (int iy=0; iy<3; iy++){
            for (int ix=0; ix<3; ix++){
                for (int iz=0; iz<3; iz++){
                    p_out[p_yind[itr][iy]][p_xind[itr][ix]][p_zind[itr][iz]] += p_yw[itr][0][iy] * p_xw[itr][0][ix] * p_zw[itr][0][iz] * p_in[itr][it];
                }
            }
        }
    }
}

void delta_m3::extract(bool add, const data_t * in, data_t ** out, int nx, int ny, int nz, int nt, int npts, int it, int itr_min, int itr_max, const int * xind, const int * yind, const int * zind, const data_t * xw, const data_t * yw, const data_t * zw) const {

    const int (* p_xind) [3] = (const int (*) [3]) xind;
    const int (* p_yind) [3] = (const int (*) [3]) yind;
    const int (* p_zind) [3] = (const int (*) [3]) zind;
    const data_t (* p_xw) [2][3] = (const data_t (*) [2][3]) xw;
    const data_t (* p_yw) [2][3] = (const data_t (*) [2][3]) yw;
    const data_t (* p_zw) [2][3] = (const data_t (*) [2][3]) zw;
    const data_t (* p_in) [nx][nz] = (const data_t (*) [nx][nz]) in;
    data_t (* p_out) [nt] = (data_t (*) [nt]) out[0];

    #pragma omp parallel for
    for (int itr = itr_min; itr < itr_max; itr++){
        p_out[itr][it] = add*p_out[itr][it];
        for (int iy=0; iy<3; iy++){
            for (int ix=0; ix<3; ix++){
                for (int iz=0; iz<3; iz++){
                    p_out[itr][it] += p_yw[itr][1][iy] * p_xw[itr][1][ix] * p_zw[itr][1][iz] * p_in[p_yind[itr][iy]][p_xind[itr][ix]][p_zind[itr][iz]];
                }
            }
        }
    }
}

void ddelta_m3::findIndicesWeights(const hypercube<data_t> &range, const std::vector<std::vector<data_t> > &loc){

    data_t dx = range.getAxis(2).d;
    data_t dy = range.getAxis(3).d;
    data_t dz = range.getAxis(1).d;
    int nx = range.getAxis(2).n;
    int ny = range.getAxis(3).n;
    int nz = range.getAxis(1).n;
    data_t ox = range.getAxis(2).o;
    data_t oy = range.getAxis(3).o;
    data_t oz = range.getAxis(1).o;

    int (* xind) [3] = (int (*) [3]) _xind.data();
    int (* yind) [3] = (int (*) [3]) _yind.data();
    int (* zind) [3] = (int (*) [3]) _zind.data();
    data_t (* xw) [2][6] = (data_t (*) [2][6]) _xw.data();
    data_t (* yw) [2][6] = (data_t (*) [2][6]) _yw.data();
    data_t (* zw) [2][6] = (data_t (*) [2][6]) _zw.data();

    data_t alpha_x=0, alpha_y=0, alpha_z=0;

    data_t Hcoef[4] = {17.0/48, 59.0/48, 43.0/48, 49.0/48};
    int npts=loc.size();

    for (int itr = 0; itr < npts; itr++){

        // compute the index of the pivot
        xind[itr][1] = std::round((loc[itr][0] - ox)/dx);
        yind[itr][1] = std::round((loc[itr][1] - oy)/dy);
        zind[itr][1] = std::round((loc[itr][2] - oz)/dz);

        // check if the pivot falls at or near the boundary
        // if so, shift it forward or backward one or two samples to make room for the other points xi-1 xi+1
        if (xind[itr][1] == 0) xind[itr][1] = 1;
        if (xind[itr][1] == nx - 1) xind[itr][1] = nx - 2;
        if (yind[itr][1] == 0) yind[itr][1] = 1;
        if (yind[itr][1] == ny - 1) yind[itr][1] = ny - 2;
        if (zind[itr][1] == 0) zind[itr][1] = 1;
        if (zind[itr][1] == nz - 1) zind[itr][1] = nz - 2;

        xind[itr][0] = xind[itr][1] - 1;
        xind[itr][2] = xind[itr][1] + 1;
        yind[itr][0] = yind[itr][1] - 1;
        yind[itr][2] = yind[itr][1] + 1;
        zind[itr][0] = zind[itr][1] - 1;
        zind[itr][2] = zind[itr][1] + 1;

        // compute alpha based on the pivot: x = xi + alpha * h ; xi is the pivot
        alpha_x = (loc[itr][0] - dx * xind[itr][1] - ox) / dx;
        alpha_y = (loc[itr][1] - dy * yind[itr][1] - oy) / dy;
        alpha_z = (loc[itr][2] - dz * zind[itr][1] - oz) / dz;

        // compute the weights based on the formulas given in the header file
        // weights for the injection phase
        xw[itr][0][0] = (pow(alpha_x,2) - alpha_x) / (2 * dx * getHcoef(Hcoef, 4, nx, xind[itr][0]));
        xw[itr][0][1] = ( -pow(alpha_x,2) + 1) / (dx * getHcoef(Hcoef, 4, nx, xind[itr][1]));
        xw[itr][0][2] = (pow(alpha_x,2) + alpha_x) / (2 * dx * getHcoef(Hcoef, 4, nx, xind[itr][2]));
        xw[itr][0][3] = (1-2*alpha_x) / (2 * dx*dx * getHcoef(Hcoef, 4, nx, xind[itr][0]));
        xw[itr][0][4] = 2*alpha_x / (dx*dx * getHcoef(Hcoef, 4, nx, xind[itr][1]));
        xw[itr][0][5] = (-1-2*alpha_x) / (2 * dx*dx * getHcoef(Hcoef, 4, nx, xind[itr][2]));

        yw[itr][0][0] = (pow(alpha_y,2) - alpha_y) / (2 * dy * getHcoef(Hcoef, 4, ny, yind[itr][0]));
        yw[itr][0][1] = ( -pow(alpha_y,2) + 1) / (dy * getHcoef(Hcoef, 4, ny, yind[itr][1]));
        yw[itr][0][2] = (pow(alpha_y,2) + alpha_y) / (2 * dy * getHcoef(Hcoef, 4, ny, yind[itr][2]));
        yw[itr][0][3] = (1-2*alpha_y) / (2 * dy*dy * getHcoef(Hcoef, 4, ny, yind[itr][0]));
        yw[itr][0][4] = 2*alpha_y / (dy*dy * getHcoef(Hcoef, 4, ny, yind[itr][1]));
        yw[itr][0][5] = (-1-2*alpha_y) / (2 * dy*dy * getHcoef(Hcoef, 4, ny, yind[itr][2]));

        zw[itr][0][0] = (pow(alpha_z,2) - alpha_z) / (2 * dz * getHcoef(Hcoef, 4, nz, zind[itr][0]));
        zw[itr][0][1] = (-pow(alpha_z,2) + 1) / (dz * getHcoef(Hcoef, 4, nz, zind[itr][1]));
        zw[itr][0][2] = (pow(alpha_z,2) + alpha_z) / (2 * dz * getHcoef(Hcoef, 4, nz, zind[itr][2]));
        zw[itr][0][3] = (1-2*alpha_z) / (2 * dz*dz * getHcoef(Hcoef, 4, nz, zind[itr][0]));
        zw[itr][0][4] = 2*alpha_z / (dz*dz * getHcoef(Hcoef, 4, nz, zind[itr][1]));
        zw[itr][0][5] = (-1-2*alpha_z) / (2 * dz*dz * getHcoef(Hcoef, 4, nz, zind[itr][2]));

        // weights for the extraction phase
        xw[itr][1][0] = (pow(alpha_x,2) - alpha_x) / 2;
        xw[itr][1][1] = ( -pow(alpha_x,2) + 1);
        xw[itr][1][2] = (pow(alpha_x,2) + alpha_x) / 2;
        xw[itr][1][3] = (1-2*alpha_x) / (2 * dx);
        xw[itr][1][4] = 2*alpha_x / (dx);
        xw[itr][1][5] = (-1-2*alpha_x) / (2 * dx);

        yw[itr][1][0] = (pow(alpha_y,2) - alpha_y) / 2;
        yw[itr][1][1] = ( -pow(alpha_y,2) + 1);
        yw[itr][1][2] = (pow(alpha_y,2) + alpha_y) / 2;
        yw[itr][1][3] = (1-2*alpha_y) / (2 * dy);
        yw[itr][1][4] = 2*alpha_y / (dy);
        yw[itr][1][5] = (-1-2*alpha_y) / (2 * dy);

        zw[itr][1][0] = (pow(alpha_z,2) - alpha_z) / 2;
        zw[itr][1][1] = ( -pow(alpha_z,2) + 1);
        zw[itr][1][2] = (pow(alpha_z,2) + alpha_z) / 2;
        zw[itr][1][3] = (1-2*alpha_z) / (2 * dz);
        zw[itr][1][4] = 2*alpha_z / (dz);
        zw[itr][1][5] = (-1-2*alpha_z) / (2 * dz);
    }
}

void ddelta_m3::inject(bool add, const data_t ** in, data_t * out, int nx, int ny, int nz, int nt, int npts, int it, int itr_min, int itr_max, const int * xind, const int * yind, const int * zind, const data_t * xw, const data_t * yw, const data_t * zw) const {

    const int (* p_xind) [3] = (const int (*) [3]) xind;
    const int (* p_yind) [3] = (const int (*) [3]) yind;
    const int (* p_zind) [3] = (const int (*) [3]) zind;
    const data_t (* p_xw) [2][6] = (const data_t (*) [2][6]) xw;
    const data_t (* p_yw) [2][6] = (const data_t (*) [2][6]) yw;
    const data_t (* p_zw) [2][6] = (const data_t (*) [2][6]) zw;
    const data_t (* p_inx) [nt] = (const data_t (*) [nt]) in[0];
    const data_t (* p_iny) [nt] = (const data_t (*) [nt]) in[1];
    const data_t (* p_inz) [nt] = (const data_t (*) [nt]) in[2];
    data_t (* p_out) [nx][nz] = (data_t (*) [nx][nz]) out;

    if (!add) {
        for (int itr = itr_min; itr < itr_max; itr++){
            for (int iy=0; iy<3; iy++){
                for (int ix=0; ix<3; ix++){
                    for (int iz=0; iz<3; iz++){
                        p_out[p_yind[itr][iy]][p_xind[itr][ix]][p_zind[itr][iz]] = 0;
                    }
                }
            }
        }
    }

    for (int itr = itr_min; itr < itr_max; itr++){
        for (int iy=0; iy<3; iy++){
            for (int ix=0; ix<3; ix++){
                for (int iz=0; iz<3; iz++){
                    p_out[p_yind[itr][iy]][p_xind[itr][ix]][p_zind[itr][iz]] += (-p_xw[itr][0][ix+3] * p_yw[itr][0][iy] * p_zw[itr][0][iz] * p_inx[itr][it] - p_xw[itr][0][ix] * p_yw[itr][0][iy+3] * p_zw[itr][0][iz] * p_iny[itr][it] - p_xw[itr][0][ix] * p_yw[itr][0][iy] * p_zw[itr][0][iz+3] * p_inz[itr][it]);
                }
            }
        }
    }
}

void ddelta_m3::extract(bool add, const data_t * in, data_t ** out, int nx, int ny, int nz, int nt, int npts, int it, int itr_min, int itr_max, const int * xind, const int * yind, const int * zind, const data_t * xw, const data_t * yw, const data_t * zw) const {

    const int (* p_xind) [3] = (const int (*) [3]) xind;
    const int (* p_yind) [3] = (const int (*) [3]) yind;
    const int (* p_zind) [3] = (const int (*) [3]) zind;
    const data_t (* p_xw) [2][6] = (const data_t (*) [2][6]) xw;
    const data_t (* p_yw) [2][6] = (const data_t (*) [2][6]) yw;
    const data_t (* p_zw) [2][6] = (const data_t (*) [2][6]) zw;
    const data_t (* p_in) [nx][nz] = (const data_t (*) [nx][nz]) in;
    data_t (* p_outx) [nt] = (data_t (*) [nt]) out[0];
    data_t (* p_outy) [nt] = (data_t (*) [nt]) out[1];
    data_t (* p_outz) [nt] = (data_t (*) [nt]) out[2];

    #pragma omp parallel for
    for (int itr = itr_min; itr < itr_max; itr++){
        p_outx[itr][it] = add*p_outx[itr][it];
        p_outy[itr][it] = add*p_outy[itr][it];
        p_outz[itr][it] = add*p_outz[itr][it];
        for (int iy=0; iy<3; iy++){
            for (int ix=0; ix<3; ix++){
                for (int iz=0; iz<3; iz++){
                    p_outx[itr][it] = p_outx[itr][it] - p_xw[itr][1][ix+3] * p_yw[itr][1][iy] * p_zw[itr][1][iz] * p_in[p_yind[itr][iy]][p_xind[itr][ix]][p_zind[itr][iz]];
                    p_outy[itr][it] = p_outy[itr][it] - p_xw[itr][1][ix] * p_yw[itr][1][iy+3] * p_zw[itr][1][iz] * p_in[p_yind[itr][iy]][p_xind[itr][ix]][p_zind[itr][iz]];
                    p_outz[itr][it] = p_outz[itr][it] - p_xw[itr][1][ix] * p_yw[itr][1][iy] * p_zw[itr][1][iz+3] * p_in[p_yind[itr][iy]][p_xind[itr][ix]][p_zind[itr][iz]];
                }
            }
        }
    }
}

void dipole_m3::findIndicesWeights(const hypercube<data_t> &range, const std::vector<std::vector<data_t> > &loc, data_t gl){
    
    data_t dx = range.getAxis(2).d;
    data_t dy = range.getAxis(3).d;
    data_t dz = range.getAxis(1).d;
    int nx = range.getAxis(2).n;
    int ny = range.getAxis(3).n;
    int nz = range.getAxis(1).n;
    data_t ox = range.getAxis(2).o;
    data_t oy = range.getAxis(3).o;
    data_t oz = range.getAxis(1).o;

    int (* xind) [6] = (int (*) [6]) _xind.data();
    int (* yind) [6] = (int (*) [6]) _yind.data();
    int (* zind) [6] = (int (*) [6]) _zind.data();
    data_t (* xw) [2][6] = (data_t (*) [2][6]) _xw.data();
    data_t (* yw) [2][6] = (data_t (*) [2][6]) _yw.data();
    data_t (* zw) [2][6] = (data_t (*) [2][6]) _zw.data();

    data_t alpha_x1=0, alpha_y1=0, alpha_z1=0, alpha_x2=0, alpha_y2=0, alpha_z2=0;

    data_t Hcoef[4] = {17.0/48, 59.0/48, 43.0/48, 49.0/48};
    int npts=loc.size();

    for (int itr = 0; itr < npts; itr++){

        // compute the index of the pivot
        // loc vector contains x,y,z,dip,azimuth in that order
        xind[itr][1] = std::round((loc[itr][0] - gl/2*cos(loc[itr][3])*cos(loc[itr][4]) - ox)/dx);
        xind[itr][4] = std::round((loc[itr][0] + gl/2*cos(loc[itr][3])*cos(loc[itr][4]) - ox)/dx);
        yind[itr][1] = std::round((loc[itr][1] - gl/2*cos(loc[itr][3])*sin(loc[itr][4]) - oy)/dy);
        yind[itr][4] = std::round((loc[itr][1] + gl/2*cos(loc[itr][3])*sin(loc[itr][4]) - oy)/dy);
        zind[itr][1] = std::round((loc[itr][2] - gl/2*sin(loc[itr][3]) - oz)/dz);
        zind[itr][4] = std::round((loc[itr][2] + gl/2*sin(loc[itr][3]) - oz)/dz);

        // check if the pivots fall at or near the boundary
        // if so, shift them forward or backward one or two samples to make room for the other points xi-1 xi+1
        if (xind[itr][1] == 0) xind[itr][1] = 1;
        if (xind[itr][1] == nx - 1) xind[itr][1] = nx - 2;
        if (xind[itr][4] == 0) xind[itr][4] = 1;
        if (xind[itr][4] == nx - 1) xind[itr][4] = nx - 2;
        if (yind[itr][1] == 0) yind[itr][1] = 1;
        if (yind[itr][1] == ny - 1) yind[itr][1] = ny - 2;
        if (yind[itr][4] == 0) yind[itr][4] = 1;
        if (yind[itr][4] == ny - 1) yind[itr][4] = ny - 2;
        if (zind[itr][1] == 0) zind[itr][1] = 1;
        if (zind[itr][1] == nz - 1) zind[itr][1] = nz - 2;
        if (zind[itr][4] == 0) zind[itr][4] = 1;
        if (zind[itr][4] == nz - 1) zind[itr][4] = nz - 2;

        xind[itr][0] = xind[itr][1] - 1;
        xind[itr][2] = xind[itr][1] + 1;
        xind[itr][3] = xind[itr][4] - 1;
        xind[itr][5] = xind[itr][4] + 1;

        yind[itr][0] = yind[itr][1] - 1;
        yind[itr][2] = yind[itr][1] + 1;
        yind[itr][3] = yind[itr][4] - 1;
        yind[itr][5] = yind[itr][4] + 1;

        zind[itr][0] = zind[itr][1] - 1;
        zind[itr][2] = zind[itr][1] + 1;
        zind[itr][3] = zind[itr][4] - 1;
        zind[itr][5] = zind[itr][4] + 1;

        // compute alpha based on the pivot: x = xi + alpha * h ; xi is the pivot
        alpha_x1 = (loc[itr][0] - gl/2*cos(loc[itr][3])*cos(loc[itr][4]) - ox) / dx - xind[itr][1];
        alpha_x2 = (loc[itr][0] + gl/2*cos(loc[itr][3])*cos(loc[itr][4]) - ox) / dx - xind[itr][4];
        alpha_y1 = (loc[itr][1] - gl/2*cos(loc[itr][3])*sin(loc[itr][4]) - oy) / dy - yind[itr][1];
        alpha_y2 = (loc[itr][1] + gl/2*cos(loc[itr][3])*sin(loc[itr][4]) - oy) / dy - yind[itr][4];
        alpha_z1 = (loc[itr][2] - gl/2*sin(loc[itr][3]) - oz) / dz - zind[itr][1];
        alpha_z2 = (loc[itr][2] + gl/2*sin(loc[itr][3]) - oz) / dz - zind[itr][4];

        // compute the weights based on the formulas given in the header file
        // weights for the injection phase
        xw[itr][0][0] = (pow(alpha_x1,2) - alpha_x1) / (2 * dx * getHcoef(Hcoef, 4, nx, xind[itr][0]));
        xw[itr][0][1] = ( -pow(alpha_x1,2) + 1) / (dx * getHcoef(Hcoef, 4, nx, xind[itr][1]));
        xw[itr][0][2] = (pow(alpha_x1,2) + alpha_x1) / (2 * dx * getHcoef(Hcoef, 4, nx, xind[itr][2]));
        xw[itr][0][3] = (pow(alpha_x2,2) - alpha_x2) / (2 * dx * getHcoef(Hcoef, 4, nx, xind[itr][3]));
        xw[itr][0][4] = ( -pow(alpha_x2,2) + 1) / (dx * getHcoef(Hcoef, 4, nx, xind[itr][4]));
        xw[itr][0][5] = (pow(alpha_x2,2) + alpha_x2) / (2 * dx * getHcoef(Hcoef, 4, nx, xind[itr][5]));

        yw[itr][0][0] = (pow(alpha_y1,2) - alpha_y1) / (2 * dy * getHcoef(Hcoef, 4, ny, yind[itr][0]));
        yw[itr][0][1] = ( -pow(alpha_y1,2) + 1) / (dy * getHcoef(Hcoef, 4, ny, yind[itr][1]));
        yw[itr][0][2] = (pow(alpha_y1,2) + alpha_y1) / (2 * dy * getHcoef(Hcoef, 4, ny, yind[itr][2]));
        yw[itr][0][3] = (pow(alpha_y2,2) - alpha_y2) / (2 * dy * getHcoef(Hcoef, 4, ny, yind[itr][3]));
        yw[itr][0][4] = ( -pow(alpha_y2,2) + 1) / (dy * getHcoef(Hcoef, 4, ny, yind[itr][4]));
        yw[itr][0][5] = (pow(alpha_y2,2) + alpha_y2) / (2 * dy * getHcoef(Hcoef, 4, ny, yind[itr][5]));

        zw[itr][0][0] = (pow(alpha_z1,2) - alpha_z1) / (2 * dz * getHcoef(Hcoef, 4, nz, zind[itr][0]));
        zw[itr][0][1] = (-pow(alpha_z1,2) + 1) / (dz * getHcoef(Hcoef, 4, nz, zind[itr][1]));
        zw[itr][0][2] = (pow(alpha_z1,2) + alpha_z1) / (2 * dz * getHcoef(Hcoef, 4, nz, zind[itr][2]));
        zw[itr][0][3] = (pow(alpha_z2,2) - alpha_z2) / (2 * dz * getHcoef(Hcoef, 4, nz, zind[itr][3]));
        zw[itr][0][4] = (-pow(alpha_z2,2) + 1) / (dz * getHcoef(Hcoef, 4, nz, zind[itr][4]));
        zw[itr][0][5] = (pow(alpha_z2,2) + alpha_z2) / (2 * dz * getHcoef(Hcoef, 4, nz, zind[itr][5]));

        // weights for the extraction phase
        xw[itr][1][0] = (pow(alpha_x1,2) - alpha_x1) / 2;
        xw[itr][1][1] = ( -pow(alpha_x1,2) + 1);
        xw[itr][1][2] = (pow(alpha_x1,2) + alpha_x1) / 2;
        xw[itr][1][3] = (pow(alpha_x2,2) - alpha_x2) / 2;
        xw[itr][1][4] = ( -pow(alpha_x2,2) + 1);
        xw[itr][1][5] = (pow(alpha_x2,2) + alpha_x2) / 2;

        yw[itr][1][0] = (pow(alpha_y1,2) - alpha_y1) / 2;
        yw[itr][1][1] = ( -pow(alpha_y1,2) + 1);
        yw[itr][1][2] = (pow(alpha_y1,2) + alpha_y1) / 2;
        yw[itr][1][3] = (pow(alpha_y2,2) - alpha_y2) / 2;
        yw[itr][1][4] = ( -pow(alpha_y2,2) + 1);
        yw[itr][1][5] = (pow(alpha_y2,2) + alpha_y2) / 2;

        zw[itr][1][0] = (pow(alpha_z1,2) - alpha_z1) / 2;
        zw[itr][1][1] = ( -pow(alpha_z1,2) + 1);
        zw[itr][1][2] = (pow(alpha_z1,2) + alpha_z1) / 2;
        zw[itr][1][3] = (pow(alpha_z2,2) - alpha_z2) / 2;
        zw[itr][1][4] = ( -pow(alpha_z2,2) + 1);
        zw[itr][1][5] = (pow(alpha_z2,2) + alpha_z2) / 2;
    }
}

void dipole_m3::inject(bool add, const data_t ** in, data_t * out, int nx, int ny, int nz, int nt, int npts, int it, int itr_min, int itr_max, const int * xind, const int * yind, const int * zind, const data_t * xw, const data_t * yw, const data_t * zw) const {

    const int (* p_xind) [6] = (const int (*) [6]) xind;
    const int (* p_yind) [6] = (const int (*) [6]) yind;
    const int (* p_zind) [6] = (const int (*) [6]) zind;
    const data_t (* p_xw) [2][6] = (const data_t (*) [2][6]) xw;
    const data_t (* p_yw) [2][6] = (const data_t (*) [2][6]) yw;
    const data_t (* p_zw) [2][6] = (const data_t (*) [2][6]) zw;
    const data_t (* p_in) [nt] = (const data_t (*) [nt]) in[0];
    data_t (* p_out) [nx][nz] = (data_t (*) [nx][nz]) out;


    if (!add) {
        for (int itr = itr_min; itr < itr_max; itr++){
            for (int iy=0; iy<3; iy++){
                for (int ix=0; ix<3; ix++){
                    for (int iz=0; iz<3; iz++){
                        p_out[p_yind[itr][iy]][p_xind[itr][ix]][p_zind[itr][iz]] = 0;
                        p_out[p_yind[itr][iy+3]][p_xind[itr][ix+3]][p_zind[itr][iz+3]] = 0;
                    }
                }
            }
        }
    }

    for (int itr = itr_min; itr < itr_max; itr++){
        for (int iy=0; iy<3; iy++){
            for (int ix=0; ix<3; ix++){
                for (int iz=0; iz<3; iz++){
                    p_out[p_yind[itr][iy]][p_xind[itr][ix]][p_zind[itr][iz]] -= p_xw[itr][0][ix] * p_yw[itr][0][iy] * p_zw[itr][0][iz] * p_in[itr][it];
                    p_out[p_yind[itr][iy+3]][p_xind[itr][ix+3]][p_zind[itr][iz+3]] += p_xw[itr][0][ix+3] * p_yw[itr][0][iy+3] * p_zw[itr][0][iz+3] * p_in[itr][it];
                }
            }
        }
    }
}

void dipole_m3::extract(bool add, const data_t * in, data_t ** out, int nx, int ny, int nz, int nt, int npts, int it, int itr_min, int itr_max, const int * xind, const int * yind, const int * zind, const data_t * xw, const data_t * yw, const data_t * zw) const {

    const int (* p_xind) [6] = (const int (*) [6]) xind;
    const int (* p_yind) [6] = (const int (*) [6]) yind;
    const int (* p_zind) [6] = (const int (*) [6]) zind;
    const data_t (* p_xw) [2][6] = (const data_t (*) [2][6]) xw;
    const data_t (* p_yw) [2][6] = (const data_t (*) [2][6]) yw;
    const data_t (* p_zw) [2][6] = (const data_t (*) [2][6]) zw;
    const data_t (* p_in) [nx][nz] = (const data_t (*) [nx][nz]) in;
    data_t (* p_out) [nt] = (data_t (*) [nt]) out[0];

    #pragma omp parallel for
    for (int itr = itr_min; itr < itr_max; itr++){
        p_out[itr][it] = add*p_out[itr][it];
        for (int iy=0; iy<3; iy++){
            for (int ix=0; ix<3; ix++){
                for (int iz=0; iz<3; iz++){
                    p_out[itr][it] += p_xw[itr][1][ix+3] * p_yw[itr][1][iy+3] * p_zw[itr][1][iz+3] * p_in[p_yind[itr][iy+3]][p_xind[itr][ix+3]][p_zind[itr][iz+3]]
                                    - p_xw[itr][1][ix] * p_yw[itr][1][iy] * p_zw[itr][1][iz] * p_in[p_yind[itr][iy]][p_xind[itr][ix]][p_zind[itr][iz]];
                }
            }
        }
    }
}
