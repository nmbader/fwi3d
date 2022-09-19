#include "operator.hpp"

// B-splines functions of order 0, 1, 2, 3
// N0
data_t N0(int i, data_t u, const std::vector<data_t> &uk);

// N1
data_t N1(int i, data_t u, const std::vector<data_t> &uk);

// N2
data_t N2(int i, data_t u, const std::vector<data_t> &uk);

// N3
data_t N3(int i, data_t u, const std::vector<data_t> &uk);

// set the knot vector from control points and multiplicity vectors
void setKnot(std::vector<data_t> &u, std::vector<data_t> &c, std::vector<int> &m);

// populate the coarse spline vector from a dense regular vector using linear interpolation
void fillin3d(std::shared_ptr<vecReg<data_t> > c, const std::shared_ptr<vecReg<data_t> > v, const std::vector<data_t> &cx, const std::vector<data_t> &cy, const std::vector<data_t> &cz);

// same as fillin3d but in 2 dimensions
void fillin2d(std::shared_ptr<vecReg<data_t> > c, const std::shared_ptr<vecReg<data_t> > v, const std::vector<data_t> &cx, const std::vector<data_t> &cz);

// same as fillin3d but in one dimension
void fillin1d(std::shared_ptr<vecReg<data_t> > c, const std::shared_ptr<vecReg<data_t> > v, const std::vector<data_t> &cz);

// duplicate x, y, z slices according to a multiplicity vector for each direction
class duplicate3d : public loper {
protected: 
    std::vector<int> _mx, _my, _mz; // multiplicity of slices
public:
    duplicate3d(){}
    ~duplicate3d(){}
    duplicate3d(const hypercube<data_t> &domain, const std::vector<int> &mx, const std::vector<int> &my, const std::vector<int> &mz){
        successCheck(domain.getNdim()>=3,"Domain for duplicate operator must have at least 3 dimensions\n");
        std::vector<axis<data_t> > axes = domain.getAxes();
        successCheck(axes[0].n == mz.size(),"Domain first axis must match the size of multiplicity vector mz\n");
        successCheck(axes[1].n == mx.size(),"Domain second axis must match the size of multiplicity vector mx\n");
        successCheck(axes[2].n == my.size(),"Domain third axis must match the size of multiplicity vector my\n");

        int nz=0, nx=0, ny=0;
        for (int i=0; i<mz.size(); i++) {successCheck(mz[i]>=1,"Multiplicity must be >=1\n") ; nz += mz[i];}
        for (int i=0; i<mx.size(); i++) {successCheck(mx[i]>=1,"Multiplicity must be >=1\n") ; nx += mx[i];}
        for (int i=0; i<my.size(); i++) {successCheck(my[i]>=1,"Multiplicity must be >=1\n") ; ny += my[i];}

        axes[0].n = nz;
        axes[1].n = nx;
        axes[2].n = ny;

        _domain = domain;
        _range = hypercube<data_t>(axes);
        _mx = mx;
        _my = my;
        _mz = mz;
    }
    duplicate3d * clone() const {
        duplicate3d * op = new duplicate3d(_domain, _mx, _my, _mz);
        return op;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        return checkCompatible(mod, dat);
    }
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        successCheck(domain.isCompatible(_domain) && range.isCompatible(_range),"Domain or range and incompatible with the operator\n");
        _domain = domain;
        _range = range;
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
};

// same as duplicate but in 2 dimensions
class duplicate2d : public loper {
protected: 
    std::vector<int> _mx, _mz; // multiplicity of slices
public:
    duplicate2d(){}
    ~duplicate2d(){}
    duplicate2d(const hypercube<data_t> &domain, const std::vector<int> &mx, const std::vector<int> &mz){
        successCheck(domain.getNdim()>=2,"Domain for duplicate operator must have at least 2 dimensions\n");
        std::vector<axis<data_t> > axes = domain.getAxes();
        successCheck(axes[0].n == mz.size(),"Domain first axis must match the size of multiplicity vector mz\n");
        successCheck(axes[1].n == mx.size(),"Domain second axis must match the size of multiplicity vector mx\n");

        int nz=0, nx=0;
        for (int i=0; i<mz.size(); i++) {successCheck(mz[i]>=1,"Multiplicity must be >=1\n") ; nz += mz[i];}
        for (int i=0; i<mx.size(); i++) {successCheck(mx[i]>=1,"Multiplicity must be >=1\n") ; nx += mx[i];}

        axes[0].n = nz;
        axes[1].n = nx;

        _domain = domain;
        _range = hypercube<data_t>(axes);
        _mx = mx;
        _mz = mz;
    }
    duplicate2d * clone() const {
        duplicate2d * op = new duplicate2d(_domain, _mx,_mz);
        return op;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        return checkCompatible(mod, dat);
    }
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        successCheck(domain.isCompatible(_domain) && range.isCompatible(_range),"Domain or range and incompatible with the operator\n");
        _domain = domain;
        _range = range;
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
};

// same as duplicate but in one dimension
class duplicate1d : public loper {
protected: 
    std::vector<int> _mz; // multiplicity vector
public:
    duplicate1d(){}
    ~duplicate1d(){}
    duplicate1d(const hypercube<data_t> &domain, const std::vector<int> &mz){
        std::vector<axis<data_t> > axes = domain.getAxes();
        successCheck(axes[0].n == mz.size(),"Domain first axis must match the size of multiplicity vector mz\n");
        int nz=0;
        for (int i=0; i<mz.size(); i++) {successCheck(mz[i]>=1,"Multiplicity must be >=1\n") ; nz += mz[i];}
        axes[0].n = nz;
        _domain = domain;
        _range = hypercube<data_t>(axes);
        _mz = mz;
    }
    duplicate1d * clone() const {
        duplicate1d * op = new duplicate1d(_domain, _mz);
        return op;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        return checkCompatible(mod, dat);
    }
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        successCheck(domain.isCompatible(_domain) && range.isCompatible(_range),"Domain or range and incompatible with the operator\n");
        _domain = domain;
        _range = range;
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
};

// Cubic B-splines model preconditioner
// The 1st and last knots are always repeated 4 times
// The control points are assumed to be the same as the knots except the 1st and last
// control points that are repeated only twice
class bsplines3d : public loper {
protected:
    std::vector<data_t> _kx, _ky, _kz; // knot vectors
    std::vector<int> _kxmin, _kymin, _kzmin; // first useful index of the knot vectors
public:
    bsplines3d(){}
    ~bsplines3d(){}
    bsplines3d(const hypercube<data_t> &domain, const hypercube<data_t> &range, const std::vector<data_t> &kx, const std::vector<data_t> &ky, const std::vector<data_t> &kz){
        successCheck(domain.getNdim()>=3,"Domain must have at least 3 dimensions\n");
        successCheck(domain.getAxis(1).n == kz.size()-4,"Domain first axis must have 4 samples less than the z knot vector\n");
        successCheck(domain.getAxis(2).n == kx.size()-4,"Domain second axis must have 4 samples less than the x knot vector\n");
        successCheck(domain.getAxis(3).n == ky.size()-4,"Domain third axis must have 4 samples less than the y knot vector\n");
        successCheck(domain.getNdim() == range.getNdim(),"Domain and range must have the same number of dimensions\n");
        for (int i=4; i<=domain.getNdim(); i++) successCheck(domain.getAxis(i).n == range.getAxis(i).n,"Domain and range must have the same size in the other dimensions\n");
        for (int i=0; i<kx.size()-1; i++){
            successCheck(kx[i]<=kx[i+1],"Ux knot vector entries must be in ascending order.\n");
        }
        for (int i=0; i<ky.size()-1; i++){
            successCheck(ky[i]<=ky[i+1],"Uy knot vector entries must be in ascending order.\n");
        }
        for (int i=0; i<kz.size()-1; i++){
            successCheck(kz[i]<=kz[i+1],"Uz knot vector entries must be in ascending order.\n");
        }
        for (int i=1; i<4; i++){
            successCheck((ky[i] == ky[0]) && (kx[i] == kx[0]) && (kz[i]==kz[0]) && (ky[ky.size()-1-i]==ky[ky.size()-1]) && (kx[kx.size()-1-i]==kx[kx.size()-1]) && (kz[kz.size()-1-i]==kz[kz.size()-1]),"The first and last 4 knots in the knot vectors must be the same.\n");
        }
        successCheck(kz[0] == range.getAxis(1).o,"The origin of z knot and range z axis must be the same\n");
        successCheck(kx[0] == range.getAxis(2).o,"The origin of x knot and range x axis must be the same\n");
        successCheck(ky[0] == range.getAxis(3).o,"The origin of y knot and range y axis must be the same\n");

        _domain = domain;
        _range = range;
        _kx = kx;
        _ky = ky;
        _kz = kz;

        _kx[kx.size()-1] += 1e-06;
        _kx[kx.size()-2] += 1e-06;
        _kx[kx.size()-3] += 1e-06;
        _kx[kx.size()-4] += 1e-06;
        _ky[ky.size()-1] += 1e-06;
        _ky[ky.size()-2] += 1e-06;
        _ky[ky.size()-3] += 1e-06;
        _ky[ky.size()-4] += 1e-06;
        _kz[kz.size()-1] += 1e-06;
        _kz[kz.size()-2] += 1e-06;
        _kz[kz.size()-3] += 1e-06;
        _kz[kz.size()-4] += 1e-06;

        axis<data_t> Z = _range.getAxis(1);
        axis<data_t> X = _range.getAxis(2);
        axis<data_t> Y = _range.getAxis(3);
        data_t x,y,z;
        _kzmin.resize(Z.n,0);
        _kxmin.resize(X.n,0);
        _kymin.resize(Y.n,0);
        int count=0;
        for (int i=1; i<Z.n; i++){
            z = Z.o + i*Z.d;
            while (z > _kz[count]) count++;
            _kzmin[i] = count-4;
        }
        count = 0;
        for (int i=1; i<X.n; i++){
            x = X.o + i*X.d;
            while (x > _kx[count]) count++;
            _kxmin[i] = count-4;
        }
        count = 0;
        for (int i=1; i<Y.n; i++){
            y = Y.o + i*Y.d;
            while (y > _ky[count]) count++;
            _kymin[i] = count-4;
        }
    }
    bsplines3d * clone() const {
        bsplines3d * op = new bsplines3d(_domain, _range, _kx,_ky,_kz);
        return op;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        return checkSame(mod, dat);
    }
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        successCheck(domain==_domain && range==_range,"Domain or range and different than the operator\n");
        _domain = domain;
        _range = range;
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
};

// same as bsplines3d but in 2 dimensions
class bsplines2d : public loper {
protected:
    std::vector<data_t> _kx, _kz; // knot vectors
    std::vector<int> _kxmin, _kzmin; // first useful index of the knot vectors
public:
    bsplines2d(){}
    ~bsplines2d(){}
    bsplines2d(const hypercube<data_t> &domain, const hypercube<data_t> &range, const std::vector<data_t> &kx, const std::vector<data_t> &kz){
        successCheck(domain.getNdim()>=2,"Domain must have at least 2 dimensions\n");
        successCheck(domain.getAxis(1).n == kz.size()-4,"Domain first axis must have 4 samples less than the z knot vector\n");
        successCheck(domain.getAxis(2).n == kx.size()-4,"Domain second axis must have 4 samples less than the x knot vector\n");
        successCheck(domain.getNdim() == range.getNdim(),"Domain and range must have the same number of dimensions\n");
        for (int i=3; i<=domain.getNdim(); i++) successCheck(domain.getAxis(i).n == range.getAxis(i).n,"Domain and range must have the same size in the other dimensions\n");
        for (int i=0; i<kx.size()-1; i++){
            successCheck(kx[i]<=kx[i+1],"Ux knot vector entries must be in ascending order.\n");
        }
        for (int i=0; i<kz.size()-1; i++){
            successCheck(kz[i]<=kz[i+1],"Uz knot vector entries must be in ascending order.\n");
        }
        for (int i=1; i<4; i++){
            successCheck((kx[i] == kx[0]) && (kz[i]==kz[0]) && (kx[kx.size()-1-i]==kx[kx.size()-1]) && (kz[kz.size()-1-i]==kz[kz.size()-1]),"The first and last 4 knots in the knot vectors must be the same.\n");
        }
        successCheck(kz[0] == range.getAxis(1).o,"The origin of z knot and range z axis must be the same\n");
        successCheck(kx[0] == range.getAxis(2).o,"The origin of x knot and range x axis must be the same\n");

        _domain = domain;
        _range = range;
        _kx = kx;
        _kz = kz;

        _kx[kx.size()-1] += 1e-06;
        _kx[kx.size()-2] += 1e-06;
        _kx[kx.size()-3] += 1e-06;
        _kx[kx.size()-4] += 1e-06;
        _kz[kz.size()-1] += 1e-06;
        _kz[kz.size()-2] += 1e-06;
        _kz[kz.size()-3] += 1e-06;
        _kz[kz.size()-4] += 1e-06;

        axis<data_t> Z = _range.getAxis(1);
        axis<data_t> X = _range.getAxis(2);
        data_t x,z;
        _kzmin.resize(Z.n,0);
        _kxmin.resize(X.n,0);
        int count=0;
        for (int i=1; i<Z.n; i++){
            z = Z.o + i*Z.d;
            while (z > _kz[count]) count++;
            _kzmin[i] = count-4;
        }
        count = 0;
        for (int i=1; i<X.n; i++){
            x = X.o + i*X.d;
            while (x > _kx[count]) count++;
            _kxmin[i] = count-4;
        }
    }
    bsplines2d * clone() const {
        bsplines2d * op = new bsplines2d(_domain, _range, _kx,_kz);
        return op;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        return checkSame(mod, dat);
    }
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        successCheck(domain==_domain && range==_range,"Domain or range and different than the operator\n");
        _domain = domain;
        _range = range;
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
};

// same as bsplines3d but in one dimension 
class bsplines1d : public loper {
protected:
    std::vector<data_t> _kz; // knot vector
    std::vector<int> _kzmin; // first useful index of the knot vector
public:
    bsplines1d(){}
    ~bsplines1d(){}
    bsplines1d(const hypercube<data_t> &domain, const hypercube<data_t> &range, const std::vector<data_t> &kz){
        successCheck(domain.getAxis(1).n == kz.size()-4,"Domain first axis must have 4 samples less than the z knot vector\n");
        successCheck(domain.getNdim() == range.getNdim(),"Domain and range must have the same number of dimensions\n");
        for (int i=2; i<=domain.getNdim(); i++) successCheck(domain.getAxis(i).n == range.getAxis(i).n,"Domain and range must have the same size in the other dimensions\n");
        for (int i=0; i<kz.size()-1; i++){
            successCheck(kz[i]<=kz[i+1],"Uz knot vector entries must be in ascending order.\n");
        }
        for (int i=1; i<4; i++){
            successCheck( (kz[i]==kz[0]) && (kz[kz.size()-1-i]==kz[kz.size()-1]),"The first and last 4 knots in the knot vector must be the same.\n");
        }
        successCheck(kz[0] == range.getAxis(1).o,"The origin of z knot and range z axis must be the same\n");

        _domain = domain;
        _range = range;
        _kz = kz;

        _kz[kz.size()-1] += 1e-06;
        _kz[kz.size()-2] += 1e-06;
        _kz[kz.size()-3] += 1e-06;
        _kz[kz.size()-4] += 1e-06;

        axis<data_t> Z = _range.getAxis(1);
        data_t z;
        _kzmin.resize(Z.n,0);
        int count=0;
        for (int i=1; i<Z.n; i++){
            z = Z.o + i*Z.d;
            while (z > _kz[count]) count++;
            _kzmin[i] = count-4;
        }
    }
    bsplines1d * clone() const {
        bsplines1d * op = new bsplines1d(_domain, _range, _kz);
        return op;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        return checkSame(mod, dat);
    }
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        successCheck(domain==_domain && range==_range,"Domain or range and different than the operator\n");
        _domain = domain;
        _range = range;
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
};

// operator doing the opposite (not the inverse) of bsplines3d; it populates the B-spline model using a gridded model with bilinear interpolation
class bsfillin3d : public loper {
protected:
    std::vector<data_t> _controlx; // control points
    std::vector<data_t> _controly; // control points
    std::vector<data_t> _controlz; // control points
public:
    bsfillin3d(){}
    ~bsfillin3d(){}
    bsfillin3d(const hypercube<data_t> &domain, const std::vector<data_t> &controlx, const std::vector<data_t> &controly, const std::vector<data_t> &controlz){

        _domain = domain;
        std::vector<axis<data_t> > axes = domain.getAxes();
        axes[0]=controlz.size();
        axes[1]=controlx.size();
        axes[2]=controly.size();
        _range = hypercube<data_t>(axes);
        _controlx = controlx;
        _controly = controly;
        _controlz = controlz;
    }
    bsfillin3d * clone() const {
        bsfillin3d * op = new bsfillin3d(_domain, _controlx, _controly, _controlz);
        return op;
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
};
