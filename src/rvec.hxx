#ifndef __RVEC_HXX__
#define __RVEC_HXX__

#include "vector.hxx"

class RVec : public Vector<Real>
{
protected:
public:
  RVec() {}
  // linux no like  RVec(uint n) : Vector<Real>(n) {}
  RVec(uint n, Real fill=0) : Vector<Real>(n, fill) {}
  RVec(const RVec & r) : Vector<Real> (r) {}
  //  RVec(const Vector<Real> & r) : Vector<Real> (r) {}
  //  RVec(Vector<Real> r) : Vector<Real> (r) {} // Need type conversion for many functions
  RVec(const Array<Real> & r) : Vector<Real> (r) {}
  ~RVec() {}
  // Real length() const {return sqrt(*this ^ *this);}
  RVec unit() const;

  Real pnorm(double p = 2) const;  // p-norm -- Golub and Van Loan p.53
  RVec house() const;

  // Real operator^(const RVec & v2) const;
};

ostream & operator<<(ostream & os, const RVec & v);

class Vector3D : public RVec
{
//friend Vector3D operator*(Real r, const Vector3D & v); // scaling
protected:
public:
  Vector3D() : RVec(3, 0.0) {}
  Vector3D(Real x, Real y, Real z);
  //  Vector3D(const Vector<Real> v);
  Vector3D(Vector<Real> v); // Need type conversion for many functions
  ~Vector3D() {}
  Vector3D(const Vector3D & a2) 
    {*this = a2;}
  Vector3D & operator=(const Vector3D & a2);
  Real & x() {return this->operator[](0);}
  Real & y() {return this->operator[](1);}
  Real & z() {return this->operator[](2);}
  const Real & x() const {return this->operator[](0);}
  const Real & y() const {return this->operator[](1);}
  const Real & z() const {return this->operator[](2);}
  Vector3D & operator%=(const Vector3D & v); // cross product
  Real distance(const Vector3D & v2) const;
};

// Making these non-member functions seems to release a lot of type
// conversion hassles which can not be solved in any other way

Vector3D operator%(const Vector3D & v1, const Vector3D & v2); // cross product
Real operator^(const Vector<Real> & v, const Vector<Real> & v2);
Real length(const Vector<Real> & v);
// RVec operator*(Real r, const RVec & v);

class Vector2D : public RVec
{
protected:
public:
  Vector2D() : RVec (2, 0.0) {}
  Vector2D(Real x, Real y); // constructor
  ~Vector2D() {}
  Real & x() {return this->operator[](0);}
  Real & y() {return this->operator[](1);}
  const Real & x() const {return this->operator[](0);}
  const Real & y() const {return this->operator[](1);}
};


Vector2D Vec(Real x, Real y);

Vector3D Vec(Real x, Real y, Real z);

#endif
