#ifndef __RMAT_HXX__
#define __RMAT_HXX__

// $Id: rmat.hxx,v 1.3 2001/12/12 18:26:27 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/rmat.hxx,v 1.3 2001/12/12 18:26:27 bruns Exp $
// $Log: rmat.hxx,v $
// Revision 1.3  2001/12/12 18:26:27  bruns
// Many changes to get compilation with gcc 3.0.2
// Removed all templated friend functions (no friends needed)
// Added .c_str() method to Mystring, and removed auto conversion to const char *
//
// Revision 1.2  2001/11/28 23:40:22  bruns
// Added cvs header tags
// Removed ^M characters
//

// #include <math.h>
#ifndef MAC_OS_X
#include <values.h> // maybe SUN only? for unit roundoff calculation
#else
#include <limits.h>
#include <float.h>
#endif

// djgpp has no FSIGNIF in values.h, so try this
#ifndef FSIGNIF
#ifdef FLT_MANT_DIG
#define FSIGNIF FLT_MANT_DIG
#endif
#endif

#ifndef FSIGNIF
#ifdef FLOATBITS
#define FSIGNIF     (FLOATBITS  - _FEXPLEN + _HIDDENBIT - 1)
#endif
#endif

#ifndef FSIGNIF
#define FSIGNIF     (BITS(float)  - _FEXPLEN + _HIDDENBIT - 1)
#endif

#include "matrix.hxx"
#include "rvec.hxx"

class RMat : public Matrix<Real>
{
// friend RMat operator*(double r, const RMat & M);
protected:
public:
  RMat(void) {}
  RMat(int m) : Matrix<Real> (m) {}
  RMat(uint m, uint n, Real val = 0) : Matrix<Real> (m, n, val) {}
  RMat(const RMat & m) : Matrix<Real> (m) {}
  RMat(const Matrix<Real> & m) : Matrix<Real> (m) {}
  RMat(const Vector< Vector<Real> > & m) : Matrix<Real> (m) {}
  ~RMat(void) {}
  Real operator^(const RMat & m2) const; // dot product

  RMat sym_eigen() const; // eigenvectors and eigenvalues
  void sym_tridiag(RMat & Q);
  void house_hess(RMat & Q);
  void sym_QRstep(uint i, uint j, RMat & Q);
//  void sym_QRstep(uint i, uint j);
  void row_rot(const uint i, const uint k, const Vector2D cs);
  void col_rot(const uint i, const uint k, const Vector2D cs);
  bool is_orthogonal() const;
  bool is_normal() const;
  bool is_orthonormal() const;
  RMat strip_trans() const; // remove translational component of a homogeneous matrix
  double angle(const RMat & m2) const; // angle between two matrices
  double angle2(const RMat & m2) const; // angle between two matrices
  double determinant() const;
  double det() const {return this->determinant();}
};

// What is wrong with these?
// How else do you get matrices by multiplying vectors?
//Real operator*(const row_vector<Real> & v1, const Vector<Real> & v2);
//RMat operator*(const RVec & v1, const row_vector<Real> & v2);
//RVec operator*(const RMat & M, const row_vector<Real> & v);
RMat eye(uint m, uint n);
RMat eye(uint m);
Vector2D givens(Real a, Real b);
RMat trans_mat(const RVec & v);
RMat row_house(const RMat & A, const RVec & v);
RMat col_house(const RMat & A, const RVec & v);

Vector3D transform(const RMat & M, const Vector3D v);
Vector3D vec3_transform(const RMat & M, const Vector3D v); 

#endif

