#ifndef _MATRIX_HXX_
#define _MATRIX_HXX_

// $Id: matrix.hxx,v 1.5 2002/06/28 17:21:54 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/matrix.hxx,v 1.5 2002/06/28 17:21:54 bruns Exp $
// $Log: matrix.hxx,v $
// Revision 1.5  2002/06/28 17:21:54  bruns
// Minor changes to avoid compiler warnings with gcc 3.1
//   changed <iostream.h> style includes to modern <iostream> style
//   removed redundant default function parameter values when already specified in
//    the header
//
// Revision 1.4  2001/12/12 18:26:27  bruns
// Many changes to get compilation with gcc 3.0.2
// Removed all templated friend functions (no friends needed)
// Added .c_str() method to Mystring, and removed auto conversion to const char *
//
// Revision 1.3  2001/12/12 03:52:06  cdputnam
// Comment out friend functions: operator<<, operator>>, Matrix<Type> operator*
// Reimplement only operator<< using public fn Matrix<Type>::print(ostream&)
// Now compiles with gcc 2.95 and 2.96.
//
// Revision 1.2  2001/11/28 23:08:17  bruns
// Added cvs header tags
// Removed ^M characters
//

#include <iostream>
// #include <math.h>
#include "vector.hxx"

#define MATRIX_TYPE double

template <class Type>
class subMatrix;

template <class Type>
class Matrix : public Vector< Vector<Type> >
{
  // <> is needed to instantiate template in absence of -fguiding-decls
// friend ostream & operator<< (ostream & os, const Matrix<Type> & M);
// friend istream & operator>> (istream & is, Matrix<Type> & M);
// friend Matrix<Type> operator*(double r, const Matrix<Type> & M);
// friend Matrix<Type> operator*(const Vector<Type> & v1, const row_vector<Type> & v2);
protected:
  void init_mat(uint i, uint j, Type val);
  void init_mat(uint i, uint j);
  Matrix(int r) : Vector< Vector<Type> >(r) {}
public:
  Matrix(void) : Vector< Vector<Type> >() {}
  Matrix(uint r, uint c, Type val) {init_mat(r,c,val);}
  Matrix(uint r, uint c) {init_mat(r,c);}
  Matrix(const Matrix & M) : Vector< Vector<Type> > (M) {}
  Matrix(const Vector< Vector<Type> > & M) : Vector< Vector<Type> > (M) {}
  Matrix(int r, int c, const Type* M); // initialize from traditional 2D array
  ~Matrix(void) {}

  void print( ostream& os ) const;
  
  uint m() const {return this->dim();}
  uint n() const 
    {
      if (!m()) return 0;
      else return this->operator[](0).dim();
    }

  Vector<Type> operator*(const Vector<Type> & v) const;
  Matrix<Type> operator*(const Matrix<Type> & M2) const;

  const Matrix<Type> csub(uint i1, uint i2, uint j1, uint j2) const;
  const Matrix<Type> sub(uint i1, uint i2, uint j1, uint j2) const
    {return csub(i1,i2,j1,j2);}
  subMatrix<Type> sub(uint i1, uint i2, uint j1, uint j2);
  const Vector<Type> row(uint i) const {return this->operator[](i);}
  Vector<Type> row(uint i) {return this->operator[](i);}
  const Vector<Type> ccol(uint i) const;
  const Vector<Type> col(uint i) const {return ccol(i);}
//  matrix_column<Type> col(uint i);

  Matrix<Type> T() const; // transpose
  Type mean() const;
  Type std_dev() const;
//  Type sum();
  const Vector<Type> diag() const;
};

template<class Type>
ostream& operator<< (ostream& os, const Matrix<Type> & M );

template<class Type>
class subMatrix
{
friend class Matrix<Type>;
private:
  Matrix<Type> * parent;
  uint i_dim, j_dim;
  uint i_offset, j_offset;
  subMatrix() {parent = NULL; i_dim = j_dim = i_offset = j_offset = 0;}
public:
  subMatrix & operator=(const Matrix<Type> & M);
  ~subMatrix() {}
};


template<class Type>
class matrix_column
{
friend class Matrix<Type>;
private:
  Matrix<Type> * parent;
  uint column;
  matrix_column() {parent = NULL; column = 0;}
public:
  matrix_column<Type> & operator=(const Vector<Type> & v);
  ~matrix_column() {}
  const Vector<Type> sub(uint i, uint j) const;
};

template<class Type>
ostream & operator<< (ostream & os, const Matrix<Type> & M);
template<class Type>
istream & operator>> (istream & is, Matrix<Type> & M);

#ifdef __TCPLUSPLUS__
 #include "matrix.cxx"
#endif

#endif

