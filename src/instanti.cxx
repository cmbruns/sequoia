#ifdef __GNUC__

// $Id: instanti.cxx,v 1.4 2002/06/28 17:21:54 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/instanti.cxx,v 1.4 2002/06/28 17:21:54 bruns Exp $
// $Log: instanti.cxx,v $
// Revision 1.4  2002/06/28 17:21:54  bruns
// Minor changes to avoid compiler warnings with gcc 3.1
//   changed <iostream.h> style includes to modern <iostream> style
//   removed redundant default function parameter values when already specified in
//    the header
//
// Revision 1.3  2001/12/12 18:26:27  bruns
// Many changes to get compilation with gcc 3.0.2
// Removed all templated friend functions (no friends needed)
// Added .c_str() method to Mystring, and removed auto conversion to const char *
//
// Revision 1.2  2001/11/28 23:08:17  bruns
// Added cvs header tags
// Removed ^M characters
//

#include "cmbmacro.cxx"
// template uint MIN(const uint &, const uint &);
template double ABS(const double &);
template int sign(const double &);
template float MAX(const float &,const float &);
template double MAX(const double &,const double &);

// #include "array.cxx"
#include "vector.cxx"

// order is important
// full instantiations should precede partial ones

template class Abstract_array<char>;
template class Abstract_array<Real>;
template class Abstract_array< Array<Real> >; // required

#include "array2d.cxx"
template class Array2d<Real>; // required

// #include "vector.cxx"
template class row_vector<Real>;
template class Vector<Real>;
template class Vector< Vector<Real> >;
template Vector<Real> operator*(Real r, const Vector<Real> & v);
// template Vector<Real> Vector<Real>::operator*(const Real r) const;
// template Vector<Real> & operator*=(Vector<Real> & v, Real r);
template class Abstract_array< Vector<Real> >; // required
template ostream & operator<<(ostream & os, const row_vector<Real> & v);
// template class row_vector< Vector<Real> >;
// template Vector< Vector<Real> > operator*(const Real r, Vector< Vector<Real> >);

#include "matrix.cxx"

// template Vector<Real> operator*(const Matrix<Real> & M, const Vector<Real> & v);
// template Matrix<Real> Matrix<Real>::operator*(const Matrix<Real> & M2) const;
template Matrix<Real> operator*(const Vector<Real> & v1, const row_vector<Real> & v2);
// template class row_vector< MATRIX_TYPE >;

// template class Abstract_array<MATRIX_TYPE>; // required
// template class Abstract_array< Vector<MATRIX_TYPE> >; // required
// template class row_vector< Vector<MATRIX_TYPE> >;
template class Matrix<MATRIX_TYPE>;
// template class Matrix<Real>;
template class subMatrix<Real>;
//template Matrix<Real> operator*(double r, const Matrix<Real> & M);
//template Matrix<Real> operator*(const Vector<Real> & v1, const row_vector<Real> & v2);
template ostream & operator<<(ostream & os, const Matrix<Real> & M);

// template class Matrix<int>;


#include "rvec.hxx"
template class Abstract_array<RVec>;

#include "brookhav.hxx"
template class Abstract_array<pdbatom>;
template class Abstract_array<pdbatom *>;
template class Abstract_array<pdbprotein>;

#include "path_mat.hxx" // path element
template class Abstract_array<path_element>; // required
template class Abstract_array< Array<path_element> >; // required
template class Array2d<path_element>; // required

#include "mystring.hxx"
template class Abstract_array<Mystring>;

#include "sequence.hxx"
template class Abstract_array<SeqRes>;

#include "parse.hxx"
template class Abstract_array<command_token>;

// IRIX5 might need this
#include "variable.hxx"
template class Abstract_array<sequoia_variable>;
#include <iomanip>
// template ostream & operator<<(ostream & os, smanip<int> const & sm);
// template ostream & operator<<(ostream & os, smanip<unsigned long> const & sm);

#endif

#ifdef __TCPLUSPLUS__

#include "array.cxx"
#include "brookhav.hxx"
// this does not help
class Abstract_array < pdbatom > ;

#endif

#ifdef __TCPLUSPLUS__

#include "mystring.hxx"
class Abstract_array<Mystring>;

#endif
