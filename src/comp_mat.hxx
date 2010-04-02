#ifndef COMP_MAT_HXX
#define COMP_MAT_HXX

#ifdef _MSC_VER
#pragma warning(disable : 4996)
#endif

// $Id: comp_mat.hxx,v 1.2 2001/11/28 23:08:17 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/comp_mat.hxx,v 1.2 2001/11/28 23:08:17 bruns Exp $
// $Log: comp_mat.hxx,v $
// Revision 1.2  2001/11/28 23:08:17  bruns
// Added cvs header tags
// Removed ^M characters
//

#ifdef IRIX
#include "stdlib.h"
#else
#include <stdlib.h>
#endif

#include "mystring.hxx"
// #include <String.h>
#include "matrix.hxx"

#define BLOSUM62 862

class Comparison_matrix : public Matrix<MATRIX_TYPE>
{
friend ostream & operator<<(ostream & os, const Comparison_matrix & s);
friend istream & operator>>(istream & is, Comparison_matrix & s);
private:
  Mystring title;
//  Matrix<MATRIX_TYPE> element;
  Mystring used; // contains list of residues used
public:
  Comparison_matrix(void);
  Comparison_matrix(int);
  ~Comparison_matrix(void);
//  const Matrix_Row<MATRIX_TYPE> & operator[](unsigned int m) const;
//  Matrix_Row<MATRIX_TYPE> & operator[](unsigned int m);
};

#endif
