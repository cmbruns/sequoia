// vector.hxx  mathematical vectors as in multivariable calculus
// Chris Bruns 10-16-93
//

// $Id: vector.cxx,v 1.2 2001/11/29 00:04:43 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/vector.cxx,v 1.2 2001/11/29 00:04:43 bruns Exp $
// $Log: vector.cxx,v $
// Revision 1.2  2001/11/29 00:04:43  bruns
// Removed ^M characters
// Added cvs tags
//

#ifndef _VECTOR_CXX_
#define _VECTOR_CXX_

#include "vector.hxx"
#include "array.cxx"

// arithmetic vector
// no simple Vector multiplication is defined, because this must be
// the domain of Matrices, where Vectors will be treated as n x 1
// in particular, dot product is '^' and cross product is '%'


template<class Type>
Vector<Type> Vector<Type>::operator-() const
{
  const Vector<Type> & v = *this;
  Vector<Type> neg = v;
  uint i;
  for (i = 0; i < dim(); ++i)
    neg[i] = - neg[i];
  return neg; 
}


template<class Type>
Vector<Type> Vector<Type>::operator-(const Vector<Type> & v2) const
{
  const Vector<Type> & v1 = *this;
  uint r = MIN(v1.dim(), v2.dim());
  Vector<Type> diff (r);
  uint i;
  for (i = 0; i < diff.dim(); ++i)
    diff[i] = v1[i] - v2[i];
  return diff; 
}


template<class Type>
Vector<Type> & Vector<Type>::operator-=(const Vector<Type> & v2)
{
  Vector<Type> & v = *this;
  if (dim() != v2.dim())
    error("Vector self-subtraction size mismatch");
  uint i;
  for (i = 0; i < dim(); ++i)
    v[i] -= v2[i];
  return *this; 
}


template<class Type>
Vector<Type> Vector<Type>::operator+(const Vector<Type> & v2) const
{
  const Vector<Type> & v1 = *this;
  if (v1.dim() != v2.dim())
    v1.error("Vector addition size mismatch");
  Vector<Type> sum = v1;
  uint i;
  for (i = 0; i < v1.dim(); ++i)
    sum[i] += v2[i];
  return sum; 
}


template<class Type>
Vector<Type> Vector<Type>::operator*(const Real r) const
{
  const Vector<Type> & v = *this;
  Vector<Type> prod (v.dim());
  uint i;
  for (i = 0; i < v.dim(); ++i) 
    prod[i] = (Type)(v[i] * r);
  return prod; 
}


template<class Type>
Vector<Type> operator*(const Real r, const Vector<Type> & v) 
{
  Vector<Type> prod (v.dim());
  uint i;
  for (i = 0; i < v.dim(); ++i) 
    prod[i] = (Type)(v[i] * r);
  return prod; 
}


template<class Type>
Vector<Type> & Vector<Type>::operator*=(const Real r)
{
  Vector<Type> & v = *this;
  uint i;
  for (i = 0; i < v.dim(); ++i)
    v[i] *= r;
  return v; 
}


template<class Type>
Vector<Type> & Vector<Type>::operator+=(const Vector<Type> & v2)
{
  Vector<Type> & v = *this;
  if (dim() != v2.dim())
    {
      v[300000000] = v2[0];
      error("Vector self-addition size mismatch");
    }
  uint i;
  for (i = 0; i < dim(); ++i)
    v[i] += v2[i];
  return *this; 
}


template<class Type>
Vector<Type> Vector<Type>::operator/(Real r) const
{
  const Vector<Type> & v = *this;
  Vector<Type> prod (dim());
  uint i;
  for (i = 0; i < dim(); ++i) 
    prod[i] = (Type)(v[i] / r);
  return prod; 
}


template<class Type>
Vector<Type> & Vector<Type>::operator/=(Real r)
{
  Vector<Type> & v = *this;
  uint i;
  for (i = 0; i < dim(); ++i)
    v[i] = (Type) (v[i] / r);
  return *this; 
}

template<class Type>
const Vector<Type> Vector<Type>::csub(uint i, uint j) const
{
  Vector<Type> answer(j-i+1);
  uint k;
  for (k=i; k <= j; ++k)
    answer[k-i] = this->operator[](k);
  return answer;
}

template<class Type>
const row_vector<Type> Vector<Type>::T() const
{
  const Array<Type> & answer = *this;
  return answer;
}

template<class Type>
ostream & operator<<(ostream & os, const row_vector<Type> & v)
{
  uint i;
  os << "(";
  for (i=0; i<v.dim(); ++i)
    {
      os.precision(3);
      os.width(7);
      os << v[i];
    }
  os << ")";
  return os;
}

// row_vector must be instantiated first
// template class row_vector<int>;
// template class Vector<int>;

#endif
