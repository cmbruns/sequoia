// vector.hxx  mathematical vectors as in multivariable calculus
// Chris Bruns 10-16-93
//
#ifndef __VECTOR_HXX__
#define __VECTOR_HXX__

// $Id: vector.hxx,v 1.4 2001/12/12 18:26:27 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/vector.hxx,v 1.4 2001/12/12 18:26:27 bruns Exp $
// $Log: vector.hxx,v $
// Revision 1.4  2001/12/12 18:26:27  bruns
// Many changes to get compilation with gcc 3.0.2
// Removed all templated friend functions (no friends needed)
// Added .c_str() method to Mystring, and removed auto conversion to const char *
//
// Revision 1.3  2001/12/12 03:52:49  cdputnam
// Remove friend function: Vector<Type> operator*(const Real,const Vector<Type>&)
// Now compiles with gcc 2.95 and 2.96
//
// Revision 1.2  2001/11/29 00:04:43  bruns
// Removed ^M characters
// Added cvs tags
//

#include <math.h>
#include <iostream>

#include "array.hxx"

using namespace std;

template<class Type>
class row_vector;

// arithmetic vector
// no simple Vector multiplication is defined, because this must be
// the domain of Matrices, where Vectors will be treated as n x 1
// in particular, dot product is '^' and cross product is '%'
template<class Type>
class Vector : public Array<Type>
{
friend class row_vector<Type>;
  // <> is needed to instantiate template in absence of -fguiding-decls
// friend Vector<Type> operator*(const Real r, const Vector<Type> & v);
protected:
  Vector(const Array<Type> & v) : Array<Type>(v) {}
public:
  Vector() {}
  Vector(uint n) : Array<Type>(n) {}
  Vector(uint n, const Type & fill) : Array<Type>(n, fill) {}
  Vector(const Vector<Type> & v) : Array<Type>(v) {}
  ~Vector() {}

  Vector<Type> operator-() const;
  Vector<Type> operator-(const Vector<Type> & v2) const;// subtraction
  Vector<Type> & operator-=(const Vector<Type> & v2);
  Vector<Type> operator+(const Vector<Type> & v2) const;// addition
  Vector<Type> & operator+=(const Vector<Type> & v2);
  Vector<Type> operator/(Real r) const;
  Vector<Type> & operator/=(Real r);

  Vector<Type> operator*(const Real r) const;
  Vector<Type> & operator*=(const Real r);

  const row_vector<Type> T() const; // transpose

  const Vector<Type> csub(uint i, uint j) const;
  const Vector<Type> sub(uint i, uint j) const {return csub(i,j);}

  // virtual Real length() const;
  // virtual Vector<Type> operator%(const Vector<Type> & v) const;
  // virtual Type operator^(const Vector<Type> & v2) const;
};

template<class Type>
Vector<Type> operator*(Real r, const Vector<Type> & v); // scaling

template<class Type>
class row_vector : public Array<Type>
{
friend class Vector<Type>;
private:
  const Vector<Type> to_column() const
    {
      const Array<Type> & phlerb = *this; // downcast
      return (Vector<Type>) phlerb; // upcast
    }
protected:
  row_vector() {}
  row_vector(const Array<Type> & v) : Array<Type>(v) {}
public:
  ~row_vector() {}

  const Vector<Type> T() const 
    {
      const Vector<Type> answer = this->to_column();
      return answer;
    }
};

template<class Type>
ostream & operator<<(ostream & os, const row_vector<Type> & v);

#ifdef __TCPLUSPLUS__
  #include "vector.cxx"
#endif

#endif
