#include "cmbmacro.hxx"

// $Id: cmbmacro.cxx,v 1.2 2001/11/28 23:08:17 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/cmbmacro.cxx,v 1.2 2001/11/28 23:08:17 bruns Exp $
// $Log: cmbmacro.cxx,v $
// Revision 1.2  2001/11/28 23:08:17  bruns
// Added cvs header tags
// Removed ^M characters
//

template<class T>
void swap(T & a, T & b)
{
	T   temp = a;
	a = b;
	b = temp;
}

template<class Type>
Type MIN(const Type & v1, const Type & v2)
{
  if (v1 > v2) return v2;
  else return v1;
}

template<class Type>
Type MAX(const Type & v1, const Type & v2)
{
  if (v1 > v2) return v1;
  else return v2;
}

template<class Type>
Type ABS(const Type & v1)
{
  if (v1 >= 0) return v1;
  else return -v1;
}

template<class Type>
int sign(const Type & v1)
{
  if (v1 >= 0) return 1;
  else return -1;
}

#ifdef __TCPLUSPLUS__
uint MIN(const uint &, const uint &);
double ABS(const double & d);
#endif

