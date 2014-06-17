#include "array2d.hxx"

// $Id: array2d.cxx,v 1.2 2001/11/28 23:08:17 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/array2d.cxx,v 1.2 2001/11/28 23:08:17 bruns Exp $
// $Log: array2d.cxx,v $
// Revision 1.2  2001/11/28 23:08:17  bruns
// Added cvs header tags
// Removed ^M characters
//

template <class Type>
Array2d<Type>::Array2d() : Array< Array<Type> > ()
{
}

template <class Type>
Array2d<Type>::Array2d(uint r, uint c, Type val) 
{
  init_2darray(r, c, val);
}

template <class Type>
Array2d<Type>::Array2d(uint r, uint c)
{
  init_2darray(r, c);
}

template <class Type>
Array2d<Type>::Array2d(const Array2d<Type> & M)
{
  *this = M;
}

template <class Type>
Array2d<Type>::~Array2d()
{
}

//template <class Type>
//Array2d<Type> & Array2d<Type>::operator=(const Array2d<Type> & M)
//{
//  *this = (const Array< Array<Type> >) M;
//  return *this;
//}

template <class Type>
void Array2d<Type>::init_2darray(uint r, uint c)
{
  Array<Type> bit(c);
  init_array(r, bit);
}

template <class Type>
void Array2d<Type>::init_2darray(uint r, uint c, Type val)
{
  Array<Type> bit(c, val);
  init_array(r, bit);
}

