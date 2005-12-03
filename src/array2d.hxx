#ifndef _ARRAY2D_HXX_
#define _ARRAY2D_HXX_

#include "array.hxx"

// $Id: array2d.hxx,v 1.2 2001/11/28 23:08:17 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/array2d.hxx,v 1.2 2001/11/28 23:08:17 bruns Exp $
// $Log: array2d.hxx,v $
// Revision 1.2  2001/11/28 23:08:17  bruns
// Added cvs header tags
// Removed ^M characters
//

template <class Type>
class Array2d : public Array< Array<Type> >
{
protected:
  Array2d(int n) : Array< Array<Type> > (n) {}
public:
  Array2d(void) ;
  Array2d(uint m, uint n);
  Array2d(uint m, uint n, Type val);
  Array2d(const Array2d<Type> & M);
  ~Array2d();
  void init_2darray(uint r, uint c);
  void init_2darray(uint r, uint c, Type val);

  uint m() const {return dim();}
  uint n() const 
    {
      if (!m()) return 0;
      else return this->operator[](0).dim();
    }
};

#ifdef __TCPLUSPLUS__
 #include "array2d.cxx"
#endif

#endif

