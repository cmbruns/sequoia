# include "array.hxx"

#ifndef _ARRAY_CXX_
#define _ARRAY_CXX_

// $Id: array.cxx,v 1.2 2001/11/28 22:38:51 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/array.cxx,v 1.2 2001/11/28 22:38:51 bruns Exp $
// $Log: array.cxx,v $
// Revision 1.2  2001/11/28 22:38:51  bruns
// Removed ^M from ends of lines
// Added cvs tags
//

template <class Type>
void Abstract_array<Type>::initialize_internals()
{
  // cerr << "Initialize internals" << endl;
  array_dim = 0;
  rep_dim = 0;
  element = NULL;
}

template <class Type>
void Abstract_array<Type>::init_array(uint n)
{
  if (!n);
  else if (!rep_dim) // new structure
    {
      rep_dim = n;
      element = new Type[rep_dim];
    }
  else if (n > rep_dim) // needs to be enlarged
    {
      Abstract_array<Type> store = *this;
      delete [] element;
      if (n <= (2 * rep_dim)) rep_dim = (2 * rep_dim);
      else rep_dim = n;
      element = new Type[rep_dim];
      *this = store;
    }
  array_dim = n;
}

template <class Type>
void Abstract_array<Type>::init_array(uint n, Type val)
{
  init_array(n);
  uint i;
  Abstract_array<Type> & a = *this;
  for (i=0; i < dim(); ++i)
    a[i] = val;
}

template<class Type>
Abstract_array<Type>::Abstract_array(const Type * pn, uint n)
{
  initialize_internals();
  init_array(n);
  uint i;
  Abstract_array<Type> & a = *this;
  for (i=0; i < dim(); ++i)
    a[i] = pn[i];
}

template<class T>
Abstract_array<T> & 
 Abstract_array<T>::operator=(const Abstract_array<T> & a)
{
  if (this == &a) return *this;
  Abstract_array & a1 = *this;
  init_array(a.dim());
  uint i;
  for (i = 0; i < dim(); ++i)
    a1[i] = a[i];
  return *this;
}

template<class T>
Abstract_array<T>::Abstract_array(const Abstract_array<T> & a)
{
  initialize_internals();
  *this = a;
}

template<class Type>
Abstract_array<Type>::~Abstract_array(void)
{
  rep_dim = 0;
  array_dim = 0;
  delete [] element;
  element = NULL;
}

template<class Type>
void  Abstract_array<Type>::error(const char * msg) const
{
  cerr << "Array error\n";
  cerr << "Array dimension = " << dim() << "\n";
  cerr << msg << "\n";
  exit(1);
}


// bounds checking function
template<class T>
uint Abstract_array<T>::check(uint i) const
{
//  if ((i < 0) || (i >= dim()))  // unsigned
  if (i >= dim())
    error("Array dimension out of bounds!!");
  return i;
}


template<class Type>
const Abstract_array<Type> Abstract_array<Type>::sub(uint i, uint j) const
{
  const Abstract_array<Type> & a = *this;
  Abstract_array<Type> answer (j - i);
  uint k;
  for (k=0; k < answer.dim(); ++k)
    answer[k] = a[k + i];
  return answer;
}


template<class Type>
subArray<Type> Abstract_array<Type>::sub(uint i, uint j)
{
  subArray<Type> answer;
  answer.offset = i;
  answer.subarray_dim = j - i;
  answer.parent = this;
  return answer;
}




template<class T>
Abstract_array<T> & 
 Abstract_array<T>::operator+=(const Abstract_array<T> & a)
{
  if (this == &a) return *this;
  Abstract_array<T> & a1 = *this;
  uint i;
  init_array(dim() + a.dim());
  for (i = (dim() - a.dim()); i < dim(); ++i)
    a1[i] = a[i - dim() + a.dim()];
  return *this;
}

template<class T>
Abstract_array<T> & 
 Abstract_array<T>::operator+=(const T & a)
{
  Abstract_array & a1 = *this;
  init_array(dim() + 1);
  a1[dim() - 1] = a;
  return *this;
}


template<class T> 
ostream & operator<<(ostream & os, const Array<T> & a)
{
  os << "(";
  int i;
  for (i = 0; i < (a.dim() - 1); ++i)
    os << a[i] << ",";
  os <<  a[a.dim() - 1] << ")";
  return os;
}


template<class Type> 
subArray<Type> & subArray<Type>::operator=(const Abstract_array<Type> & a)
{
  Abstract_array<Type> & a1 = *parent;
  uint i;
  for (i=0; i < subarray_dim; ++i)
    a1[i + offset] = a[i];
  return *this;
}

// template class Abstract_array<int>;

#endif
