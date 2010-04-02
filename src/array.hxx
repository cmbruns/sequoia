#ifndef ARRAY_HXX
#define ARRAY_HXX

// $Id: array.hxx,v 1.3 2001/12/12 18:26:27 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/array.hxx,v 1.3 2001/12/12 18:26:27 bruns Exp $
// $Log: array.hxx,v $
// Revision 1.3  2001/12/12 18:26:27  bruns
// Many changes to get compilation with gcc 3.0.2
// Removed all templated friend functions (no friends needed)
// Added .c_str() method to Mystring, and removed auto conversion to const char *
//
// Revision 1.2  2001/11/28 23:08:17  bruns
// Added cvs header tags
// Removed ^M characters
//

#include <iostream>
// exit()...
#ifdef IRIX
#include "stdlib.h"
#else
#include <stdlib.h>
#endif
#include "cmbmacro.hxx"

template<class Type>
class subArray;

template<class Type> // template so any type may be used
class Abstract_array  // abstract container base class
{
private:
  uint array_dim; // actual size used
  uint rep_dim; // size of representation array
  void Abstract_array<Type>::initialize_internals();
protected:
  uint check(uint n) const; // for use in bounds checking
  Type * element;
  void init_array(uint n);
  void init_array(uint n, Type val);
  const Type * array_ptr() const {return element;}
  Type * array_ptr() {return element;}
  void error(const char * msg) const;
public:
  Abstract_array<Type>(void) {initialize_internals();}
  Abstract_array<Type>(uint n) {initialize_internals(); init_array(n);}
  Abstract_array<Type>(uint n, Type val) {initialize_internals(); init_array(n, val);}
  Abstract_array<Type>(const Type * pn, uint n);
  Abstract_array<Type>(const Abstract_array<Type> & a);
  virtual ~Abstract_array();

  const Type & operator[](uint i) const {return element[i];}
  Type & operator[](uint i) {return element[i];}
  const Type &  operator()(uint i) const {return element[i];}
  Type &  operator()(uint i) {return element[i];}
  subArray<Type> sub(uint i, uint j);
  const Abstract_array<Type> sub(uint i, uint j) const;

  virtual Abstract_array<Type> & operator=(const Abstract_array<Type> & a);
  Abstract_array<Type> & operator+=(const Abstract_array<Type> & a);
  virtual Abstract_array<Type> & operator+=(const Type & a);

  uint dim() const {return array_dim;}
  void deletelast() {if (array_dim > 0) --array_dim;}
  void clear_array() {init_array(0);}
};


template<class T>
class Array : public Abstract_array<T>
{
public:
  Array(void) : Abstract_array<T>() {}
  Array(uint n) : Abstract_array<T>(n) {}
  Array(uint n, T val) : Abstract_array<T>(n, val) {}
  Array(const T * pn, uint n) : Abstract_array<T>(pn, n) {}
  Array(const Array<T> & a) : Abstract_array<T>(a) {}
};


template<class Type>
class subArray
{
friend class Abstract_array<Type>;
private:
  Abstract_array<Type> * parent;
  uint subarray_dim;
  uint offset;
  subArray() {parent = NULL; subarray_dim = 0; offset = 0;}
public:
  subArray<Type> & operator=(const Abstract_array<Type> & a);
  ~subArray() {}
};

// instaniation for Turbo C++
#ifdef __TCPLUSPLUS__
  #include "array.cxx"
#endif

#endif

