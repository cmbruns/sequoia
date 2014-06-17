#ifndef CMBMACRO_HXX
#define CMBMACRO_HXX

// $Id: cmbmacro.hxx,v 1.2 2001/11/28 23:08:17 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/cmbmacro.hxx,v 1.2 2001/11/28 23:08:17 bruns Exp $
// $Log: cmbmacro.hxx,v $
// Revision 1.2  2001/11/28 23:08:17  bruns
// Added cvs header tags
// Removed ^M characters
//

typedef unsigned char uchar;
typedef unsigned char ubyte;
// assumes float is 4 bytes!!!
typedef float real4;
typedef double Real;
typedef unsigned int uint;

#ifdef __TCPLUSPLUS__
 typedef int bool;
#endif

# ifndef TRUE
#  define TRUE 1
#  define FALSE 0
# endif

# ifndef true
#  define true 1
#  define false 0
# endif

# ifndef NULL
#  define NULL 0
# endif

#include<math.h>
#ifndef PI
# define PI 3.1415926536
#endif

#ifdef __TCPLUSPLUS__
#include "cmbmacro.cxx"
#endif

template<class T>
void swap(T & a, T & b);

template<class Type>
Type MIN(const Type & v1, const Type & v2);

template<class Type>
Type MAX(const Type & v1, const Type & v2);

template<class Type>
Type ABS(const Type & v1);

template<class Type>
int sign(const Type & v1);

#ifndef DEG2RAD
#define DEG2RAD 0.017453293
#endif

#ifndef RAD2DEG
#define RAD2DEG 57.295779513
#endif

#endif
// cmbmacro.hxx
