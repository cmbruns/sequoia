#ifndef _LTABLE_HXX_
#define _LTABLE_HXX_

// $Id: ltable.hxx,v 1.2 2001/11/28 23:08:17 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/ltable.hxx,v 1.2 2001/11/28 23:08:17 bruns Exp $
// $Log: ltable.hxx,v $
// Revision 1.2  2001/11/28 23:08:17  bruns
// Added cvs header tags
// Removed ^M characters
//

#include <math.h>

#define USE_TABLES 1

#define SQR2 1.414213562

#ifndef PI
#define PI 3.141592654
#endif

double table_asin(double x);
double table_acos(double x);
double table_atan2(double y, double x);
double table_log(double x);
double table_sqrt(double x);

double zprob(double x);
double table_zprob(double x);
double angprob(double x);
double table_angprob(double x);
double angprobC(double x);
double table_angprobC(double x);
double radprob(double x);
double table_radprob(double x);

#ifdef USE_TABLES
#define ZPROB(x) table_zprob(x)
#define SQRT(x) table_sqrt(x)
#define ASIN(x) table_asin(x)
#define ACOS(x) table_acos(x)
#define ATAN2(x,y) table_atan2(x,y)
#define LOG(x) table_log(x)
#define ANGPROB(x) table_angprob(x)
#define ANGPROBC(x) table_angprobC(x)
#else
#define ZPROB(x) zprob(x)
#define SQRT(x) sqrt(x)
#define ASIN(x) asin(x)
#define ACOS(x) acos(x)
#define ATAN2(x,y) atan2(x,y)
#define LOG(x) log(x)
#define ANGPROB(x) angprob(x)
#define ANGPROBC(x) angprobC(x)
#endif

#endif
