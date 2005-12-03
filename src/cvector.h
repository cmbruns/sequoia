#ifndef _CVECTOR_H_
#define _CVECTOR_H_

/* $Id: cvector.h,v 1.2 2001/11/28 23:08:17 bruns Exp $ */
/* $Header: /usr/data/cvs/sequoia1/cvector.h,v 1.2 2001/11/28 23:08:17 bruns Exp $ */
/* $Log: */

#include <stdio.h>
#include <math.h>

#ifndef PI
#define PI 3.14159265359
#endif

#define DEG2RAD PI/180
#define RAD2DEG 180/PI

typedef struct vector_struct
{
  float x;
  float y;
  float z;
} vector;

typedef struct matrix_struct
{
  vector row1;
  vector row2;
  vector row3;
} matrix;

extern matrix EYE;

typedef struct homog_struct
{
  matrix rot;
  vector trans;
}homog;

vector cross(const vector v1, const vector v2);
double dot(const vector v1, const vector v2);
double length(const vector v);
double distance(const vector v1, const vector v2);
vector unit(const vector v);
double vel(const vector v, int j);
void vecprint(const vector v);
vector vecscale(double r, const vector v);
double rotangle(const vector v);
vector get_vector();
vector vectoradd(const vector v1, const vector v2);
vector vecsubtract(const vector v1, const vector v2);

matrix matmult(const matrix m1, const matrix m2);
matrix matadd(const matrix m1, const matrix m2);
matrix matsubtract(const matrix m1, const matrix m2);
matrix matscale(double r, const matrix m2);
vector matvmult(const matrix M, const vector v);
matrix axis_angle(const vector v, double angle);
vector axis(const matrix m);
vector to_xeuler(const matrix m1);
matrix from_xeuler(const vector v1);
vector to_ceuler(const matrix m1);
matrix from_ceuler(const vector v1);
double mel(const matrix m, int i, int j);
vector mvel(const matrix m, int i);
void matprint(const matrix m);
matrix transpose(const matrix M);

homog homogmult(const homog h1, const homog h2);
homog v2h(const vector v);
homog m2h(const matrix M);
homog inverse(const homog H);

#endif
