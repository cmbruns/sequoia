#include "cvector.h"

/* $Id: cvector.c,v 1.2 2001/11/28 23:08:17 bruns Exp $ */
/* $Header: /usr/data/cvs/sequoia1/cvector.c,v 1.2 2001/11/28 23:08:17 bruns Exp $ */
/* $Log: cvector.c,v $
/* Revision 1.2  2001/11/28 23:08:17  bruns
/* Added cvs header tags
/* Removed ^M characters
/* */

matrix EYE = {
{ 1,0,0 },
{ 0,1,0 },
{ 0,0,1 }
};

vector cross(const vector v1, const vector v2)
{
  vector answer;
  answer.x = v1.y * v2.z - v1.z * v2.y;
  answer.y = v1.z * v2.x - v1.x * v2.z;
  answer.z = v1.x * v2.y - v1.y * v2.x;
  return answer;
}

double dot(const vector v1, const vector v2)
{
  return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

double length(const vector v)
{
  return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}
 
double distance(const vector v1, const vector v2)
{
  return length(vecsubtract(v1,v2));
}
 
vector unit(const vector v)
{
  vector answer;
  double d = length(v);
  answer.x = v.x / d;
  answer.y = v.y / d;
  answer.z = v.z / d;
  return answer;
}

double vel(const vector v, int j)
{
  if (j == 0) return v.x;
  else if (j == 1) return v.y;
  else if (j == 2) return v.z;
  else 
    {
      printf("vector index out of bounds\n");
      vecprint(v);
      printf(", index = %d\n",j);
      exit(1);
    }
}

void vecprint(const vector v)
{
  printf("( %6.3f %6.3f %6.3f )",v.x,v.y,v.z);
}

vector vecscale(double r, const vector v)
{
  vector answer;
  answer.x = r * v.x;
  answer.y = r * v.y;
  answer.z = r * v.z;
  return answer;
}

vector vecsubtract(const vector v1, const vector v2)
{
  vector answer;
  answer.x = v1.x - v2.x;
  answer.y = v1.y - v2.y;
  answer.z = v1.z - v2.z;
  return answer;
}

/* get rotation angle from vector length from "axis" routine */
double rotangle(const vector v)
{
  double answer;
  answer = length(v) - 10;
  answer *= 180/PI;
  return answer;
}

vector get_vector()
{
  vector answer;
  float t1, t2, t3;
  scanf("%f %f %f",&t1,&t2,&t3);
  answer.x = t1;
  answer.y = t2;
  answer.z = t3;
  return answer;
}

vector vectoradd(const vector v1, const vector v2)
{
  vector answer;
  answer.x = v1.x + v2.x;
  answer.y = v1.y + v2.y;
  answer.z = v1.z + v2.z;
  return answer;
}

double mel(const matrix M, int i, int j)
{
  vector v = mvel(M,i);
  return vel(v,j);
}

vector mvel(const matrix M, int i)
{
  if (i == 0) return M.row1;
  else if (i == 1) return M.row2;
  else if (i == 2) return M.row3;
  else 
    {
      printf("matrix index out of bounds\n");
      matprint(M);
      printf(", index = %d\n",i);
      exit(1);
    }
}

vector mat_column(const matrix M, int col)
{
  vector answer;
  answer.x = mel(M,0,col);
  answer.y = mel(M,1,col);
  answer.z = mel(M,2,col);
  return answer;
}

void matprint(const matrix m)
{
  vecprint(m.row1);
  printf("\n");
  vecprint(m.row2);
  printf("\n");
  vecprint(m.row3);
  printf("\n");
}

matrix transpose(const matrix M)
{
  matrix answer;
  answer.row1 = mat_column(M,0);
  answer.row2 = mat_column(M,1);
  answer.row3 = mat_column(M,2);
  return answer;
}

matrix matmult(const matrix m1, const matrix m2)
{
  matrix answer;
  answer.row1.x = dot(m1.row1,mat_column(m2,0));
  answer.row1.y = dot(m1.row1,mat_column(m2,1));
  answer.row1.z = dot(m1.row1,mat_column(m2,2));
  answer.row2.x = dot(m1.row2,mat_column(m2,0));
  answer.row2.y = dot(m1.row2,mat_column(m2,1));
  answer.row2.z = dot(m1.row2,mat_column(m2,2));
  answer.row3.x = dot(m1.row3,mat_column(m2,0));
  answer.row3.y = dot(m1.row3,mat_column(m2,1));
  answer.row3.z = dot(m1.row3,mat_column(m2,2));
  return answer;
}

matrix matscale(double r, const matrix m2)
{
  matrix answer;
  answer.row1 = vecscale(r, m2.row1);
  answer.row2 = vecscale(r, m2.row2);
  answer.row3 = vecscale(r, m2.row3);
  return answer;
}

vector matvmult(const matrix M, const vector v)
{
  vector answer;
  answer.x = dot(M.row1,v);
  answer.y = dot(M.row2,v);
  answer.z = dot(M.row3,v);
  return answer;
}

/* return a rotation matrix equivalent to a rotation about a 
   particular axis */
matrix axis_angle(const vector v, double angle)
{
  matrix R1, R2, R3;
  vector test1, test2;

  /* first make a matrix that rotates v onto the x-axis */
  R1.row1 = unit(v);
  test1 = cross(R1.row1,EYE.row1);
  test2 = cross(R1.row1,EYE.row2);
  if (length(test1) > length(test2)) R1.row2 = unit(test1);
  else R1.row2 = unit(test2);
  R1.row3 = unit(cross(R1.row1,R1.row2));
  /* then rotate about x */
  R2 = EYE;
  R2.row2.y = cos(angle);
  R2.row3.z = R2.row2.y;
  R2.row2.z = sin(angle);
  R2.row3.y = -R2.row2.z;
  /* now rotate back */
  R3 = transpose(R1);
  return matmult(matmult(R3,R2),R1);
}

/* convert a rotation matrix into a rotation about a vector */
vector axis(const matrix M)
{
  vector v; /* rotation axis */
  vector diff[3],ihat[3],ihatprime[3];
  vector a,b;
  double len[3], angle, sinus, cosine;
  int big,littl,i;

  /* find displacements of the coordinate axes */
  for(i=0;i<3;++i)
    {
      ihat[i] = mat_column(EYE,i);
      ihatprime[i] = mat_column(M,i);
      diff[i] = vecsubtract(ihatprime[i],ihat[i]);
      len[i] = length(diff[i]);
    }
  /* find longest and shortest changes */
  big = 0; littl = 0;
  if (len[1] > len[0]) big = 1;
  else littl = 1;
  if (len[2] > len[big]) big = 2;
  if (len[2] < len[littl]) littl = 2;
  if (length(diff[big]) < 0.0001)
    {
       v.x = v.y = 0;
       v.z = 10;
       return v;
    }
  /* find rotation axis using two longest differences */
  if (littl == 0) v = cross(diff[1],diff[2]);
  else if (littl == 1) v = cross(diff[2],diff[0]);
  else v = cross(diff[0],diff[1]);
  v = unit(v);
  /* find rotation angle using the longest difference */
  /* a and b are projections of untransformed and transformed
     axis onto the plane perpendicular to rotation axis */
  a = unit(cross(v, cross(ihat[big],v)));
  b = unit(cross(v, cross(ihatprime[big],v)));
  sinus = dot(cross(b,a),v);
  cosine = dot(a,b);
  angle = acos(cosine);
  if (sinus < 0) angle = -angle;
  /* store angle in the vector itself */
  v = vecscale(10 + angle,v);
  return v;
}

matrix from_xeuler(const vector v1)
{
  matrix R1, R2, R3;
  R1 = axis_angle(EYE.row3, DEG2RAD * v1.x);
  R2 = axis_angle(EYE.row1, DEG2RAD * v1.y);
  R3 = axis_angle(EYE.row3, DEG2RAD * v1.z);
  return matmult(matmult(R3,R2),R1);
}

matrix from_ceuler(const vector v1)
{
  matrix R1, R2, R3;
  R1 = axis_angle(EYE.row3, DEG2RAD * -v1.x);
  R2 = axis_angle(EYE.row2, DEG2RAD * -v1.y);
  R3 = axis_angle(EYE.row3, DEG2RAD * -v1.z);
  return matmult(matmult(R1,R2),R3);
}

/* X-PLOR Euler angles */
vector to_xeuler(const matrix M)
{
  vector answer;

  answer.x = atan2(M.row3.x,-M.row3.y) * 180/PI;
  answer.y = acos(M.row3.z) * 180/PI;
  answer.z = atan2(M.row1.z,M.row2.z) * 180/PI;
  if (answer.x < 0) answer.x += 360;
  if (answer.y < 0) answer.y += 360;
  if (answer.z < 0) answer.z += 360;

  return answer;
}
 
/* CCP4 Euler angles */
vector to_ceuler(const matrix M)
{
  vector answer;

  if (M.row3.z > 0.9999) 
    {
       answer.y = 0;
       answer.z = 0; 
       answer.x = atan2(M.row2.x,M.row1.x) * 180/PI; 
    }
  else
    {
      answer.x = atan2(M.row2.z,M.row1.z) * 180/PI;
      answer.y = acos(M.row3.z) * 180/PI; /* OK */
      answer.z = atan2(M.row3.y,-M.row3.x) * 180/PI;
    }
  if (answer.x < 0) answer.x += 360;
  if (answer.y < 0) answer.y += 360;
  if (answer.z < 0) answer.z += 360;

  return answer;
}
 
homog homogmult(const homog h1, const homog h2)
{
  homog answer;
  answer.rot = matmult(h1.rot,h2.rot);
  answer.trans = matvmult(h1.rot,h2.trans);
  answer.trans = vectoradd(answer.trans,h1.trans);
  return answer;
}

homog v2h(const vector t)
{
  homog answer;
  answer.rot = EYE;
  answer.trans = t;
  return answer;
}

homog m2h(const matrix M)
{
  homog answer;
  answer.rot = M;
  answer.trans.x = 0;
  answer.trans.y = 0;
  answer.trans.z = 0;
  return answer;
}
 
homog inverse(const homog H)
{
  homog answer;
  answer.rot = transpose(H.rot);
  answer.trans = vecscale(-1, H.trans);
  answer.trans = matvmult(answer.rot,answer.trans);
  return answer;
}
