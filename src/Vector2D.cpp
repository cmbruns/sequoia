/* Copyright (c) 2005 Christopher M. Bruns
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
 
// 
// $Id$
//
// $Header$
//
// $Log$
// Revision 1.1  2004/06/04 19:34:48  cmbruns
// Imported structure related sources from archive on baxter
// Debugged simple conversion of structures to sequences.
// Implemented computation of solvent accessible surface areas
// Created target residue_area, for output of residue solvent accessible surfaces areas
// Updated GPL headers
//
// Revision 1.3  2002/09/14 00:02:51  bruns
// Added license header to most .cc files
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#include "Vector2D.h"

#define EPSILON FLT_EPSILON
#define ABS(x) (((x) < 0) ? (-(x)) : (x))

Vector2D::Vector2D() {
  int i;
  for(i = 0; i < 2; i ++)
    private_element[i] = 0;
}
int Vector2D::get_size() const {return 2;}

const double & Vector2D::operator[](const int index) const {return private_element[index];}
double & Vector2D::operator[](const int index) {return private_element[index];}
Vector2D Vector2D::operator+(const Vector2D & v2) const {
  Vector2D answer = *this;
  int i;
  for(i = 0; i < get_size(); i ++)
    answer[i] += v2[i];
  return answer;
}
Vector2D Vector2D::operator-(const Vector2D & v2) const {
  Vector2D answer = *this;
  int i;
  for(i = 0; i < get_size(); i ++)
    answer[i] -= v2[i];
  return answer;
}
Vector2D Vector2D::operator*(double r) const {
  Vector2D answer = *this;
  int i;
  for(i = 0; i < get_size(); i ++)
    answer[i] *= r;
  return answer;
}
double Vector2D::dot(const Vector2D & v2) const {
  const Vector2D & v1 = *this;
  double answer = 0.0;
  int i;
  for(i = 0; i < get_size(); i ++)
    answer += v1[i] * v2[i];
  return answer;
}

Vector2D Vector2D::unit() const {
  const Vector2D v = *this;
  double scale = v[0] * v[0] +
    v[1] * v[1] +
    v[2] * v[2];
  if (ABS(1.0 - scale) < EPSILON) return v;
  else {
    scale = sqrt(1.0 / scale);
    Vector2D answer = v * scale;
    return answer;
  }
}

const double Vector2D::distance(const Vector2D & v2) const {
  const Vector2D & v1 = *this;
  Vector2D delta = v2 - v1;
  return sqrt(delta.dot(delta));
}

Vector2D givens(const double a, const double b) {
  Vector2D answer;
  double tau, c, s;

  if (b == 0) {
    c = 1.0;
    s = 0.0;
  }
  else {
    if (ABS(b) > ABS(a)) {
      tau = -a / b;
      s = 1.0 / sqrt(1 + tau*tau);
      c = s * tau;
    }
    else {
      tau = -b / a;
      c = 1.0 / sqrt(1.0 + tau*tau);
      s = c * tau;
    }
  }
  
  answer[0] = c;
  answer[1] = s;
  return answer;
}

