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
// Revision 1.5  2002/09/13 23:28:08  bruns
// Added license header to all header files
//
// Revision 1.4  2002/09/13 22:30:59  bruns
// Added negation operator
//
// Revision 1.3  2001/12/14 19:21:55  bruns
// Changes to compile with gcc 3.0.2
//  - "using namespace std;"  to use standard classes
//  - "ios::fmtflags" instead of "fmtflags"
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#ifndef _QUATERNION_H_
#define _QUATERNION_H_

using namespace std;

class Quaternion;

#include "Matrix3D.h"
#include "Matrix4D.h"
#include "Vector4D.h"

// Try to make sure that the Quaternion is a unit quaternion

class Quaternion {
private:
  double private_l; /* lambda, or x */
  double private_m; /* mu, or y */
  double private_n; /* nu, or z */
  double private_s; /* sigma, or w */
public:
  Quaternion();
  Quaternion(double l, double m, double n, double s);
  Quaternion(const Matrix3D & M);
  Quaternion(const Matrix4D & M);

  double get_l() const {return private_l;}
  double get_m() const {return private_m;}
  double get_n() const {return private_n;}
  double get_s() const {return private_s;}
  double get_lambda() const {return private_l;}
  double get_mu() const {return private_m;}
  double get_nu() const {return private_n;}
  double get_sigma() const {return private_s;}
  double get_x() const {return private_l;}
  double get_y() const {return private_m;}
  double get_z() const {return private_n;}
  double get_w() const {return private_s;}

  void set(double l, double m, double n, double s) {
    private_l = l;
    private_m = m;
    private_n = n;
    private_s = s;
  }
  void set(const Vector4D & v) {
    private_l = v[0];
    private_m = v[1];
    private_n = v[2];
    private_s = v[3];
  }

  double & set_l()  {return private_l;}
  double & set_m()  {return private_m;}
  double & set_n()  {return private_n;}
  double & set_s()  {return private_s;}
  double & set_lambda()  {return private_l;}
  double & set_mu()  {return private_m;}
  double & set_nu()  {return private_n;}
  double & set_sigma()  {return private_s;}
  double & set_x()  {return private_l;}
  double & set_y()  {return private_m;}
  double & set_z()  {return private_n;}
  double & set_w()  {return private_s;}

  double & set_l(double d)  {private_l = d; return private_l;}
  double & set_m(double d)  {private_m = d; return private_m;}
  double & set_n(double d)  {private_n = d; return private_n;}
  double & set_s(double d)  {private_s = d; return private_s;}
  double & set_lambda(double d)  {private_l = d; return private_l;}
  double & set_mu(double d)  {private_m = d; return private_m;}
  double & set_nu(double d)  {private_n = d; return private_n;}
  double & set_sigma(double d)  {private_s = d; return private_s;}
  double & set_x(double d)  {private_l = d; return private_l;}
  double & set_y(double d)  {private_m = d; return private_m;}
  double & set_z(double d)  {private_n = d; return private_n;}
  double & set_w(double d)  {private_s = d; return private_s;}

  Quaternion operator-() const { // inverse - only works for unit quaternions now TODO
    Quaternion answer = *this;
    answer.set_l(- answer.get_l());
    answer.set_m(- answer.get_m());
    answer.set_n(- answer.get_n());
    return answer;
  }
  Quaternion operator*(const Quaternion & q2) const; // multiplication
  Quaternion unit() const;
  bool is_unit() const;
  Matrix3D rotmat() const;
};

ostream & operator<<(ostream & os, const Quaternion & q);

#endif


