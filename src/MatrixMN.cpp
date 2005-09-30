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
// Revision 1.4  2002/09/14 00:02:51  bruns
// Added license header to most .cc files
//
// Revision 1.3  2002/09/13 22:46:43  bruns
// removed redundant parameter default
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#include "MatrixMN.h"
#include "Vector2D.h"

#define SIGN(x) (((x) < 0) ? -1.0 : 1.0)

int MatrixMN::get_n_rows() const {return private_element.size();}

int MatrixMN::get_n_columns() const {
  if (get_n_rows() == 0) return 0;
  else {
    const MatrixMN & M = *this;
    return M[0].get_size();
  }
}

VectorND & MatrixMN::operator[](int index) { // Non-const version  
  return private_element[index];
}

const VectorND & MatrixMN::operator[](int index) const { // const version  
  return private_element[index];
}

VectorND MatrixMN::operator*(const VectorND & v) const {
  int i, j;
  VectorND answer(get_n_rows());
  const MatrixMN & M = *this;
  
  for (i = 0; i < M.get_n_rows(); ++i) {
    answer[i] = 0;
    for (j = 0; j < M.get_n_columns(); ++j)
      answer[i] += M[i][j] * v[j];
  }
  return answer;
}

MatrixMN MatrixMN::operator+(const MatrixMN & m2) const {
  MatrixMN answer = *this;
  int i, j;

  for (i = 0; i < get_n_rows(); ++i)
    for (j = 0; j < get_n_columns(); ++j)
      answer[i][j] += m2[i][j];
  
  return answer;
}

MatrixMN MatrixMN::transpose() const {
  MatrixMN answer(get_n_columns(), get_n_rows());
  int i, j;
  const MatrixMN m = *this;
  
  for (i = 0; i < get_n_rows(); ++i)
    for (j = 0; j < get_n_columns(); ++j)
      answer[j][i] = m[i][j];
  
  return answer;
}

MatrixMN MatrixMN::sym_QRstep(int p, BaseMatrix * Q) const {
  const MatrixMN & T = *this;
  MatrixMN answer = *this;
  Vector2D cs;
  double d, mu, x, z;
  int k;
  int n;

  n = T.get_n_columns();

  d = (answer[n-2][n-2] - answer[n-1][n-1]) / 2.0;

  mu = answer[n-1][n-1] - (answer[n-1][n-2]) * (answer[n-1][n-2]) /
    (d + SIGN(d) * sqrt(d*d + (answer[n-1][n-2]) * (answer[n-1][n-2])));

  x = answer[0][0] - mu;
  z = answer[1][0];
  
  for (k = 0; k < (n-1); ++k) {
    cs = givens(x,z);
    answer.row_rot(k, k+1, cs);
    answer.col_rot(k, k+1, cs);
    if (Q != NULL) {
      Q->col_rot(k+p, k+1+p, cs);
    }

    if (k < (n - 2)) {
      x = answer[k+1][k];
      z = answer[k+2][k];
    }
  }

  return answer;
}

/* left Given's rotation - G&vL alg. 5.1.6 p. 203 */
void MatrixMN::row_rot(const int i, const int k, const Vector2D & cs) {
  const MatrixMN & A = *this;
  MatrixMN answer = *this;
  double c, s;
  double tau1, tau2;
  int j;

  c = cs[0];
  s = cs[1];
  for (j = 0; j < A.get_n_columns(); ++j) {
    tau1 = answer[i][j];
    tau2 = answer[k][j];
    answer[i][j] = c*tau1 - s*tau2;
    answer[k][j] = s*tau1 + c*tau2;
  }

  // return answer;
  *this = answer;
}

/* right Given's rotation - G&vL alg. 5.1.7 p. 203 */
void MatrixMN::col_rot(const int j, const int k, const Vector2D & cs) {
  const MatrixMN & A = *this;
  MatrixMN answer = *this;
  double c, s;
  double tau1, tau2;
  int i;

  c = cs[0];
  s = cs[1];
  for (i = 0; i < A.get_n_rows(); ++i) {
    tau1 = answer[i][j];
    tau2 = answer[i][k];
    answer[i][j] = c*tau1 - s*tau2;
    answer[i][k] = s*tau1 + c*tau2;
  }

  // return answer;
  *this = answer;
}

