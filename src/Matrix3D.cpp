// This file is part of the Sequoia package for macromolecular 
//  sequence/structure analysis
// Copyright (C) 2004  Christopher M. Bruns, Ph.D.
// 
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
// 
// See the accompanying file 'LICENSE' for details
// 
// To contact the author, write to cmbruns@comcast.net or bruns@scripps.edu
// In publications please cite: Bruns et al (1999), J.Mol.Biol. 288:427-439
// Please submit bug reports at http://bruns.homeip.net/bugzilla/
// 
// To obtain a non-GPL version of this program, see http://bruns.homeip.net/sequoia.html
// 
// $Id$
//
// $Header$
//
// $Log$
// Revision 1.1  2004/06/04 19:34:47  cmbruns
// Imported structure related sources from archive on baxter
// Debugged simple conversion of structures to sequences.
// Implemented computation of solvent accessible surface areas
// Created target residue_area, for output of residue solvent accessible surfaces areas
// Updated GPL headers
//
// Revision 1.5  2002/09/14 00:02:51  bruns
// Added license header to most .cc files
//
// Revision 1.4  2002/09/13 22:46:43  bruns
// removed redundant parameter default
//
// Revision 1.3  2002/05/15 19:28:46  bruns
// Added axis_angle method
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#include "Matrix3D.h"
#include "Quaternion.h"
#include "MatrixMN.h"
#include "Vector2D.h"

/* need an estimate of unit roundoff error for real */
#include <climits>
#define EPSILON FLT_EPSILON
#define ABS(x) (((x) < 0) ? (-(x)) : (x))

Matrix3D::Matrix3D() {
  // Identity matrix
  int i, j;
  for (i = 0; i < 3; i ++) 
    for (j = 0; j < 3; j ++) {
      if (i == j) private_row[i][j] = 1.0; // Diagonal
      else private_row[i][j] = 0.0; // Off-diagonal
    }
}

Matrix3D::Matrix3D(double e00, double e01, double e02,
	 double e10, double e11, double e12,
	 double e20, double e21, double e22) {
  // Identity matrix
  private_row[0][0] = e00;
  private_row[0][1] = e01;
  private_row[0][2] = e02;
  private_row[1][0] = e10;
  private_row[1][1] = e11;
  private_row[1][2] = e12;
  private_row[2][0] = e20;
  private_row[2][1] = e21;
  private_row[2][2] = e22;
}

Matrix3D & Matrix3D::set(const Vector3D & v1, 
			 const Vector3D & v2,
			 const Vector3D & v3) {
  private_row[0][0] = v1[0];
  private_row[0][1] = v1[1];
  private_row[0][2] = v1[2];
  private_row[1][0] = v2[0];
  private_row[1][1] = v2[1];
  private_row[1][2] = v2[2];
  private_row[2][0] = v3[0];
  private_row[2][1] = v3[1];
  private_row[2][2] = v3[2];
  return *this;
}

const Vector3D & Matrix3D::operator[](int index) const {
  return private_row[index];
}

Vector3D & Matrix3D::operator[](int index) {
  return private_row[index];
}

Matrix3D Matrix3D::operator*(const Matrix3D & m2) const {
  int i, j, k;
  Matrix3D answer;
  const Matrix3D & m1 = *this;

  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j) {
      answer[i][j] = 0;
      for (k = 0; k < 3; ++k)
	answer[i][j] += m1[i][k] * m2[k][j];
    }
  return answer;
}

Vector3D Matrix3D::operator*(const Vector3D & v) const {
  int i, j;
  Vector3D answer;
  const Matrix3D & M = *this;

  for (i = 0; i < 3; ++i) {
    answer[i] = 0;
    for (j = 0; j < 3; ++j)
      answer[i] += M[i][j] * v[j];
  }
  return answer;
}

Matrix3D Matrix3D::operator*(const double scale) const {
  int i, j;
  Matrix3D answer = *this;

  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      answer[i][j] *= scale;
  return answer;
}

Matrix3D Matrix3D::operator+(const Matrix3D & m2) const {
  Matrix3D answer;
  int i, j;
  const Matrix3D m1 = *this;

  answer = m1;
  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      answer[i][j] += m2[i][j];
  
  return answer;
}

Matrix3D Matrix3D::operator-(const Matrix3D & m2) const {
  Matrix3D answer;
  int i, j;
  const Matrix3D & m1 = *this;

  answer = m1;
  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      answer[i][j] -= m2[i][j];
  
  return answer;
}

Matrix3D Matrix3D::scale(double scale) const {
  Matrix3D answer;
  int i, j;
  const Matrix3D & M = *this;

  answer = M;
  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      answer[i][j] = scale * M[i][j];

  return answer;
}

double Matrix3D::trace() const {
  const Matrix3D & M = *this;
  return M[0][0] + M[1][1] + M[2][2];
}

Matrix3D Matrix3D::transpose() const {
  Matrix3D answer;
  int i, j;
  const Matrix3D m = *this;

  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      answer[i][j] = m[j][i];
  
  return answer;
}

Quaternion Matrix3D::quaternion() const {
  Quaternion answer;
  double  tr, s, q[4];
  int    i, j, k;
  const Matrix3D & M = *this;
  
  int nxt[3] = {1, 2, 0};
  
  tr = M.trace();
  
  /* check the diagonal */
  if (tr > 0.0) {
    s = sqrt (tr + 1.0);
    answer.set_s() = s / 2.0;
    s = 0.5 / s;
    answer.set_l() = (M[1][2] - M[2][1]) * s;
    answer.set_m() = (M[2][0] - M[0][2]) * s;
    answer.set_n() = (M[0][1] - M[1][0]) * s;
  } 
  else {                
    /* diagonal is negative */
    i = 0;
    if (M[1][1] > M[0][0]) i = 1;
    if (M[2][2] > M[i][i]) i = 2;
    j = nxt[i];
    k = nxt[j];
    
    s = sqrt ((M[i][i] - (M[j][j] + M[k][k])) + 1.0);
    
    q[i] = s * 0.5;
    
    if (s != 0.0) s = 0.5 / s;
    
    q[3] = (M[j][k] - M[k][j]) * s;
    q[j] = (M[i][j] + M[j][i]) * s;
    q[k] = (M[i][k] + M[k][i]) * s;
    
    answer.set_l() = q[0];
    answer.set_m() = q[1];
    answer.set_n() = q[2];
    answer.set_s() = -q[3];
  }

  return answer;
}

MatrixMN Matrix3D::get_submatrix(const int m1, const int m2, const int n1, const int n2) const {
	MatrixMN answer(m2 - m1 + 1, n2 - n1 + 1);
	int i;
	// int j;
	const Matrix3D & m = *this;
  
	int n_rows = m2 - m1 + 1;
  
	for (i = 0; i < n_rows; i++ )
		answer[i] = m[i + m1].subvector(n1,n2);
  
	return answer;
}

Matrix3D & Matrix3D::set_submatrix(const MatrixMN & M2, int m, int n) {
  int i, j;
  Matrix3D & M1 = *this;
  
  for (i = 0; i < M2.get_n_rows(); ++i)
    for (j = 0; j < M2.get_n_columns(); ++j)
      M1[i + m][j + n] = M2[i][j];

  return *this;
}

Vector3D Matrix3D::top_eigenvector() const {
  const Matrix3D & P = *this;
  Matrix3D eigenvectors;
  Vector3D eigenvalues;
  Vector3D answer;
  int i, j;

  // eigenvectors = eye(P.n);
  eigenvalues = P.sym_eigen(&eigenvectors);

  /* find top eigenvalue */
  j = 0; /* guess */
  for (i = 1; i < 3; ++i)
    if (eigenvalues[i] > eigenvalues[j]) j = i;
  eigenvectors = eigenvectors.transpose();
  answer = eigenvectors[j]; /* top eigenvector */
  
  return answer;
}

/* symmetric matrix eigenvector/value decomposition */
  // FIXME
/* eigenvalues and eigenvectors of a symmetric matrix */
/* algorithm 8.2.3 of G&vL (p. 423) - symmetric QR algorithm */
Vector3D Matrix3D::sym_eigen(Matrix3D * Q) const {
  const Matrix3D & A = *this;
  Vector3D answer;
  /* Q is the matrix of eigenvectors */
  /* D is a diagonal matrix of eigenvalues */
  /* epsilon is a measure of the unit roundoff */

  Matrix3D T;
  MatrixMN T22;
  int i, j;
  int p, q, n;
  long double epsilon;
  
  /* 1 - form tridiagonal matrix by G&vL 8.2.1 (p. 420) */
  n = 3;
  /*  (Householder tridiagonalization) */
  T = A.sym_tridiag(Q);

  /* 2 - zero the off-diagonal elements */
  do {
    /* set near-zeros to zero */
    for (i = 0; i < n - 1; ++i) {
      epsilon = EPSILON * (ABS(T[i][i]) + ABS(T[i+1][i+1]));
      if (ABS(T[i+1][i]) < epsilon) {
	T[i+1][i] = 0.0;
	T[i][i+1] = 0.0;
      }
      if (ABS(T[i][i+1]) < epsilon) {
	T[i+1][i] = 0.0;
	T[i][i+1] = 0.0;
      }
      /* keep it strictly pentadiagonal */
      for (j = i + 3; j < n; ++j) {
	T[i][j] = 0.0;
	T[j][i] = 0.0;
      }
    }
    /* Find largest q... */
    q = 0;
    i = n - 1;
    if ((i >= 1) && (T[i-1][i] == 0)) {
      ++q;
      --i;
      while ((i >= 1) && (T[i-1][i] == 0)) {
	++q;
	--i;
      }
      if (i <= 0) {q = n; i = 0;}
    }

    /* and the smallest p... */
    p = n - q - 1;
    if ((i >= 1) && (T[i-1][i] != 0)) {
      --p;
      --i;
      while ((i >= 1) && (T[i-1][i] != 0)) {
	--p;
	--i;
      }
      if (i <= 0) {p = 0; i = 0;}
    }
    if (p < 0) p = 0;

    if (q < n) {
      /* Apply Algorithm 8.2.2 to T22 */
      T22 = T.get_submatrix(p, n-1-q, p, n-1-q);
      T22 = T22.sym_QRstep(p, Q);
      T.set_submatrix(T22, p, p);
    }
  } while (q < n);
  
  for (i = 0; i < T.get_n_columns(); ++i)
    answer[i] = T[i][i];
  
  return answer;
}

/* left Given's rotation - G&vL alg. 5.1.6 p. 203 */
void Matrix3D::row_rot(const int i, const int k, const Vector2D & cs) {
  Matrix3D answer;
  const Matrix3D & A = *this;
  double c, s;
  double tau1, tau2;
  int j;

  answer = A;
  c = cs[0];
  s = cs[1];
  for (j = 0; j < 3; ++j) {
    tau1 = answer[i][j];
    tau2 = answer[k][j];
    answer[i][j] = c*tau1 - s*tau2;
    answer[k][j] = s*tau1 + c*tau2;
  }

  // return answer;
  *this = answer;
}

/* right Given's rotation - G&vL alg. 5.1.7 p. 203 */
void Matrix3D::col_rot(const int j, const int k, const Vector2D & cs) {
  Matrix3D answer;
  const Matrix3D & A = *this;
  double c, s;
  double tau1, tau2;
  int i;

  answer = A;
  c = cs[0];
  s = cs[1];
  for (i = 0; i < 3; ++i) {
    tau1 = answer[i][j];
    tau2 = answer[i][k];
    answer[i][j] = c*tau1 - s*tau2;
    answer[i][k] = s*tau1 + c*tau2;
  }

  // return answer;
  *this = answer;
}

/* Householder tridiagonalization
   Golub & Van Loan alg. 8.2.1 p. 420 
   accumulate the orthogonal matrix in Q */
/* returns tridiagonal matrix T = transpose(newQ)*A*newQ */
/* updates Q with newQ * Q */
/* Q should be initialized to the identity matrix, or some matrix that
   you are accumulating, before calling this function */
Matrix3D Matrix3D::sym_tridiag(Matrix3D * Q) const {
  const Matrix3D & A = *this;
  Matrix3D answer;
  MatrixMN A2;
  VectorND v, w;
  double beta;
  int r, k, i;

  answer = A;
  r = 3;

  for (k = 0; k < (r - 2); ++k) {
    A2 = answer.get_submatrix(k+1,r-1,k,k).transpose();
    v = A2[0];
    
    v = v.householder();
    beta = -2.0 / v.dot(v);

    w = (answer.get_submatrix(k,r-1,k+1,r-1) * v) * beta;
    answer.set_submatrix(answer.get_submatrix(k+1,r-1,k,r-1) + v.outer_product(w),
			 k+1,k);
    
    w = (answer.get_submatrix(k,r-1,k+1,r-1) * v) * beta;
    answer.set_submatrix(
			 answer.get_submatrix(k,r-1,k+1,r-1) + w.outer_product(v),
			 k,k+1);
    
    /* put real zero's into the now small elements */
    if (1) {
      for (i = k + 2; i < r; ++i) {
	answer[k][i] = answer[i][k] = 0.0;
      }
    }
    
    if (Q != NULL) {
      /* Update Q, the orthogonal matrix */
      w = (Q->get_submatrix(k,r-1,k+1,r-1) * v) * beta;
      Q->set_submatrix(Q->get_submatrix(k,r-1,k+1,r-1) + w.outer_product(v),k,k+1);
    }
  }

  return answer;
}

// *** Non-member functions ***

Matrix3D zeros3() {
  Matrix3D answer;
  int i,j;
  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      answer[i][j] = 0;
  return answer;
}

Matrix3D eye3() {
  Matrix3D answer;
  int i,j;
  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j) {
      if (i == j) answer[i][j] = 1.0;
      else answer[i][j] = 0.0;
    }
  return answer;
}

Matrix3D axis_angle(const Vector3D & axis, const double angle) {
  Matrix3D R1, R2, R3;
  Vector3D test1, test2;

  /* first make a matrix that rotates v onto the x-axis */
  // First row
  R1[0] = axis.unit();

  // Second row
  test1 = R1[0].cross(eye3()[0]);
  test2 = R1[0].cross(eye3()[1]);
  if (test1.length() > test2.length()) R1[1] = test1.unit();
  else R1[1] = test2.unit();

  // Third row
  R1[2] = R1[0].cross(R1[1]).unit();

  /* then rotate about x */
  R2 = eye3();
  R2[1][1] = cos(angle);
  R2[2][2] = R2[1][1];
  R2[1][2] = sin(angle);
  R2[2][1] = -R2[1][2];

  /* now rotate back */
  R3 = R1.transpose();

  return R3 * R2 * R1;
}
