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
// Revision 1.4  2002/09/14 00:02:51  bruns
// Added license header to most .cc files
//
// Revision 1.3  2002/09/13 22:46:43  bruns
// removed redundant parameter default
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#include "Matrix4D.h"
#include "Matrix3D.h"
#include "Vector3D.h"
#include "MatrixMN.h"
#include "Vector2D.h"

#define EPSILON FLT_EPSILON
#define ABS(x) (((x) < 0) ? (-(x)) : (x))

// Make homogeneous matrix from rotation matrix and translation
Matrix4D & Matrix4D::set(const Matrix3D & R, const Vector3D & t) {
  Matrix4D & M = *this;
  M.set_submatrix(R, 0, 0);
  M[0][3] = t[0];
  M[1][3] = t[1];
  M[2][3] = t[2];
  M[3][0] = 0.0;
  M[3][1] = 0.0;
  M[3][2] = 0.0;
  M[3][3] = 1.0;
  return *this;
}

const Vector4D & Matrix4D::operator[](int index) const {
  return private_row[index];
}

Vector4D & Matrix4D::operator[](int index) {
  return private_row[index];
}

Matrix4D Matrix4D::operator*(const Matrix4D & m2) const {
  int i, j, k;
  Matrix4D answer = *this;
  const Matrix4D & m1 = *this;

  for (i = 0; i < 4; ++i)
    for (j = 0; j < 4; ++j) {
      answer[i][j] = 0;
      for (k = 0; k < 4; ++k)
	answer[i][j] += m1[i][k] * m2[k][j];
    }
  return answer;
}

Vector4D Matrix4D::operator*(const Vector4D & v) const {
  int i, j;
  Vector4D answer;
  const Matrix4D & M = *this;

  for (i = 0; i < 4; ++i) {
    answer[i] = 0;
    for (j = 0; j < 4; ++j)
      answer[i] += M[i][j] * v[j];
  }
  return answer;
}

Matrix3D Matrix4D::get_rotation() const {
  const Matrix4D & M = *this;
  Matrix3D rotation(M[0][0], M[0][1], M[0][2],
		    M[1][0], M[1][1], M[1][2],
		    M[2][0], M[2][1], M[2][2]);
  return rotation;
}
Vector3D Matrix4D::get_translation() const {
  const Matrix4D & M = *this;
  Vector3D translation(M[0][3], M[1][3], M[2][3]);
  return translation;
}
void Matrix4D::set_rotation(const Matrix3D & R) {
  set_submatrix(R, 0, 0);
}
void Matrix4D::set_translation(const Vector3D & t) {
  set_submatrix(t.to_column(), 0, 3);
}

Matrix4D Matrix4D::inverse_homog() const { // Inverse of homogeneous 3D transform matrix
  const Matrix4D & M = *this;
  Vector3D translation = M.get_translation();
  Matrix3D rotation = M.get_rotation();

  // Here is the inversion
  rotation = rotation.transpose();
  translation = rotation * translation * -1;

  Matrix4D answer;
  answer = *this;
  answer.set_rotation(rotation);
  answer.set_translation(translation);
  
  return answer;
}

// Homogeneous transformation
Vector3D Matrix4D::operator*(const Vector3D & v) const {
  Vector4D long_answer;
  long_answer.set(v[0], v[1], v[2], 1.0);
  const Matrix4D & M = *this;
  long_answer = M * long_answer;

  Vector3D answer;
  answer[0] = long_answer[0];
  answer[1] = long_answer[1];
  answer[2] = long_answer[2];

  return answer;
}

Matrix4D Matrix4D::operator+(const Matrix4D & m2) const {
  Matrix4D answer = *this;
  int i, j;
  const Matrix4D m1 = *this;

  for (i = 0; i < 4; ++i)
    for (j = 0; j < 4; ++j)
      answer[i][j] += m2[i][j];
  
  return answer;
}

Matrix4D Matrix4D::operator-(const Matrix4D & m2) const {
  Matrix4D answer = *this;
  int i, j;
  // const Matrix4D & m1 = *this;

  for (i = 0; i < 4; ++i)
    for (j = 0; j < 4; ++j)
      answer[i][j] -= m2[i][j];
  
  return answer;
}

Matrix4D Matrix4D::scale(double scale) const {
  Matrix4D answer = *this;
  int i, j;
  const Matrix4D & M = *this;

  for (i = 0; i < 4; ++i)
    for (j = 0; j < 4; ++j)
      answer[i][j] = scale * M[i][j];

  return answer;
}

double Matrix4D::trace() const {
  const Matrix4D & M = *this;
  int i;
  double answer = 0.0;
  for (i = 0; i < 4; ++i)
    answer += M[i][i];
  return answer;
}

Matrix4D Matrix4D::transpose() const {
  Matrix4D answer;
  int i, j;
  const Matrix4D m = *this;

  for (i = 0; i < 4; ++i)
    for (j = 0; j < 4; ++j)
      answer[j][i] = m[i][j];
  
  return answer;
}

Vector4D Matrix4D::top_eigenvector() const {
  const Matrix4D & P = *this;
  Matrix4D eigenvectors;
  Vector4D eigenvalues;
  Vector4D answer;
  int i, j;

  eigenvectors = eye4();
  eigenvalues = P.sym_eigen(&eigenvectors);

  /* find top eigenvalue */
  j = 0; /* guess */
  for (i = 1; i < 4; ++i)
    if (eigenvalues[i] > eigenvalues[j]) j = i;
  eigenvectors = eigenvectors.transpose();
  answer = eigenvectors[j]; /* top eigenvector */
  
  return answer;
}

/* symmetric matrix eigenvector/value decomposition */
  // FIXME
/* eigenvalues and eigenvectors of a symmetric matrix */
/* algorithm 8.2.3 of G&vL (p. 423) - symmetric QR algorithm */
Vector4D Matrix4D::sym_eigen(Matrix4D * Q) const {
  const Matrix4D & A = *this;
  Vector4D answer;
  /* Q is the matrix of eigenvectors */
  /* D is a diagonal matrix of eigenvalues */
  /* epsilon is a measure of the unit roundoff */

  Matrix4D T;
  MatrixMN T22;
  int i, j;
  int p, q, n;
  long double epsilon;
  
  /* 1 - form tridiagonal matrix by G&vL 8.2.1 (p. 420) */
  n = 4;
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


/* Householder tridiagonalization
   Golub & Van Loan alg. 8.2.1 p. 420 
   accumulate the orthogonal matrix in Q */
/* returns tridiagonal matrix T = transpose(newQ)*A*newQ */
/* updates Q with newQ * Q */
/* Q should be initialized to the identity matrix, or some matrix that
   you are accumulating, before calling this function */
Matrix4D Matrix4D::sym_tridiag(Matrix4D * Q) const {
  const Matrix4D & A = *this;
  Matrix4D answer;
  MatrixMN A2;
  VectorND v, w;
  double beta;
  int r, k, i;

  answer = A;
  r = 4;

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

MatrixMN Matrix4D::get_submatrix(const int m1, const int m2, const int n1, const int n2) const {
	MatrixMN answer(m2 - m1 + 1, n2 - n1 + 1);
	int i;
	// int j;
	const Matrix4D & m = *this;
	
	int n_rows = m2 - m1 + 1;
	
	for (i = 0; i < n_rows; i++ )
		answer[i] = m[i + m1].subvector(n1,n2);
	
	return answer;
}

Matrix4D & Matrix4D::set_submatrix(const BaseMatrix & M2, int m, int n) {
  int i, j;
  Matrix4D & M1 = *this;
  
  for (i = 0; i < M2.get_n_rows(); ++i)
    for (j = 0; j < M2.get_n_columns(); ++j)
      M1[i + m][j + n] = M2[i][j];
  
  return *this;
}

Matrix4D eye4() {
  Matrix4D answer;
  int i, j;
  for (i = 0; i < 4; i ++) 
    for (j = 0; j < 4; j ++) {
      if (i == j) answer[i][j] = 1.0;
      else answer[i][j] = 0.0;
    }
  return answer;
}

Matrix4D zeros4() {
  Matrix4D answer;
  int i, j;
  for (i = 0; i < 4; i ++) 
    for (j = 0; j < 4; j ++)
      answer[i][j] = 0.0;
  return answer;
}

/* left Given's rotation - G&vL alg. 5.1.6 p. 203 */
void Matrix4D::row_rot(const int i, const int k, const Vector2D & cs) {
  Matrix4D answer;
  const Matrix4D & A = *this;
  double c, s;
  double tau1, tau2;
  int j;

  answer = A;
  c = cs[0];
  s = cs[1];
  for (j = 0; j < 4; ++j) {
    tau1 = answer[i][j];
    tau2 = answer[k][j];
    answer[i][j] = c*tau1 - s*tau2;
    answer[k][j] = s*tau1 + c*tau2;
  }

  // return answer;
  *this = answer;
}

/* right Given's rotation - G&vL alg. 5.1.7 p. 203 */
void Matrix4D::col_rot(const int j, const int k, const Vector2D & cs) {
  Matrix4D answer;
  const Matrix4D & A = *this;
  double c, s;
  double tau1, tau2;
  int i;

  answer = A;
  c = cs[0];
  s = cs[1];
  for (i = 0; i < 4; ++i) {
    tau1 = answer[i][j];
    tau2 = answer[i][k];
    answer[i][j] = c*tau1 - s*tau2;
    answer[i][k] = s*tau1 + c*tau2;
  }

  // return answer;
  *this = answer;
}


