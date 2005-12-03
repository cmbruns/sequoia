#include "rmat.hxx"
#include "cmbmacro.cxx"
#include "matrix.cxx"
#include "rvec.cxx"

// $Id: rmat.cxx,v 1.6 2005/12/03 03:27:48 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/rmat.cxx,v 1.6 2005/12/03 03:27:48 bruns Exp $
// $Log: rmat.cxx,v $
// Revision 1.6  2005/12/03 03:27:48  bruns
// Commit latest files from baxter.
// I don't know how old they are
//
// Revision 1.5  2003/01/04 17:22:41  bruns
// type cast second argument of pow() function, to compile on osf1
//
// Revision 1.4  2002/06/28 17:18:39  bruns
// Added explicit type casts in pow() funtion to avoid compiler warning with gcc 3.1
//
// Revision 1.3  2001/12/12 18:26:27  bruns
// Many changes to get compilation with gcc 3.0.2
// Removed all templated friend functions (no friends needed)
// Added .c_str() method to Mystring, and removed auto conversion to const char *
//
// Revision 1.2  2001/11/28 23:40:22  bruns
// Added cvs header tags
// Removed ^M characters
//

RMat operator*(double r, const RMat & M)
{
  return r * ((Matrix<Real>)M);
}

// return a matrix containing the eigenvectors of the symmetric input
// matrix, plus an extra row containing the corresponding eigenvalues.

RMat RMat::sym_eigen() const
{
  uint r = MIN(m(), n());
  RMat answer (r+1, r);
  RMat Q = eye(r);
  //  cerr << 2 << ", " << 1 - FSIGNIF << endl;
  double epsilon = pow(2.0, (int)(1 - FSIGNIF)); // pretend float accuracy

  // symmetric QR algorithm
  // Golub and Van Loan 8.2.3
  RMat TD = *this;
  TD.sym_tridiag(Q); // reduce to tridiagonal form
  uint i;
  uint q = 0;
  while (q < r)
    {
      for (i=0; i < r-1; ++i)
	if (ABS(TD[i+1][i]) <= epsilon * (TD[i][i] + TD[i+1][i+1]))
	  TD[i+1][i] = TD[i][i+1] = 0;
      // find largest q such that lower block of size q is diagonal
      q = 0;
      while ((q < (r-1)) && (! TD[r-1-q][r-2-q]))
	++q;
      if (q == (r-1)) q = r;
      // find smallest p such that central block is unreduced
      int p;
      if (q == r) p = 0;
      else
	{
	  p = r - q - 2;
	  while ((p > 0) && TD[p-1][p])
	    --p;
	}
      // perform QR step
      if (q < r)
	TD.sym_QRstep(p,r-1-q, Q);
    }
  answer.sub(0,r-1, 0,r-1) = Q.T();
  answer[r] = TD.diag(); // eigenvalues
  return answer;
}


// Householder tridiagonalization
// Golub & Van Loan 8.2.1
// accumulate the orthogonal matrix in Q
// should also have a version without this
void RMat::sym_tridiag(RMat & Q)
{
  uint r = MIN(m(), n());
  uint k, i;
  RVec v, w, tmp;
  RMat & TD = *this;
  Real beta;
  for (k=0; k < r-2; ++k)
    {
// redid algorithm less efficiently, but to actually zero the elements
      tmp = TD.ccol(k).csub(k+1,r-1);
      v = tmp.house();
      beta = -2 / (v^v);
      w = beta * RMat(TD.csub(k,r-1,k+1,r-1)) * v;
      TD.sub(k+1,r-1, k,r-1) = TD.csub(k+1,r-1, k,r-1) + (v * w.T());
      w = beta * RMat(TD.csub(k,r-1,k+1,r-1)) * v;
      TD.sub(k,r-1, k+1,r-1) = TD.csub(k,r-1, k+1,r-1) + (w * v.T());
      // put real zeros into the now small elements
      for (i=k+2; i < r; ++i)
	TD[k][i] = TD[i][k] = 0;
      // Update Q, the orthogonal matrix
      w = beta * RMat(Q.csub(k,r-1,k+1,r-1)) * v;
      Q.sub(k,r-1, k+1,r-1) = Q.csub(k,r-1, k+1,r-1) + (w * v.T());
    }
}


// Householder reduction to upper Hessenberg form
// In the symmetric case, this should reduce the matrix to
// tridiagonal form.  I cannot get the tridiagonalization
// algorithm to work, so I am trying this
// Golub & Van Loan 7.4.2
//void RMat::house_hess(RMat & Q)
//{
//  uint k;
//  RVec v, tmp; // Householder vector
//  for (k=0; k < n()-2; ++k)
//    {
//      tmp = ccol(k).csub(k+1,n()-1);
//      v = tmp.house();
//      sub(k+1,n()-1,k,n()-1) = row_house(csub(k+1,n()-1,k,n()-1), v);
//      sub(0,n()-1,k+1,n()-1) = col_house(csub(0,n()-1,k+1,n()-1), v);
//      // store essential part of v(k) in col(k)(k+2,n()-1)
//      // col(k).sub(k+2,n()-1) = v.sub(1,n()-1);
//    }
//}


// Householder premultiplication
// Golub & Van Loan 5.1.2
RMat row_house(const RMat & A, const RVec & v)
{
  double beta = -2/(v^v);
  RVec w = beta * RMat(A.T()) * v;
  return A + (v * w.T());
}


// Householder premultiplication
// Golub & Van Loan 5.1.3
RMat col_house(const RMat & A, const RVec & v)
{
  double beta = -2/(v^v);
  RVec w = beta * A * v;
  return A + (w * v.T());
}


// Golub and Van Loan 8.2.2
// Implicit Symmetric QR Step with Wilkinson shift
// i and j are the extents of the submatrix
// Q is the accumulated orthogonal matrix
void RMat::sym_QRstep(uint i, uint j, RMat & Q)
{
  RMat & t = *this;
  RMat G;
  Real d = (t[j-1][j-1] - t[j][j])/2;
  Real t2 = t[j][j-1]; t2 *= t2;
  Real mu = t[j][j] - t2/(d + sign(d) * sqrt(d*d + t2));
  Real x = t[i][i] - mu;
  Real z = t[i+1][i];
  uint k;
  Vector2D cs;
  for (k=i; k < j; ++k)
    {
      cs = givens(x,z);
      G = eye(m());
      G[k][k] = G[k+1][k+1] = cs.x();
      G[k][k+1] = cs.y();
      G[k+1][k] = -cs.y();
      t.row_rot(k, k+1, cs);
      t.col_rot(k, k+1, cs);
      Q.col_rot(k, k+1, cs);
      if (k < j-1)
	{
	  x = t[k+1][k];
	  z = t[k+2][k];
	}
    }
}

// Golub and Van Loan 8.2.2
// Implicit Symmetric QR Step with Wilkinson shift
// i and j are the extents of the submatrix
// void RMat::sym_QRstep(uint i, uint j)
// {
//   RMat & t = *this;
//   RMat G;
//   Real d = (t[j-1][j-1] - t[j][j])/2;
//   Real t2 = t[j][j-1]; t2 *= t2;
//   Real mu = t[j][j] - t2/(d + sign(d) * sqrt(d*d + t2));
//   Real x = t[i][i] - mu;
//   Real z = t[i+1][i];
//   uint k;
//   Vector2D cs;
//   for (k=i; k < j; ++k)
//     {
//       cs = givens(x,z);
//       t.row_rot(k, k+1, cs);
//       t.col_rot(k, k+1, cs);
//       if (k < j-1)
//      {
//        x = t[k+1][k];
//        z = t[k+2][k];
//      }
//     }
// }

// Real operator*(const row_vector<Real> & v1, const Vector<Real> & v2)
// {return RVec(v1.T()) ^ RVec(v2);}

// calculate c and s for givens rotation
// Golub and Van Loan 5.1.5
Vector2D givens(Real a, Real b)
{
  Vector2D cs;
  Real tau, c, s;
  if (b == 0) {c = 1; s = 0;}
  else
    {
      if (ABS(b) > ABS(a))
	{
	  tau = -a/b;
	  s = 1/sqrt(1 + tau*tau);
	  c = s * tau;
	}
      else
	{
	  tau = -b/a;
	  c = 1/sqrt(1 + tau*tau);
	  s = c * tau;
	}
    }
  cs.x() = c;
  cs.y() = s;
  return cs;
}


// apply left Givens rotation
// Golub and Van Loan 5.1.6
void RMat::row_rot(const uint i, const uint k, const Vector2D cs)
{
  RMat & A = *this;
  Real c = cs.x();
  Real s = cs.y();
  Real tau1, tau2;
  uint q = n();
  uint j;
  for (j=0; j < q; ++j)
    {
      tau1 = A[i][j];
      tau2 = A[k][j];
      A[i][j] = c*tau1 - s*tau2;
      A[k][j] = s*tau1 + c*tau2;
    }
}

// apply right Givens rotation
// Golub and Van Loan 5.1.7
void RMat::col_rot(const uint j, const uint k, const Vector2D cs)
{
  RMat & A = *this;
  Real c = cs.x();
  Real s = cs.y();
  Real tau1, tau2;
  uint q = m();
  uint i;
  for (i=0; i < q; ++i)
    {
      tau1 = A[i][j];
      tau2 = A[i][k];
      A[i][j] = c*tau1 - s*tau2;
      A[i][k] = s*tau1 + c*tau2;
    }
}


// Is the matrix orthogonal?
bool RMat::is_orthogonal() const
{
  const RMat & M = *this;
  uint space = M.dim(); // number of dimensions
  double small = 0.0001; // cutoff for zero value

  uint i, j;
  for (i = 0; i < space - 1; ++i)
    for (j = i + 1; j < space; ++j)
      {
	// dot products should be small, if orthogonal
	if ((ABS(M[i] ^ M[j])) > small) return false;
      }
  return true;
}


// Is the matrix normal?
bool RMat::is_normal() const
{
  const RMat & M = *this;
  uint space = M.dim(); // number of dimensions
  double small = 0.0001; // cutoff for zero value

  uint i;
  for (i = 0; i < space; ++i)
    {
      // dot products should be one, if normal
      if ((ABS((M[i] ^ M[i]) - 1)) > small) return false;
    }
  return true;
}


// Is the matrix orthonormal?
bool RMat::is_orthonormal() const
{
  const RMat & M = *this;
  if ((M.is_orthogonal()) && (M.is_normal())) return true;
  else return false;
}


RMat RMat::strip_trans() const // remove translational component of a homogeneous matrix
{
  const RMat & M = *this;

  if (M.dim() < 1) exit(1);

  uint space = M.dim() - 1;
  RMat answer(space, space);
  uint i, j;
  for (i = 0; i < space; ++i)
    for (j = 0; j < space; ++j)
      answer[i][j] = M[i][j];
  return answer;
}

// attempt to port from vector.c routines assumes that matrix is
// not homogeneous (i.e. does not includes an extra column for
// translation)
// Should return a value between zero and PI; -100 on error
double RMat::angle(const RMat & m2) const // angle between two rotation matrices
{
  const RMat & m1 = *this;
  double answer = 0;
  Vector3D a, b;
  float len[3];
  Vector3D rotaxis;
  RMat I(3);
  RMat Iprime(3);
  RMat diff(3);
  uint big, littl;

  // 1) Is it orthonormal?  If not, return an error
  if ((!m1.is_orthonormal()) || (!m2.is_orthonormal()))
    {
      cerr << "Rotation matrices are not orthonormal!!\n";
      return -100;
    }

  // 2) Find displacements of the coordinate axes
  I = eye(3);
  Iprime = m1 * m2.T(); // Actual difference rotation matrix
  diff = Iprime - I;
  len[0] = length(diff[0]);
  len[1] = length(diff[1]);
  len[2] = length(diff[2]);

  /* find longest and shortest axial changes */
  big = 0; littl = 0;
  if (len[1] > len[0]) big = 1;
  else littl = 1;
  if (len[2] > len[big]) big = 2;
  if (len[2] < len[littl]) littl = 2;

  // Are the two matrices the same?
  double small = 0.0001; // approximately zero
  if (length(diff[big]) < small) return 0;

  /* find rotation axis using two longest differences */
  if (littl == 0) rotaxis = diff[1] % diff[2];
  else if (littl == 1) rotaxis = diff[2] % diff[0];
  else rotaxis = diff[0] % diff[1];

  rotaxis = rotaxis.unit();
  /* find rotation angle using the longest difference */
  /* a and b are projections of untransformed and transformed
     axis onto the plane perpendicular to rotation axis */
  a = (rotaxis % (I[big] % rotaxis)).unit();
  b = (rotaxis % (Iprime[big] % rotaxis)).unit();

  answer = acos(a ^ b); // angle between them
  return answer;
}

// try to do this faster
double RMat::angle2(const RMat & m2) const // angle between two rotation matrices
{
  const RMat & m1 = *this;
  double answer = 0;
  Vector3D a, b;
  float len[3];
  Vector3D rotaxis;
  RMat I(3);
  RMat Iprime(3);
  RMat diff(3);
  uint big, littl;

  // *** first, for speed, assume they are orthonormal (i.e. don't check)

  // 2) Find displacements of the coordinate axes
  // *** don't take square roots just to check size
  I = eye(3);
  Iprime = m1 * m2.T(); // Actual difference rotation matrix
  diff = Iprime - I;
  len[0] = diff[0]^diff[0];
  len[1] = diff[1]^diff[1];
  len[2] = diff[2]^diff[2];

  /* find longest and shortest axial changes */
  big = 0; littl = 0;
  if (len[1] > len[0]) big = 1;
  else littl = 1;
  if (len[2] > len[big]) big = 2;
  if (len[2] < len[littl]) littl = 2;

  // Are the two matrices the same?
  double small = 0.00001; // approximately zero
  if ((diff[big]^diff[big]) < small) return 0;

  /* find rotation axis using two longest differences */
  if (littl == 0) rotaxis = diff[1] % diff[2];
  else if (littl == 1) rotaxis = diff[2] % diff[0];
  else rotaxis = diff[0] % diff[1];

  // *** don't need this square root : rotaxis = rotaxis.unit();
  /* find rotation angle using the longest difference */
  /* a and b are projections of untransformed and transformed
     axis onto the plane perpendicular to rotation axis */
  a = (rotaxis % (I[big] % rotaxis)).unit();
  b = (rotaxis % (Iprime[big] % rotaxis)).unit();

  answer = acos(a ^ b); // angle between them
  return answer;
}

RMat eye(uint r, uint c)
{
  RMat answer(r,c,0);
  uint i;
  for (i=0; i < MIN(r,c); ++i)
    answer[i][i] = 1;
  return answer;
}

RMat eye(uint r)
{
  return eye(r,r);
}

// RMat operator*(const RVec & v1, const row_vector<Real> & v2)
// {
//   return (Vector<Real>) v1 * v2;
// }

// RVec operator*(const RMat & M, const row_vector<Real> & v)
// {
//   RVec answer(M.m());
//   uint i;
//   for (i=0; i < answer.dim();  ++i)
//     answer[i] = (RVec)M[i] ^ (RVec)v;
//   return answer;
// }

RMat trans_mat(const RVec & v)
{
  RMat answer (v.dim()+1, v.dim()+1, 0);
  uint i;
  for (i=0; i < v.dim(); ++i)
    {
      answer[i][v.dim()] = v[i];
      answer[i][i] = 1;
    }
  answer[v.dim()][v.dim()] = 1;
  return answer;
}

Vector3D transform(const RMat & M, const Vector3D v)
{
  if ((M.m() == 3) && (M.n() == 3)) return (M * v);
  else if ((M.m() == 1) && (M.n() >= 3)) return M[0] + v;
  else if ((M.m() == 1) && (M.n() >= 3)) return M[0] + v;
  else if ((M.m() >= 3) && (M.n() == 1)) return M.ccol(0) + v;
  // Homogeneous coordinates - ignore 4th row for now - cmb
  else if ((M.m() >= 3) && (M.n() >= 4))
    return ((M.csub(0,2,0,2) * v) + M.col(3).csub(0,2));
  else return v;
}

Vector3D vec3_transform(const RMat & M, const Vector3D v)
{return transform(M,v);}

double RMat::determinant() const
{
  const RMat & t = *this;
  // must be symmetric ?
  if ( m() != n() ) return 0;
  else if (n() == 0) return 0;
  else if (n() == 1) return t[0][0];
  else
    {
      double answer = 0;
      unsigned int j;
      for (j = 0; j < n(); ++j)
	{
	  // delete first row and jth column of matrix
	  RMat A1j(n()-1, n()-1);
	  unsigned int i2, j2;
	  for (i2 = 0; i2 < n()-1; ++i2)
	    for (j2 = 0; j2 < n()-1; ++j2)
	      {
		if (j2 < j) A1j[i2][j2] = t[i2+1][j2];
		else        A1j[i2][j2] = t[i2+1][j2+1];
	      }
	  answer += pow((double) -1.0, (int) j) * t[0][j] * A1j.determinant();
	}
      return answer;
    }
}
