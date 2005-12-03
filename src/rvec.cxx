#include "rvec.hxx"

// $Id: rvec.cxx,v 1.2 2001/11/28 23:40:22 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/rvec.cxx,v 1.2 2001/11/28 23:40:22 bruns Exp $
// $Log: rvec.cxx,v $
// Revision 1.2  2001/11/28 23:40:22  bruns
// Added cvs header tags
// Removed ^M characters
//

Real length(const Vector<Real> & v)
{return sqrt(v ^ v);}

Real operator^(const Vector<Real> & v, const Vector<Real> & v2)
{
  Real dotprod = 0;
  uint r = MIN(v.dim(), v2.dim());
  uint i;
  for (i = 0; i < r; ++i)
    dotprod += (v[i] * v2[i]);
  return dotprod; 
}

Vector3D operator%(const Vector3D & v1, const Vector3D & v2) // cross product
{
  Vector3D a;
  a.x() = v1[1]*v2[2] - v1[2]*v2[1];
  a.y() = v1[2]*v2[0] - v1[0]*v2[2];
  a.z() = v1[0]*v2[1] - v1[1]*v2[0];
  return a;
}

RVec RVec::unit() const
{
  double l = length(*this);
  if (l == 0) error("Cannot normalize zero length vector");
  return (1/l) * (*this);
}


// RVec RVec::operator-() const
// {
//   RVec answer = (RVec) (- (Vector<Real>)(*this));
//   return answer;
// }


Real RVec::pnorm(double p) const
{
  const RVec & x = *this;
  Real answer = 0;
  uint i;
  if (p == 2)
    {
      answer = sqrt(x ^ x);
    } 
  else
    {
      for (i=0; i < dim(); ++i)
	{
	  answer += pow(ABS(x[i]), p);
	}
      answer = pow(answer, 1/p);
    }
  return answer;
}


// compute Householder vector
// Golub and Van Loan 5.1.1
RVec RVec::house() const
{
  const RVec & x = *this;
  Real mu = pnorm(2);
  RVec v = *this;
  Real beta;
  if (mu)
    {
      beta = x[0] + sign(x[0]) * mu;
      v *= 1/beta;
//      v.sub(1,dim()-1) = v.sub(1,dim()-1) * (1/ beta);
    }
  v[0] = 1;
  return v;
}

RVec operator*(Real r, const RVec & v)
{
   return v*r;
}

ostream & operator<<(ostream & os, const RVec & v)
{
  uint i;
  for (i=0; i<v.dim(); ++i)
    {
      os.precision(3);
      os.width(7);
      os.setf(ios::fixed);
      os << v[i];
      os << "\n";
    }
  return os;
}

// Vector3D::Vector3D(const Vector<Real> v) : RVec(3, 0.0)
// {
//   Vector3D & t = *this;
//   t.x() = v[0];
//   t.y() = v[1];
//   t.z() = v[2];
// }

Vector3D::Vector3D(Vector<Real> v) : RVec(3, 0.0)
{
  Vector3D & t = *this;
  t.x() = v[0];
  t.y() = v[1];
  t.z() = v[2];
}

Vector3D & Vector3D::operator=(const Vector3D & a2) 
{
  if (this == &a2) return *this;
  Vector3D & v = *this;
  init_array(3);
  v[0] = a2[0];
  v[1] = a2[1];
  v[2] = a2[2];
  return *this; 
}


Real Vector3D::distance(const Vector3D & v2) const 
{
  Real t = x() - v2.x();
  Real d = (t * t);
  t = y() - v2.y();
  d += (t * t);
  t = z() - v2.z();
  d += (t * t);
  return sqrt(d); 
}

Vector3D & Vector3D::operator%=(const Vector3D & v)
{
  *this = *this % v;
  return *this;
}

Vector2D Vec(Real x, Real y)
{
  Vector2D answer;
  answer.x() = x;
  answer.y() = y;
  return answer;
}

Vector3D Vec(Real x, Real y, Real z)
{
  Vector3D answer;
  answer.x() = x;
  answer.y() = y;
  answer.z() = z;
  return answer;
}


