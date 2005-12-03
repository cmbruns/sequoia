#include "ltable.hxx"
#include "geoscore.hxx"
#include "variable.hxx"

// $Id: geoscore.cxx,v 1.2 2001/11/28 23:08:17 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/geoscore.cxx,v 1.2 2001/11/28 23:08:17 bruns Exp $
// $Log: geoscore.cxx,v $
// Revision 1.2  2001/11/28 23:08:17  bruns
// Added cvs header tags
// Removed ^M characters
//

// geoscore.cxx, log-odds scoring for geometric comparison
// of two coordinate files

// from path_mat.cxx, to satisfy old functions
#define MAX_RATIO 100.0

#define MIN_DISTANCE_SIGMA 0.001
// probability ratio score based on distance
// as in paper (Bruns, SEQUOIA)
// distance is in Angstroms
double distance_score2(double distance, double sigma, double sphere_rad)
{
  double Eij, Rij; // num. and den. of log-odds expression
  if (sigma <= MIN_DISTANCE_SIGMA) sigma = MIN_DISTANCE_SIGMA;

  // Numerator of log-odds expression (Eij):
  Eij = ZPROB(distance/sigma);

  // Denominator of log-odds expression (Rij):
  // sigma added to prevent divide-by-zero in return value
  Rij = radprob((distance+sigma)/sphere_rad);

  return Eij/Rij;
}

// comparison score between two residues in PDB files
// DISTANCE INFORMATION ONLY
double distance_score1(double distance, double sigma)
{
  return sigma*sigma/((distance+sigma)*(distance+sigma));
}

#define MIN_ANGLE_SIGMA 0.0001
// log-odds score based on orientation
// as in paper (Bruns, SEQUOIA)
// angle is in radians
double orientation_score2(double angle, double sigma)
{
  double Eij, Rij; // num. and den. of log-odds expression

  if (sigma <= MIN_ANGLE_SIGMA) sigma = MIN_ANGLE_SIGMA;

  // Numerator of log-odds expression (Eij):
  Eij = ZPROB(angle/sigma);

  // Denominator of log-odds expression (Rij):
  // DONT WAN'T TO DIVIDE BY ZERO, SO INCLUDE SIGMA TO MAKE IT POSITIVE
  Rij = ANGPROB(angle) + ANGPROB(sigma);

  double answer = Eij/Rij;
  if (!(answer >= 0)) // NaN
    {
	  cerr << "Internal error: orientation score NaN" << endl;
	  cerr << "answer = " << answer << endl;	
	  cerr << "Eij = " << Eij << endl;
	  cerr << "Rij = " << Rij << endl;
	  cerr << "angle = " << angle << endl;
	  cerr << "sigma = " << sigma << endl;
    }
  return answer;
}


// comparison score between two residues in PDB files
// ORIENTATION INFORMATION ONLY
// this was the original, slow, working way
double orientation_score1(double angle, double sigma)
{
  return sigma*sigma/((angle+sigma)*(angle+sigma));
}


// using COSINE as argument
// comparison score between two residues in PDB files
// ORIENTATION INFORMATION ONLY
// this was the original, slow, working way
double orientation_score1C(double cosangle, double sigma)
{
  return sigma*sigma/((cosangle+sigma)*(cosangle+sigma));
}

// actual fast way to get cosine of kappa angle
double kappa_cos_angle(const RMat & o1, const RMat & o2)
{
  double answer = 0;
  // only need diagonal elements of o1*o2.T();
  double a11 = o1[0][0]*o2[0][0] + o1[0][1]*o2[0][1] + o1[0][2]*o2[0][2];
  double a22 = o1[1][0]*o2[1][0] + o1[1][1]*o2[1][1] + o1[1][2]*o2[1][2];
  double a33 = o1[2][0]*o2[2][0] + o1[2][1]*o2[2][1] + o1[2][2]*o2[2][2];
  answer = (a11 + a22 + a33 - 1.0)/2.0;
  // avoid fatal roundoff errors, while still choking on fundamental errors
  if ((answer < -1.0) && (answer > -1.01)) answer = -1.0;
  else if ((answer > 1.0) && (answer < 1.01)) answer = 1.0;
  return answer;
}

  double answer = 0;
// Try to make cheapest way to get cosine
double fast_cos_angle(const RMat & o1, const RMat & o2)
{
  double answer = 0;
  Vector3D a, b;
  float len2[3];
  Vector3D rotaxis;
  RMat I(3);
  RMat Iprime(3);
  RMat diff(3);
  uint big, littl;

  // 2) Find displacements of the coordinate axes
  I = eye(3);
  Iprime = o1 * o2.T(); // Actual difference rotation matrix
  diff = Iprime - I;

  // find longest and shortest axial changes */
  // Can order them faster using squares
  len2[0] = diff[0] ^ diff[0];
  len2[1] = diff[1] ^ diff[1];
  len2[2] = diff[2] ^ diff[2];
  big = 0; littl = 0;
  if (len2[1] > len2[0]) big = 1;
  else littl = 1;
  if (len2[2] > len2[big]) big = 2;
  if (len2[2] < len2[littl]) littl = 2;

  // Are the two matrices the same?
  double small = 0.0001; // approximately zero
  if ((diff[big] ^ diff[big]) < small) return 0;

  /* find rotation axis using two longest differences */
  if (littl == 0) rotaxis = diff[1] % diff[2];
  else if (littl == 1) rotaxis = diff[2] % diff[0];
  else rotaxis = diff[0] % diff[1];

  /* find rotation angle using the longest difference */
  /* a and b are projections of untransformed and transformed
     axis onto the plane perpendicular to rotation axis */
  a = rotaxis % (I[big] % rotaxis);
  b = rotaxis % (Iprime[big] % rotaxis);

  answer = (a ^ b)/(a ^ a); // cosine of angle between them
  return answer;

}


// probability score that two orientations are equivalent
double p_orient1(const RMat & o1, const RMat & o2)
{
  double angle = acos(kappa_cos_angle(o1, o2)) * 180 / PI;

  double acutoff = sq_var_ptr("ACUTOFF")->value();
  double answer = acutoff*acutoff/((angle+acutoff)*(angle+acutoff));

  return answer;
}

