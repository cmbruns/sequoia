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
// Revision 1.1  2004/06/04 19:34:48  cmbruns
// Imported structure related sources from archive on baxter
// Debugged simple conversion of structures to sequences.
// Implemented computation of solvent accessible surface areas
// Created target residue_area, for output of residue solvent accessible surfaces areas
// Updated GPL headers
//
// Revision 1.6  2002/09/14 00:02:51  bruns
// Added license header to most .cc files
//
// Revision 1.5  2002/09/13 22:40:49  bruns
// Added constructor from Matrix4D data type
//
// Revision 1.4  2002/05/22 17:19:18  bruns
// Added comment about complexity of computing orientation difference angles
//
// Revision 1.3  2001/12/14 19:21:55  bruns
// Changes to compile with gcc 3.0.2
//  - "using namespace std;"  to use standard classes
//  - "ios::fmtflags" instead of "fmtflags"
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#include "Quaternion.h"
#include <cmath>

using namespace std;

#define EPSILON FLT_EPSILON
#define ABS(x) (((x) < 0) ? (-(x)) : (x))

Quaternion::Quaternion() {
  // Set to identity quaternion
  private_l = 0.0;
  private_m = 0.0;
  private_n = 0.0;
  private_s = 1.0;
}

Quaternion::Quaternion(double l, double m, double n, double s) {
  private_l = l;
  private_m = m;
  private_n = n;
  private_s = s;
}

Quaternion::Quaternion(const Matrix3D & M) {
  *this = M.quaternion();
}

Quaternion::Quaternion(const Matrix4D & M) {
  *this = M.get_rotation().quaternion();
}

// Quaternion multiplication
// 13 multiplications, 12 additions/subtractions
// Compare to 27 mulitiplications, 27 additions/subtractions for Matrix3D multiply
// (You could get just the "s" portion of the result with just 4 multiplications
// and 3 addition/subtractions)
Quaternion Quaternion::operator*(const Quaternion & q2) const { // multiplication
  Quaternion answer;
  const Quaternion & q1 = *this;

  /* s3 = s1s2 - v1.v2 */
  answer.set_s() = q1.get_s() * q2.get_s() - 
    (q1.get_l() * q2.get_l() + q1.get_m() * q2.get_m() + q1.get_n() * q2.get_n());
  /* v3 = s1v2 + s2v1 + v1xv2 */
  answer.set_l() = q1.get_s() * q2.get_l() + q2.get_s() * q1.get_l() + 
    (q1.get_m() * q2.get_n() - q1.get_n() * q2.get_m());
  answer.set_m() = q1.get_s() * q2.get_m() + q2.get_s() * q1.get_m() + 
    (q1.get_n() * q2.get_l() - q1.get_l() * q2.get_n());
  answer.set_n() = q1.get_s() * q2.get_n() + q2.get_s() * q1.get_n() + 
    (q1.get_l() * q2.get_m() - q1.get_m() * q2.get_l());

  return answer;
}

/* assumes a unit quaternion */
Matrix3D Quaternion::rotmat() const {
  const Quaternion & q = *this;

  Matrix3D M(
	     1.0 - 2.0*(q.get_m()*q.get_m() + q.get_n()*q.get_n()), 
	     2.0*(q.get_l()*q.get_m() - q.get_s()*q.get_n()), 
	     2.0*(q.get_l()*q.get_n() + q.get_s()*q.get_m()),
	     2.0*(q.get_l()*q.get_m() + q.get_s()*q.get_n()),
	     1.0 - 2.0*(q.get_l()*q.get_l() + q.get_n()*q.get_n()), 
	     2.0*(q.get_m()*q.get_n() - q.get_s()*q.get_l()),
	     2.0*(q.get_l()*q.get_n() - q.get_s()*q.get_m()),
	     2.0*(q.get_m()*q.get_n() + q.get_s()*q.get_l()), 
	     1.0 - 2.0*(q.get_l()*q.get_l() + q.get_m()*q.get_m())
	     );
	     
  return M;
}

Quaternion Quaternion::unit() const {
  Quaternion answer;
  double scale;
  const Quaternion & q1 = *this;

  scale = q1.get_s()*q1.get_s() + 
    q1.get_l()*q1.get_l() + 
    q1.get_m()*q1.get_m() + 
    q1.get_n()*q1.get_n();
  if (ABS(scale - 1.0) < EPSILON) return q1;

  answer = q1;
  scale = 1.0 / sqrt(scale);
  answer.set_s() *= scale;
  answer.set_l() *= scale;
  answer.set_m() *= scale;
  answer.set_n() *= scale;

  return answer;
}

bool Quaternion::is_unit() const {
  double scale;
  const Quaternion & q1 = *this;
  scale = q1.get_s()*q1.get_s() + 
    q1.get_l()*q1.get_l() + 
    q1.get_m()*q1.get_m() + 
    q1.get_n()*q1.get_n();
  if (ABS(scale - 1.0) < EPSILON) return true;
  else return false;
}

ostream & operator<<(ostream & os, const Quaternion & q) {
  // Remember initial format flags settings
  ios::fmtflags old_format = os.flags(); // new iostream way
  // long old_format = os.flags(); // old iostream way

  // int i;
  os << "(";

  os.setf(ios::fixed); // Makes precision be what I want?
  os.precision(5); os.width(9);
  os << q.get_l();

  os << "   ";

  os.setf(ios::fixed); // Makes precision be what I want?
  os.precision(5); os.width(9);
  os << q.get_m();

  os << "   ";

  os.setf(ios::fixed); // Makes precision be what I want?
  os.precision(5); os.width(9);
  os << q.get_n();

  os << "   ";

  os.setf(ios::fixed); // Makes precision be what I want?
  os.precision(5); os.width(9);
  os << q.get_s();

  os << ")";

  // Restore initial settings
  os.flags(old_format);

  return os;
}


