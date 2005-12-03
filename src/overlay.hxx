#ifndef _OVERLAY_HXX_
#define _OVERLAY_HXX_

// $Id: overlay.hxx,v 1.3 2002/02/19 18:39:13 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/overlay.hxx,v 1.3 2002/02/19 18:39:13 bruns Exp $
// $Log: overlay.hxx,v $
// Revision 1.3  2002/02/19 18:39:13  bruns
// Added ill_conditioned_rotation_error exception class for overlay functions
//
// Revision 1.2  2001/11/28 23:40:22  bruns
// Added cvs header tags
// Removed ^M characters
//

#include "rvec.hxx"
#include "rmat.hxx"
#include "sequence.hxx"
#include "brookhav.hxx"

class ill_conditioned_rotation_error {
public:
  ill_conditioned_rotation_error(void) {
    // cerr << "Ill conditioned rotation" << endl;
  }
  ill_conditioned_rotation_error(const SeqRes & r) {
    // cerr << "Ill conditioned determinant for rotation matrix!" << endl;
    // cerr << "Determinant = " << r.orientation().determinant() << endl;
    // cerr << "Residue = " << r.letter() << r.res_num() << endl;
    // cerr << r.orientation() << endl;
  }
};

// Use Kabsch method
RMat overlay(const Array<RVec> & X, 
	     const Array<RVec> & Y, 
	     const Array<Real> & w) throw(ill_conditioned_rotation_error);
// find transformation
RMat overlay(const pdbprotein & p1, 
		  const pdbprotein & p2, 
		  const Alignment & a);
RMat overlay(const pdboverlay & o1, 
		  const pdboverlay & o2, 
		  const Alignment & a);
RMat overlay(const SeqRes & r1, const SeqRes & r2) throw(ill_conditioned_rotation_error);
// reassign capital letters
Alignment struct_equiv(const pdbprotein & p1, 
		       const pdbprotein & p2,
		       const Alignment & a,
		       Real cutoff, uint runlength);
// WITHIN an overlay
Alignment struct_equiv(const pdboverlay & o, 
		       const Alignment & a,
		       Real cutoff, uint runlength);
// BETWEEN overlays
Alignment struct_equiv(const pdboverlay & o1, 
		       const pdboverlay & o2,
		       const Alignment & a,
		       Real cutoff, uint runlength);
Alignment struct_equiv(const Alignment & a,
		       Real dcutoff, Real acutoff, uint runlength);

#endif
