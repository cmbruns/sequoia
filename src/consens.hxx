#ifndef __CONSENS_HXX__
#define __CONSENS_HXX__

// $Id: consens.hxx,v 1.2 2001/11/28 23:08:17 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/consens.hxx,v 1.2 2001/11/28 23:08:17 bruns Exp $
// $Log: consens.hxx,v $
// Revision 1.2  2001/11/28 23:08:17  bruns
// Added cvs header tags
// Removed ^M characters
//

#include "sequence.hxx"

Sequence consensus(const Alignment & a, 
		   const Comparison_matrix & cm, 
		   double cutoff);

#endif
