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
// $Header$
// $Log$
// Revision 1.5  2004/06/04 19:15:07  cmbruns
// Updated GPL header
//
// Moved conservidue pair score function from here to be a member function of conservidue
//

#include "Score.h"
// #include "MutationMatrix.h"
// #include "SequenceAlignment.h"

AlignmentScore operator*(const double r, const AlignmentScore & score) {
	AlignmentScore answer = score;
	answer.match() *= r;
	answer.extension() *= r;
	answer.opening() *= r;
	answer.closing() *= r;
	answer.deletion() *= r;
	answer.transition() *= r;
	return answer;
}

