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
// Revision 1.4  2004/06/04 19:07:26  cmbruns
// Updated GPL header
//
// Merged in class based exception code from archive on baxter, replacing integer exceptions.
//
// Revision 1.3  2002/09/13 23:28:08  bruns
// Added license header to all header files
//
//

#ifndef __EXCEPTIONS_H__
#define __EXCEPTIONS_H__

#include <iostream>
#include <string>

using namespace std;

// generic Exception
class Exception {
 private:
   string private_message;
 public:
	Exception() {
		private_message = "no message provided";
	}
	Exception(const char * message) {
    private_message = message;
  }
  void report() {
    cerr << private_message << endl;
  }
  string get_message() {
    return private_message;
  }
};

class no_such_atom_error {
 private:
  string private_atom_name;
 public:
  no_such_atom_error(string atom_name) {
    private_atom_name = atom_name;
  }
  string & get_atom_name() {
    return private_atom_name;
  }
};

class AREA_NOT_YET_COMPUTED_EXCEPTION : public Exception {};
class no_such_torsion_error : public Exception {};
class GAP_MODEL_PROBLEM_EXCEPTION : public Exception {};
class NO_SUCH_FILE_EXCEPTION : public Exception {};
class ALIGNMENT_LENGTH_MISMATCH_EXCEPTION : public Exception {};
class STRANGE_ELEMENT_NAME_EXCEPTION : public Exception {};

#endif
