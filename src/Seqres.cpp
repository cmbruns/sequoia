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
// Revision 1.1  2004/06/04 19:34:46  cmbruns
// Imported structure related sources from archive on baxter
// Debugged simple conversion of structures to sequences.
// Implemented computation of solvent accessible surface areas
// Created target residue_area, for output of residue solvent accessible surfaces areas
// Updated GPL headers
//
// Revision 1.4  2002/09/14 00:02:51  bruns
// Added license header to most .cc files
//
// Revision 1.3  2002/09/13 22:30:18  bruns
// removed redundant parameter default
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#include "Seqres.h"
#include <cstdio>

Seqres::Seqres(const char & chain_id) {
  private_chain_id = chain_id;
}

const vector<string> & Seqres::get_sequence() const {return private_sequence;}

const char Seqres::get_chain_id() const {
  return private_chain_id;
}

const int Seqres::get_num_residues() const {
  return private_num_residues;
}

void Seqres::add_residue(const string & residue) {
  private_sequence.push_back(residue);

  // Add residue to counter for residue types
  if (private_seqres_types.find(residue) == private_seqres_types.end())
    private_seqres_types[residue] = 0;
  private_seqres_types[residue] ++;
}

// How many times does this residue occur in the sequence?
const int Seqres::get_residue_count(const string & residue) const {
  if (private_seqres_types.find(residue) == private_seqres_types.end())
    return 0;
  else return (*(private_seqres_types.find(residue))).second;
}

void Seqres::set_length(int num_residues) {
  private_num_residues = num_residues;
}

bool Seqres::error() const {
  // vector<string>::const_iterator i;
  // const vector<string> & s = private_sequence;
  // for (i = s.begin(); i != s.end();  i ++)
  //   cerr << *i << endl;
  if (private_sequence.size() != private_num_residues)
    return true;
  return false;
}
