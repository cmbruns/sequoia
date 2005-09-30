/* Copyright (c) 2005 Christopher M. Bruns
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
 
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
