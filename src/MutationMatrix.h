// $Id$
// This file is part of the Sequoia package for macromolecular 
//  sequence/structure analysis
// Copyright (C) 2002  Christopher M. Bruns, Ph.D.
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
// To contact the author, write to cmbruns@attbi.com or bruns@scripps.edu
// In publications please cite: Bruns et al (1999), J.Mol.Biol. 288:427-439
// Please submit bug reports at http://bruns.homeip.net/bugzilla/
// 
// $Id$
// $Header$
// $Log$
// Revision 1.2  2004/05/19 00:58:31  cmbruns
// Add scale_factor argument to constructor
// add prototype for scale routine
//
// Revision 1.1.1.1  2004/05/11 20:26:12  cmbruns
// Initial Max repository for latest sequoia
//
// Revision 1.4  2002/09/13 23:28:08  bruns
// Added license header to all header files
//
// Revision 1.3  2001/12/14 19:21:55  bruns
// Changes to compile with gcc 3.0.2
//  - "using namespace std;"  to use standard classes
//  - "ios::fmtflags" instead of "fmtflags"
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#ifndef _MUTATION_MATRIX_H_
#define _MUTATION_MATRIX_H_

#include <iostream>
#include <vector>
#include <map>
#include <string>

using namespace std;

/* amino acid comparison matrix */
class MutationMatrix {
private:
  vector<string> private_headers;
  vector<char> private_row_order;
  vector<char> private_col_order;
  map<string, double> private_scores;
public:
  MutationMatrix() {}
  MutationMatrix(const char * table, float scale_factor = 1.0);
  void clear();

  double get_score(const char r1, const char r2) const;
  double & set_score(const char r1, const char r2);
  double & set_score(const char r1, const char r2, double value);
  void scale_score(float scale_factor);

  const vector<char> & get_row_order() const {return private_row_order;}
  const vector<char> & get_col_order() const {return private_col_order;}
  vector<char> & set_row_order() {return private_row_order;}
  vector<char> & set_col_order() {return private_col_order;}

  void add_header(const string & line) {private_headers.push_back(line);}
  const vector<string> & get_header() const {return private_headers;}
};

istream & operator>>(istream & is, MutationMatrix & M);
ostream & operator<<(ostream & os, const MutationMatrix & M);

extern MutationMatrix protein_identity;
extern MutationMatrix blosum62;

#endif
