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
//
// $Header$
//
// $Log$
// Revision 1.1  2004/05/28 19:32:52  cmbruns
// Match and extension components of sum-of-pairs are correct in profile alignment.
// Open, close and Delete components are still lower bounds
// Changed name from cc to cpp to match other files
//
// Revision 1.2  2004/05/19 00:57:25  cmbruns
// Added scale score function
// scale blosum matrix by 1/3 to make it actual bits
//
// Revision 1.1.1.1  2004/05/11 20:26:12  cmbruns
// Initial Max repository for latest sequoia
//
// Revision 1.5  2002/09/14 00:02:51  bruns
// Added license header to most .cc files
//
// Revision 1.4  2002/09/13 22:45:23  bruns
// Upgraded strstream to sstream
//
// Revision 1.3  2002/05/15 21:48:28  bruns
// Added lines to convert lower case to upper case for scoring
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#include "MutationMatrix.h"
#include <sstream>

MutationMatrix::MutationMatrix(const char * table, float scale_factor) {
  MutationMatrix & M = *this;
  istringstream is(table);
  is >> M;
  scale_score(scale_factor);
}

// Multiply all scores by a scale factor
void MutationMatrix::scale_score(float scale_factor) {
	map<string, double>::iterator scoreit;
	for (scoreit = private_scores.begin();
		 scoreit != private_scores.end();
		 scoreit ++) {
		scoreit->second *= scale_factor;
	}

}

void MutationMatrix::clear() {
  private_headers.clear();
  private_row_order.clear();
  private_col_order.clear();
    private_scores.clear();
}

double & MutationMatrix::set_score(const char row, const char col, double value) {
  string key = "";
  key += row;
  key += col;
  private_scores[key] = value;
  return private_scores[key];
}

double MutationMatrix::get_score(const char row, const char col) const {
  double answer = 0.0; // default score

  string key = "";
  // convert to upper case
  key += toupper(row);
  key += toupper(col);

  if (private_scores.find(key) != private_scores.end())
    answer = (*(private_scores.find(key))).second;
  return answer;
}

istream & operator>>(istream & is, MutationMatrix & M) {
  if (!is.good()) return is;
  char line_buffer[1001];
  M.clear();

  char this_row_char;
  unsigned int i;
  unsigned int line_pos;
  unsigned int this_column;
  string line;

  // Header comments
  bool in_header = true;
  while (is.good() && in_header) {
    is.getline(line_buffer, 1000);
    line = line_buffer;

    if (line.length() < 1) continue;
    if (
	(line[0] == '#')
	&& in_header) {
      M.add_header(line);
    }
    else in_header = false;
  }

  // Index row
  /* current line should be column index */
  for (i = 0; i < line.length(); ++i) {
    this_row_char = line[i];
    if (isgraph(this_row_char)) {
      M.set_col_order().push_back(this_row_char);
    }
  }

  /* read in the data lines */
  while (is.good()) {
    is.getline(line_buffer, 1000);
    line = line_buffer;
    if (line.length() < 1) continue;

    /* first find the row character */
    line_pos = 0;
    line_pos = line.find_first_not_of(" \t", line_pos);
    this_row_char = line[line_pos];
    M.set_row_order().push_back(this_row_char);

    line_pos = line.find_first_of(" \t", line_pos);
    line_pos = line.find_first_of("0123456789-+.", line_pos);
    this_column = 0;
    while (line_pos != string::npos) {
    /* read in the values */
      char this_col_char = M.get_col_order()[this_column];
      M.set_score(this_row_char, this_col_char, strtod(line.c_str() + line_pos, NULL));

      this_column ++;
      line_pos = line.find_first_of(" \t", line_pos);
      line_pos = line.find_first_of("0123456789-+.", line_pos);
    }
  }
  return is;
}

ostream & operator<<(ostream & os, const MutationMatrix & M) {

  // Header comments
  const vector<string> & head = M.get_header();
  vector<string>::const_iterator h;
  for (h = head.begin(); h != head.end(); h ++)
    os << *h << endl;

  // Column header
  const vector<char> & col = M.get_col_order();
  vector<char>::const_iterator i, j;
  os << " ";
  for (j = col.begin(); j != col.end(); j++) {
    os << "  ";
    os << *j;
  }
  os << endl;

  // Data lines
  const vector<char> & row = M.get_row_order();
  for (i = row.begin(); i != row.end(); i++) {
    os << *i;
    for (j = col.begin(); j != col.end(); j++) {
      os.width(3);
      os << (int)(M.get_score(*i,*j));
    }
    os << " " << endl;
  }
  return os;
}

MutationMatrix protein_identity = \
"   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  U  X  *\n\
A 10 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  0 -4 \n\
R -4 10 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  0 -4 \n\
N -4 -4 10 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4  0 -4 \n\
D -4 -4 -4 10 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4  0 -4 \n\
C -4 -4 -4 -4 10 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  0 -4 \n\
Q -4 -4 -4 -4 -4 10 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4  0 -4 \n\
E -4 -4 -4 -4 -4 -4 10 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4  0 -4 \n\
G -4 -4 -4 -4 -4 -4 -4 10 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  0 -4 \n\
H -4 -4 -4 -4 -4 -4 -4 -4 10 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  0 -4 \n\
I -4 -4 -4 -4 -4 -4 -4 -4 -4 10 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  0 -4 \n\
L -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 10 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  0 -4 \n\
K -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 10 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  0 -4 \n\
M -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 10 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  0 -4 \n\
F -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 10 -4 -4 -4 -4 -4 -4 -4 -4 -4  0 -4 \n\
P -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 10 -4 -4 -4 -4 -4 -4 -4 -4  0 -4 \n\
S -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 10 -4 -4 -4 -4 -4 -4 -4  0 -4 \n\
T -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 10 -4 -4 -4 -4 -4 -4  0 -4 \n\
W -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 10 -4 -4 -4 -4 -4  0 -4 \n\
Y -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 10 -4 -4 -4 -4  0 -4 \n\
V -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 10 -4 -4 -4  0 -4 \n\
B -4 -4  5  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 10 -4 -4  0 -4 \n\
Z -4 -4 -4 -4 -4  5  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 10 -4  0 -4 \n\
U -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 10  0  5 \n\
X  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 \n\
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5  0 10 \n";

// Entries are 1/3 bits, so multiply by 0.3333 to get bits -> see final argument
MutationMatrix blosum62(\
"   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  U  X  *\n\
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1 -2  0 -4 \n\
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -2  0 -4 \n\
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -2  0 -4 \n\
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -2  0 -4 \n\
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3  2  0 -4 \n\
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -2  0 -4 \n\
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -2  0 -4 \n\
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -2  0 -4 \n\
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -2  0 -4 \n\
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -2  0 -4 \n\
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -2  0 -4 \n\
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -2  0 -4 \n\
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -2  0 -4 \n\
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -2  0 -4 \n\
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2  0 -4 \n\
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0 -2  0 -4 \n\
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1 -2  0 -4 \n\
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2  0 -4 \n\
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -2  0 -4 \n\
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -2  0 -4 \n\
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -2  0 -4 \n\
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -2  0 -4 \n\
U -2 -2 -2 -2  2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2  6  0  2 \n\
X  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -4 \n\
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  2 -4  1 \n",
0.3333);

