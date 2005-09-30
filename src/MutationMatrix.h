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
// $Header$
// $Log$
// Revision 1.3  2004/06/04 19:10:36  cmbruns
// Updated GPL header
//
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
