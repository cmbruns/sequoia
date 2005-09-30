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
// Revision 1.3  2001/12/14 19:21:55  bruns
// Changes to compile with gcc 3.0.2
//  - "using namespace std;"  to use standard classes
//  - "ios::fmtflags" instead of "fmtflags"
//
// Revision 1.2  2001/11/15 20:19:55  bruns
// Added cvs tags up through FastaSequence.h
//
#include "BaseVector.h"
#include "MatrixMN.h"

MatrixMN BaseVector::to_row() const {
  MatrixMN answer(1, get_size());
  const BaseVector & v = *this;
  int i;
  for (i = 0; i < get_size(); i++)
    answer[0][i] = v[i];
    return answer;
}
MatrixMN BaseVector::to_column() const {
  MatrixMN answer(get_size(), 1);
  const BaseVector & v = *this;
  int i;
  for (i = 0; i < get_size(); i++)
    answer[i][0] = v[i];
  return answer;
}
	
ostream & operator<<(ostream & os, const BaseVector & v) {
  // Remember initial format flags settings
  ios::fmtflags old_format = os.flags(); // new iostream way
  // long old_format = os.flags(); // old iostream way

  int i;
  os << "(";
  for (i = 0; i < v.get_size(); ++i) {
    os.setf(ios::fixed); // Makes precision be what I want?
    os.precision(5); os.width(9);
    os << v[i];
    if (i < (v.get_size() - 1)) os << "   ";
  }
  os << ")";

  // Restore initial settings
  os.flags(old_format);

  return os;
}
