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
// Revision 1.1  2004/06/04 19:34:48  cmbruns
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
#include "BaseMatrix.h"

BaseMatrix & BaseMatrix::operator*=(double scale) {
  BaseMatrix & answer = *this;
  int i,j;
  for (i = 0; i < get_n_rows(); i ++)
    for (j = 0; j < get_n_columns(); j ++)
      answer[i][j] *= scale;
  return *this;
}

ostream & operator<<(ostream & os, const BaseMatrix & M) {
  // Remember initial format flags settings
  ios::fmtflags old_format = os.flags(); // new iostream way
  // long old_format = os.flags(); // old iostream way

  int i;
  for (i = 0; i < M.get_n_columns(); ++i) {
    os << M[i] << endl;
  }

  // Restore initial settings
  os.flags(old_format);

  return os;
}
