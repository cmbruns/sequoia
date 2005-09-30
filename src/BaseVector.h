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
// Revision 1.1  2004/06/04 19:34:47  cmbruns
// Imported structure related sources from archive on baxter
// Debugged simple conversion of structures to sequences.
// Implemented computation of solvent accessible surface areas
// Created target residue_area, for output of residue solvent accessible surfaces areas
// Updated GPL headers
//
// Revision 1.5  2002/09/13 23:28:08  bruns
// Added license header to all header files
//
// Revision 1.4  2001/12/14 19:21:55  bruns
// Changes to compile with gcc 3.0.2
//  - "using namespace std;"  to use standard classes
//  - "ios::fmtflags" instead of "fmtflags"
//
// Revision 1.3  2001/11/15 20:52:20  bruns
// included <cfloat> so that FLT_EPSILON will be defined for the numerical routines
//
// Revision 1.2  2001/11/15 20:19:55  bruns
// Added cvs tags up through FastaSequence.h
//
#ifndef _BASE_VECTOR_H_
#define _BASE_VECTOR_H_

using namespace std;

#include <iostream>
// FLT_EPSILON
#include <cfloat>

class MatrixMN;

// Common base class for algebraic Vector types
class BaseVector {
protected:
public:
  virtual const double & operator[](const int index) const = 0;
  virtual double & operator[](const int index) = 0;
  virtual int get_size() const = 0;

  MatrixMN to_row() const;
  MatrixMN to_column() const;
};

ostream & operator<<(ostream & os, const BaseVector & v);

#endif
