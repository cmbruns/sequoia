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
