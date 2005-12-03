#ifndef _MYSTRING_HXX_
#define _MYSTRING_HXX_

// $Id: mystring.hxx,v 1.5 2002/06/28 17:21:54 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/mystring.hxx,v 1.5 2002/06/28 17:21:54 bruns Exp $
// $Log: mystring.hxx,v $
// Revision 1.5  2002/06/28 17:21:54  bruns
// Minor changes to avoid compiler warnings with gcc 3.1
//   changed <iostream.h> style includes to modern <iostream> style
//   removed redundant default function parameter values when already specified in
//    the header
//
// Revision 1.4  2001/12/12 18:26:27  bruns
// Many changes to get compilation with gcc 3.0.2
// Removed all templated friend functions (no friends needed)
// Added .c_str() method to Mystring, and removed auto conversion to const char *
//
// Revision 1.3  2001/12/04 23:41:06  bruns
// Removed non-const char * cast, to avoid compiler warnings
//
// Revision 1.2  2001/11/28 23:40:22  bruns
// Added cvs header tags
// Removed ^M characters
//

#include <iostream>

// extern "C" {
#include <ctype.h>
// };

// extern "C" {
#include <string.h>
// };
#include "array.hxx"

class Mystring : public Array<char>
{
private:
public:
  Mystring();
  // operator char *() {return array_ptr();} // dangerous!
  // operator const char *() const {return array_ptr();}
  Mystring(const char c);
  Mystring(const Mystring & s);
  Mystring(const char * s);
  ~Mystring();

  const char * c_str() const {return array_ptr();}

  Mystring & operator=(const Mystring & s);
  Mystring & operator=(const char * s);
  Abstract_array<char> & operator+=(const char & c);
  uint length() const {if (dim() < 2) return 0; else return (dim() - 1);}
  bool operator==(const Mystring & s2) const;
  bool operator==(const char * c) const;
  bool operator!=(const Mystring & s2) const;
  bool operator!=(const char * c) const;
  Mystring upcase() const;
  Mystring downcase() const;
  Mystring at(uint start, uint len) const;
  Mystring after(const uint & start) const;
  Mystring after(const char & c) const;
  Mystring & remove(const char c);

  friend ostream & operator<<(ostream & os, const Mystring & s);
  friend istream & operator>>(istream & is, Mystring & s);
  Array<Mystring> Mystring::split(char delim = ' ') const;
};

Mystring upcase(const Mystring & s);

#endif
