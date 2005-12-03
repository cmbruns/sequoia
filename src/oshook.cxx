#include "oshook.hxx"

// $Id: oshook.cxx,v 1.2 2001/11/28 23:40:22 bruns Exp $
// $Log: oshook.cxx,v $
// Revision 1.2  2001/11/28 23:40:22  bruns
// Added cvs header tags
// Removed ^M characters
//
// $Header: /usr/data/cvs/sequoia1/oshook.cxx,v 1.2 2001/11/28 23:40:22 bruns Exp $

#ifdef USE_BUFFERED_OUTPUT

// must use statically allocated buffer to rewind and reuse (I think)
char * osbuf = new char[OSBUFSIZE];
ostrstream status_os(osbuf, OSBUFSIZE);

// must use as "myoshook(os)", for "os << myoshook" dumps core
ostrstream & myoshook(ostrstream & os)
{
  // Add a null character to ensure the character array is terminated
  os << '\0';

  // output character array to standard output stream
  char * cp = os.str();
  cout << cp;

  // rewind character array to beginning, to accept further input
  os.seekp(0, ios::beg);
  return os;
}

#else

ostream & status_os = cout;
ostream & myoshook(ostream & os) {
  return os;
}

#endif
