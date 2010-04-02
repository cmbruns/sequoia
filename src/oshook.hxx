#ifndef OSHOOK1_HXX
#define OSHOOK1_HXX

// $Id: oshook.hxx,v 1.3 2002/06/28 17:21:54 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/oshook.hxx,v 1.3 2002/06/28 17:21:54 bruns Exp $
// $Log: oshook.hxx,v $
// Revision 1.3  2002/06/28 17:21:54  bruns
// Minor changes to avoid compiler warnings with gcc 3.1
//   changed <iostream.h> style includes to modern <iostream> style
//   removed redundant default function parameter values when already specified in
//    the header
//
// Revision 1.2  2001/11/28 23:40:22  bruns
// Added cvs header tags
// Removed ^M characters
//

#ifdef USE_BUFFERED_OUTPUT

#include <sstream>

// header file to define myoshook for output to motif window
// size of character buffer for output
#define OSBUFSIZE 20000

extern ostrstream status_os;
ostrstream & myoshook(ostrstream & os);

#else

#include <iostream>
using namespace std;

extern ostream & status_os;
ostream & myoshook(ostream & os);

#endif

#endif
