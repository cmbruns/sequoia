#ifndef _TIMELICENSE_HXX_
#define _TIMELICENSE_HXX_

// $Id: timelice.hxx,v 1.5 2002/06/28 17:17:46 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/timelice.hxx,v 1.5 2002/06/28 17:17:46 bruns Exp $
// $Log: timelice.hxx,v $
// Revision 1.5  2002/06/28 17:17:46  bruns
// Changed to avoid warnings with gcc 3.1
//   changed <iostream.h> type includes to <iostream>
//   removed redundant function parameter default values from .cxx files
//     (when they are already defined in the header)
//
// Revision 1.4  2002/02/19 17:26:49  bruns
// Set time license start to Jan, 2000, to avoid some "not yet active" errors
//
// Revision 1.3  2001/11/29 00:04:43  bruns
// Removed ^M characters
// Added cvs tags
//

extern "C" {
#include <sys/time.h>
}
#include<iostream>

// expire on December 1, 2005
#define START_YEAR 2000
#define START_MONTH 1
#define END_YEAR 2005
#define END_MONTH 11

int timed_license();

#endif
