#ifndef _SORT_H_
#define _SORT_H_

// $Id: sort.h,v 1.3 2002/06/28 17:17:46 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/sort.h,v 1.3 2002/06/28 17:17:46 bruns Exp $
// $Log: sort.h,v $
// Revision 1.3  2002/06/28 17:17:46  bruns
// Changed to avoid warnings with gcc 3.1
//   changed <iostream.h> type includes to <iostream>
//   removed redundant function parameter default values from .cxx files
//     (when they are already defined in the header)
//
// Revision 1.2  2001/11/29 00:04:43  bruns
// Removed ^M characters
// Added cvs tags
//

#include <iostream> // to define NULL on Linux

template <class Type>
void quicksort(Type * inarray, int size, int * order = NULL);

#endif
