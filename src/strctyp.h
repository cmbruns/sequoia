#ifndef STRCTYP_H
#define STRCTYP_H

/* $Id: strctyp.h,v 1.2 2001/11/29 00:04:43 bruns Exp $ */
/* $Header: /usr/data/cvs/sequoia1/strctyp.h,v 1.2 2001/11/29 00:04:43 bruns Exp $ */
/* $Log: strctyp.h,v $
/* Revision 1.2  2001/11/29 00:04:43  bruns
/* Removed ^M characters
/* Added cvs tags
/* */

/* data structure info for passing C++ coordinates to C motif windows */

#include "cvector.h"
typedef vector *coordarray;

typedef char reslabel[20];
typedef reslabel * reslabelp;

#endif
