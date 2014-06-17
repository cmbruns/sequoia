#ifndef __READLINE_H_CMB__
#define __READLINE_H_CMB__

/* $Id: readline.h,v 1.4 2001/12/04 23:46:53 bruns Exp $ */
/* $Header: /usr/data/cvs/sequoia1/readline.h,v 1.4 2001/12/04 23:46:53 bruns Exp $ */
/* 
 * $Log: readline.h,v $
 * Revision 1.4  2001/12/04 23:46:53  bruns
 * Another crack at preventing compiler warnings by changing the comment
 * structure of the cvs log header
 *
 * Revision 1.3  2001/12/04 23:42:13  bruns
 * Removed comment characters inside of header comments, to avoid compiler warnings
 * 
 * Revision 1.2  2001/11/28 23:40:22  bruns
 * Added cvs header tags
 * Removed ^M characters
 * */

char * rl_gets (const char * prompt);

#endif
