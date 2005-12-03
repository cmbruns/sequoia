// Use gnu READLINE interface for getting command lines
#include <stdio.h> // defines FILE for readline library
#include <readline/readline.h> // readline()

/* $Id: readline.c,v 1.2 2001/11/28 23:40:22 bruns Exp $ */
/* $Header: /usr/data/cvs/sequoia1/readline.c,v 1.2 2001/11/28 23:40:22 bruns Exp $ */
/* $Log: */

#ifndef NO_READLINE_HISTORY
#include <readline/history.h> // add_history()
#endif

#include <stdlib.h> // free()
#include "readline.h"

/* A static variable for holding the line. */
static char *line_read = (char *)NULL;

/* Read a string, and return a pointer to it.  Returns NULL on EOF. */
char * rl_gets (char * prompt)
{
  /* If the buffer has already been allocated, return the memory
     to the free pool. */
  if (line_read)
    {
      free (line_read);
      line_read = (char *)NULL;
    }
  
  /* Get a line from the user. */
  line_read = readline (prompt);
  
#ifndef NO_READLINE_HISTORY
  /* If the line has any text in it, save it on the history. */
  if (line_read && *line_read)
    add_history (line_read);
#endif
  
  return (line_read);
}

