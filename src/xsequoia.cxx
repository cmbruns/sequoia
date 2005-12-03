#include "strstream.h"
#include "parse.hxx"
#include "xsequoia.hxx"

// $Id: xsequoia.cxx,v 1.2 2001/11/29 00:04:43 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/xsequoia.cxx,v 1.2 2001/11/29 00:04:43 bruns Exp $
// $Log: xsequoia.cxx,v $
// Revision 1.2  2001/11/29 00:04:43  bruns
// Removed ^M characters
// Added cvs tags
//

extern "C" {
#include "mtwin.h"
extern void init_windows(int argc, char *argv[]);
extern XtAppContext app;
}

ostrstream & os = status_os;

cli cmd_line("XSEQUOIA> ");

main(int argc, char *argv[])
{

  os << "XSEQUOIA multiple sequence alignment tool\n";
  os << "version " << VERSION << "\n";
  os << "copyright (c) 1995-1999\nby Chris Bruns, Ph.D.\n";
  os << "The Scripps Research Institute, La Jolla, CA\n";
  os << "All rights reserved\n";

  init_windows(argc, argv);
  
  myoshook(os);

  XtAppMainLoop(app);

  exit(0);
}

