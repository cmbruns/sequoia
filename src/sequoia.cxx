#include "parse.hxx"
#include "oshook.hxx"
#include "timelice.hxx"
#include "timestmp.h"
extern "C" {
#include "readline.h"
}

// $Id: sequoia.cxx,v 1.5 2003/01/04 17:24:05 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/sequoia.cxx,v 1.5 2003/01/04 17:24:05 bruns Exp $
// $Log: sequoia.cxx,v $
// Revision 1.5  2003/01/04 17:24:05  bruns
// Changed timestmp.h to meet 8.3 filename for MS-DOS
// Changed copyright to 2003
//
// Revision 1.4  2002/06/21 00:45:42  bruns
// Added printing of compilation time/host stamp to program output header
//
// Revision 1.3  2001/12/04 23:42:59  bruns
// Updated copyright years and documentation URL
//
// Revision 1.2  2001/11/29 00:04:43  bruns
// Removed ^M characters
// Added cvs tags
//

cli cmd_line("SEQUOIA> ");

int main()
{

  status_os << "\nSEQUOIA multiple sequence alignment tool\n";
  status_os << "version " << VERSION << "\n";
  status_os << "compiled on " << timestamp << "\n";
  status_os << "copyright (c) 1995-2003\nby Chris Bruns, Ph.D.\n";
  status_os << "All rights reserved\n";
  status_os << "\n*** Please see http://bruns.homeip.net/~bruns/scripps/sequoia.html for ***\n";
  status_os << "*** proper citation information. ***\n\n";
  myoshook(status_os);

  int license_status;
  if ((license_status = timed_license()))
    {
      if (license_status == 1)
        {
          cerr << "Temporary license expired.  Please get latest version." << endl;
          exit(1);
        }
      else if (license_status == -1)
        {
          cerr << "Temporary license not yet active." << endl;
          exit(1);
        }
      else
        {
          cerr << "License verification problem." << endl;
          exit(1);
        }
    }

  make_command_hierarchy();

  // char buffer[2002];
  char * buffer;
  while(1)
    {
      status_os << endl;
      myoshook(status_os);
      buffer = rl_gets("SEQUOIA> ");
      // if stream input gives EOF, exit gracefully
      if (buffer == NULL) buffer = "QUIT";
      // cin.getline(buffer, 2000);
      cmd_line.parse_and_run2(buffer, status_os);
      myoshook(status_os);
    }
}
