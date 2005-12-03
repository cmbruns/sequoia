#include <iostream.h>
#include "parse.hxx"
#include "variable.hxx"

// $Id: vartest.cxx,v 1.2 2001/11/29 00:04:43 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/vartest.cxx,v 1.2 2001/11/29 00:04:43 bruns Exp $
// $Log: vartest.cxx,v $
// Revision 1.2  2001/11/29 00:04:43  bruns
// Removed ^M characters
// Added cvs tags
//

main()
{
  cout << "Hello world\n";
  assign_variable("phlerb", 3);
  sequoia_variable * vp = sq_var_ptr("phlerb");
  if (vp == NULL) cout << "NULL\n";
  else cout << vp->value() << "\n";
}
