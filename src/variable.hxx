#ifndef _VARIABLE_HXX_
#define _VARIABLE_HXX_

#include <iostream>
#include "mystring.hxx"

// $Id: variable.hxx,v 1.4 2002/06/28 17:17:46 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/variable.hxx,v 1.4 2002/06/28 17:17:46 bruns Exp $
// $Log: variable.hxx,v $
// Revision 1.4  2002/06/28 17:17:46  bruns
// Changed to avoid warnings with gcc 3.1
//   changed <iostream.h> type includes to <iostream>
//   removed redundant function parameter default values from .cxx files
//     (when they are already defined in the header)
//
// Revision 1.3  2001/12/04 23:43:35  bruns
// Made first argument to assign_variable "const"
//
// Revision 1.2  2001/11/29 00:04:43  bruns
// Removed ^M characters
// Added cvs tags
//

class sequoia_variable
{
friend int add_sq_variable(const char * s, double v=0.0);
protected:
  double internal_value;
  Mystring internal_name;
  double default_value;
  double * value_ptr; // use this, so we can also point to other stuff
  sequoia_variable(Mystring n, double dv=0.0)
    {
      value_ptr = &internal_value; 
      name() = n; 
      value() = dv;
    }
public:
  sequoia_variable(double dv=0.0)
    {
      value_ptr = &internal_value; 
      value() = dv;
    }
  // need a copy constructor to handle value_ptr properly
  sequoia_variable(const sequoia_variable & sv1);
  sequoia_variable & operator=(const sequoia_variable & sv1);

  Mystring & name() {return internal_name;}
  const Mystring & name() const {return internal_name;}
  double & value() {return *value_ptr;};
  const double & value() const {return *value_ptr;};
  void set_value_pointer(double * dp) {value_ptr = dp;}
};

int list_variables(ostream & os);
sequoia_variable * sq_var_ptr(const char * s);
double assign_variable(const char * s, double v);
bool variable_exists(const char * s);

#endif
