#include "variable.hxx"

// $Id: variable.cxx,v 1.4 2002/06/28 17:17:46 bruns Exp $
// $Log: variable.cxx,v $
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
// $Header: /usr/data/cvs/sequoia1/variable.cxx,v 1.4 2002/06/28 17:17:46 bruns Exp $

#include "array.cxx"
template class Array<sequoia_variable>;
static Array<sequoia_variable> var; // array of all variables

int variable_matches(const sequoia_variable & v, const Mystring s);
// int variable_matches(const sequoia_variable & v, const char * s);
int add_sq_variable(const char * s, double v);

int list_variables(ostream & os)
{
  unsigned int i;
  for (i=0; i < var.dim(); ++i)
    os << var[i].name() << " = " << var[i].value() << "\n";
  return 0;
}


double assign_variable(const char * s, double v)
{
  unsigned int i;
  for (i=0; i < var.dim(); ++i) // see if variable exists
    {
      if (variable_matches(var[i], s))
	{
	  var[i].value() = v;
	  return v;
	}
    }
  add_sq_variable(s, v);
  return v;
}

int variable_matches(const sequoia_variable & v, const Mystring s)
{
  if (v.name().upcase() == s.upcase()) return TRUE;
  else return FALSE;
}

bool variable_exists(const char * s)
{
  unsigned int i;
  for (i=0; i < var.dim(); ++i) // see if variable exists
    {
      if (variable_matches(var[i], s))
	{
	  return true;
	}
    }
  return false;
}

sequoia_variable * sq_var_ptr(const char * s)
{
  unsigned int i;
  for (i=0; i < var.dim(); ++i)
    {
      if (variable_matches(var[i], s))
	{
	  return & var[i];
	}
    }
  return NULL;
}

int add_sq_variable(const char * s, double v)
{
  sequoia_variable sv(s, v);
  var += sv;
  var[var.dim()-1].name() = s;
  var[var.dim()-1].value() = v;
  return 0;
}


// member functions

// copy constructor
sequoia_variable::sequoia_variable(const sequoia_variable & sv1)
{
  sequoia_variable & sv = *this;
  sv.internal_value = sv1.internal_value;
  sv.internal_name = sv1.internal_name;
  sv.default_value = sv1.default_value;
  // if value pointer is internal, keep it internal
  if (sv1.value_ptr == &(sv1.internal_value))
    sv.value_ptr = &(sv.internal_value);
  else sv.value_ptr = sv1.value_ptr;
}


// assignment operator
sequoia_variable & sequoia_variable::operator=(const sequoia_variable & sv1)
{
  if (this == &sv1) return *this;
  sequoia_variable & sv = *this;
  sv.internal_value = sv1.internal_value;
  sv.internal_name = sv1.internal_name;
  sv.default_value = sv1.default_value;
  // if value pointer is internal, keep it internal
  if (sv1.value_ptr == &(sv1.internal_value))
    sv.value_ptr = &(sv.internal_value);
  else sv.value_ptr = sv1.value_ptr;
  return sv;
}

