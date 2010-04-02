#ifndef PARSE_HXX
#define PARSE_HXX

// $Id: parse.hxx,v 1.4 2002/06/28 17:21:54 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/parse.hxx,v 1.4 2002/06/28 17:21:54 bruns Exp $
// $Log: parse.hxx,v $
// Revision 1.4  2002/06/28 17:21:54  bruns
// Minor changes to avoid compiler warnings with gcc 3.1
//   changed <iostream.h> style includes to modern <iostream> style
//   removed redundant default function parameter values when already specified in
//    the header
//
// Revision 1.3  2001/12/12 18:26:27  bruns
// Many changes to get compilation with gcc 3.0.2
// Removed all templated friend functions (no friends needed)
// Added .c_str() method to Mystring, and removed auto conversion to const char *
//
// Revision 1.2  2001/11/28 23:40:22  bruns
// Added cvs header tags
// Removed ^M characters
//

#ifdef IRIX
#include "stdlib.h"
#else
#include <cstdlib>
#endif

#include <fstream>
#include <cstdio>
#include <sstream>
#include "mystring.hxx"
#include "matrix.hxx"
#include "sequence.hxx"
#include "comp_mat.hxx"
#include "path_mat.hxx"
#include "brookhav.hxx"
#include "overlay.hxx"
#include "consens.hxx"

class command_token;

// Command Line Interface
class cli
{
friend class command_token;
friend ostream & operator<<(ostream & os, cli & c);
friend istream & operator>>(istream & is, cli & c);
private:
  Mystring prompt_string;
  Comparison_matrix cm;
public:
  Score_matrix sm;
  Alignment align[3];
  pdboverlay strct[3];
  cli(const Mystring & s);
  const Mystring & prompt(void) const;
  int parse_and_run(const char * buffline, ostream & os);
  int parse_and_run2(const char * buffline, ostream & os);
};


class command_token
{
protected:
  Mystring name_string;
  int (command_token::*myjob)
    (ostream & os); // actual function of the token
  command_token * parent;
public:
  command_token(); // default constructor
  int run(const char * line, ostream & os); // how to start it
  const char * rest_of_line; // rest of line after this token
  command_token & 
    add_subtoken(const char *,
		 int (command_token::*fn)(ostream & os));
  int add_synonym(const char *);
  const char * name() const {return name_string.c_str();}
  Mystring & name() {return name_string;}
  Array<Mystring> synonym;
  Mystring help_string;
  Array<command_token> subtoken; // tokens that can come after this
  bool takes_an_input_filename_argument;
  bool takes_an_output_filename_argument;
  bool prepend_token; // can it be directly prepended to next token?
  istream * infile;
  ostream * outfile;

  int evaluate_expression(ostream & os); // SET something
  
  // Prototypes for token action functions
  // 
  // Linux does not like setting 'myjob' to NULL, so this function
  // is there to replace NULL
  int has_no_function(ostream & os); // no assigned function
  // 
  int align_sequences(ostream & os); // ALIGN
  int do_nothing(ostream & os); // COMMENT
  int make_consensus(ostream & os); // CONSENSUS
  int capitalize_equivalences(ostream & os); // EQUIV
  int capitalize_equivalences2(ostream & os); // EQUIV2
  int start_help(ostream & os); // HELP
  int optimize(ostream & os); // OPTIMIZE
  int overlay_structures(ostream & os); // OVERLAY
  int print_sequence_file(ostream & os); // PRINT <register>
  int print_id(ostream & os); // PRINT ID
  int print_matrix(ostream & os); // PRINT MATRIX
  int quit_program(ostream & os); // QUIT
  int random(ostream & os); // RANDOM
  int read_sequence_file(ostream & os); // READ
  int run_command(ostream & os); // RUN
  int salign_sequences(ostream & os); // SALIGN
  int set_variable(ostream & os); // SET
  int split_sequences(ostream & os); // SPLIT
  int stabulate(ostream & os); // STABULATE
  int superpose(ostream & os); // SUPERPOSE
  int system_command(ostream & os); // SYSTEM
  int tabulate(ostream & os); // TABULATE
  int weight(ostream & os); // WEIGHT
  int write_sequence_file(ostream & os); // WRITE

};

command_token & make_command_hierarchy(void);
const char * advance_to_next_token(const char * s);

#endif
