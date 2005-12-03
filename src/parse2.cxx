// parse2.cxx
// object oriented command parsing

// $Id: parse2.cxx,v 1.7 2002/06/28 17:21:54 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/parse2.cxx,v 1.7 2002/06/28 17:21:54 bruns Exp $
// $Log: parse2.cxx,v $
// Revision 1.7  2002/06/28 17:21:54  bruns
// Minor changes to avoid compiler warnings with gcc 3.1
//   changed <iostream.h> style includes to modern <iostream> style
//   removed redundant default function parameter values when already specified in
//    the header
//
// Revision 1.6  2001/12/12 18:26:27  bruns
// Many changes to get compilation with gcc 3.0.2
// Removed all templated friend functions (no friends needed)
// Added .c_str() method to Mystring, and removed auto conversion to const char *
//
// Revision 1.5  2001/12/05 00:46:48  bruns
// Added average percent ID column to print ID output
//
// Revision 1.4  2001/12/05 00:26:21  bruns
// Added additional syntax for weight command
// WEIGHT REGISTER INDEX VALUE
// weight seq1 2 3.0
// sets the weight of a particular sequence manually
//
// Revision 1.3  2001/12/05 00:01:49  bruns
// Added line of output to beginning of PRINT output, to say the number
// of sequences present.
//
// Revision 1.2  2001/11/28 23:40:22  bruns
// Added cvs header tags
// Removed ^M characters
//

#include <iomanip>
#include "array.hxx"
#include "mystring.hxx"
#include "oshook.hxx"
#include "parse.hxx"
#include "variable.hxx"

extern cli cmd_line;
// cli cmd_line("vartest> ");

static double pretty_length = 50;
static command_token toplevel; // the ur-command_token
static bool is_help = FALSE; // TRUE when we want help rather than execution

// expressions for variable parsing
// primary expressions, can be constants, variables, or assignments 
class sequoia_expression
{
protected:
public:
  Mystring expression_string;
  bool is_constant(); // e.g. a number
  bool is_variable();
  bool is_lvalue(); // can its value be assigned?  (e.g. a variable)
  bool is_assignment(); // e.g. a = 3
  double evaluate();
  int strip_parentheses();  
};


// Other function prototypes
// PARSING
const char * discard_leading_whitespace(const char * s);
int command_matches(command_token & c, const char * s);
Mystring get_one_token(const char * s);
int print_help(const command_token & c, ostream & os);

// This function defines the parsing behavior
command_token & make_command_hierarchy(void)
{

  command_token *new_token, *new_token2;

  // The "toplevel" token precedes the first actual command.
  // In other words, it is the root of the command tree.
  toplevel.help_string = "SEQUOIA is a tool for the alignment of "
    "biological sequences and atomic coordinates.";
  toplevel.name() = "toplevel";

  // ALIGN
  new_token = & toplevel.add_subtoken("ALIGN",
				      &command_token::align_sequences);
  new_token->help_string = "ALIGN aligns the sequences in register "
    "SEQ1 with those in SEQ2, and places the result in the ALIGN "
      "register\n  Usage: ALIGN";

  // COMMENT
  new_token = & toplevel.add_subtoken("COMMENT",
				      &command_token::do_nothing);
  new_token->add_synonym("#");
  new_token->add_synonym("[");
  new_token->add_synonym("{");
  new_token->add_synonym("/");
  new_token->add_synonym("REM");
  new_token->prepend_token = TRUE;
  new_token->help_string = "A line starting with COMMENT is ignored.";

  // CONSENSUS
  new_token = & toplevel.add_subtoken("CONSENSUS",
				      &command_token::make_consensus);
  new_token->help_string = "CONSENSUS <sigma> generates a consensus "
    "sequence based on the contents of ALIGN, and places the result "
      "in SEQ2.  The <sigma> determines how many poorly conserved "
	"residues will be represented by an `X` character.";

  // EQUIVALENCE
  new_token = & toplevel.add_subtoken("EQUIVALENCE",
				      &command_token::capitalize_equivalences2);
  new_token->help_string = "EQUIVALENCE capitalizes those residue codes in "
    "the ALIGN register whose atomic coordinates are structurally "
    "equivalent.  Others are converted to lower case.  EQUIVALENCE is one "
    "step in the OVERLAY process.\n"
    "Usage: EQUIVALENCE";

  // EQUIVALENCE2
  new_token = & toplevel.add_subtoken("EQUIVALENCE2",
				      &command_token::capitalize_equivalences2);
  new_token->help_string = "EQUIVALENCE capitalizes those residue codes in "
    "the ALIGN register whose atomic coordinates are structurally "
    "equivalent.  Others are converted to lower case.  EQUIVALENCE is one "
    "step in the OVERLAY process.\n"
    "Usage: EQUIVALENCE2";

  // HELP
  new_token = & toplevel.add_subtoken("HELP", &command_token::start_help);
  new_token->add_synonym("?");
  new_token->add_synonym("MAN");
  new_token->help_string = "HELP prints documentation about SEQUOIA " 
    "commands.\n  Usage: HELP <command>\n"
      "  EXAMPLE: HELP ALIGN\n";

  // OPTIMIZE
  new_token = & toplevel.add_subtoken("OPTIMIZE",
				      &command_token::optimize);
  new_token->help_string = "OPTIMIZE removes some progressive "
    "alignment artifacts from the ALIGN register by reALIGNing "
      "each individual sequence in turn\n";

  // OVERLAY
  new_token = & toplevel.add_subtoken("OVERLAY",
				      &command_token::overlay_structures);
  new_token->help_string = "OVERLAY iteratively superposes and "
    "realigns two coordinate files to convergence.  OVERLAY calls the "
      "EQUIVALENCE, SUPERPOSE, and SALIGN commands.\n";


  // PRINT
  new_token = & toplevel.add_subtoken
    ("PRINT", &command_token::print_sequence_file); 
  new_token->help_string = "PRINT is used to display data from "
    "SEQUOIA.\nPRINT <register> displays the sequence or structure in "
      "the register.\nPRINT MATRIX displays the current comparison matrix";
  new_token2 = & new_token->add_subtoken
    ("ID", &command_token::print_id); 
  new_token2->help_string = "PRINT ID <register> displays a "
    "table of sequence identity.";
  new_token2 = & new_token->add_subtoken
    ("MATRIX", &command_token::print_matrix); 
  new_token2->help_string = "PRINT MATRIX displays the current "
    "comparison matrix"; 

  // QUIT
  new_token = & toplevel.add_subtoken("QUIT", &command_token::quit_program);
  new_token->add_synonym("BYE");
  new_token->add_synonym("END");
  new_token->add_synonym("EXIT");
  new_token->add_synonym("HALT");
  new_token->add_synonym("KILL");
  new_token->add_synonym("Q");
  new_token->add_synonym("STOP");
  new_token->help_string = "QUIT causes SEQUOIA to terminate.";

  // READ
  new_token = & toplevel.add_subtoken
    ("READ", &command_token::read_sequence_file); 
  new_token->help_string = "READ is used to load data into "
    "SEQUOIA from files.\n"
    " READ <register> <filename> loads structure data into a register.\n"
    " READ MATRIX <filename> loads a comparison matrix.\n";

  // RANDOM
  new_token = & toplevel.add_subtoken
    ("RANDOM", &command_token::random); 
  new_token->help_string = "RANDOM <iterations> calculates alignment "
    "scores with one register randomly shuffled.  This is used to "
      "provide one estimate of alignment significance\n";
  
  // RUN
  new_token = & toplevel.add_subtoken("RUN", 
                 &command_token::run_command);
  new_token->add_synonym("@");
  new_token->prepend_token = TRUE;
  // For recursion reasons, we must open the file in the action
  // function
  // Therefore, "run" should not try to get the filename
  // new_token2->takes_an_input_filename_argument = TRUE;
  new_token->help_string = "RUN <filename> executes SEQUIOA commands "
    "contained in a file.";

  // SALIGN
  new_token = & toplevel.add_subtoken("SALIGN",
				      &command_token::salign_sequences);
  new_token->help_string = "SALIGN aligns the sequences in registers "
    "SEQ1 with those in SEQ2 based upon the superposed atomic "
      "coordinates , and places the result in the ALIGN "
      "register\n  Usage: SALIGN";

  // SET
  new_token = & toplevel.add_subtoken("SET", 
                 &command_token::set_variable);
  new_token->add_synonym("LET");
  new_token->help_string = "SET <variable> <value> assigns a value to "
    "a SEQUOIA variable.\n  EXAMPLE: SET GPEN 5.0\n"
    "SET <register1> <register2> copies the contents of <register1> "
    "into <register2>\n  EXAMPLE: SET SEQ1 ALIGN\n"
    "SET with no arguments lists the values of all variables.\n";

  // SPLIT
  new_token = & toplevel.add_subtoken("SPLIT",
				      &command_token::split_sequences);
  new_token->help_string = "SPLIT divides the sequences in the ALIGN "
    "register between SEQ1 and SEQ2\n"
      "register\n  Usage: SPLIT <sequence number>";

  // STABULATE
  new_token = & toplevel.add_subtoken("STABULATE", 
                 &command_token::stabulate);
  new_token->help_string = "STABULATE computes the score matrix used by "
    "the SALIGN command.  (This is ordinarily done automatically.)";

  // SUPERPOSE
  new_token = & toplevel.add_subtoken("SUPERPOSE", 
                 &command_token::superpose);
  new_token->help_string = "SUPERPOSE transforms the coordinates in "
    "STRUCT2 onto STRUCT1 using the capitalized residues in ALIGN "
      "and places the result into OVERLAY\n";

  // SYSTEM
  new_token = & toplevel.add_subtoken("SYSTEM", 
                 &command_token::system_command);
  new_token->add_synonym("!");
  new_token->prepend_token = TRUE;
  new_token->help_string = "SYSTEM <command> executes a command in the "
    "local operating system.\n  EXAMPLE: SYSTEM ls";

  // TABULATE
  new_token = & toplevel.add_subtoken("TABULATE", 
                 &command_token::tabulate);
  new_token->help_string = "TABULATE computes the score matrix used by "
    "the ALIGN command.  (This is ordinarily done automatically.)";

  // WEIGHT
  new_token = & toplevel.add_subtoken("WEIGHT", 
                 &command_token::weight);
  new_token->help_string = "WEIGHT applies a weight to each sequence "
    "to compensate for unequal patterns of sequence similarity\n";

  // WRITE
  new_token = & toplevel.add_subtoken("WRITE", 
                 &command_token::write_sequence_file);
  new_token->help_string = "WRITE is used to create files from "
    "SEQUOIA.\nWRITE <register> <filename> creates a sequence or "
      "coordinate file.\n";

  // ********************************* //

  // INITIALIZE VARIABLES
  // these initial values that get re-pointered won't be used
  sequoia_variable * sv;
  // ACUTOFF
  assign_variable("acutoff", 15.0); // create variable
  sv = sq_var_ptr("acutoff"); // get pointer
  sv->set_value_pointer(&cmd_line.sm.acutoff()); // set value to
  // DCUTOFF
  assign_variable("dcutoff", 4.5); // create variable
  sv = sq_var_ptr("dcutoff"); // get pointer
  sv->set_value_pointer(&cmd_line.sm.dcutoff()); // set value to
  // ECHO
  assign_variable("echo", 0.0); // create variable
  // EPEN
  assign_variable("epen", 1.0); // create variable
  sv = sq_var_ptr("epen"); // get pointer
  sv->set_value_pointer(&cmd_line.sm.epen()); // set value to external
  // GPEN
  assign_variable("gpen", 10.0); // create variable
  sv = sq_var_ptr("gpen"); // get pointer
  sv->set_value_pointer(&cmd_line.sm.gpen()); // set value to external
  // PRETTY_LENGTH
  assign_variable("pretty_length", 0); // create variable
  sv = sq_var_ptr("pretty_length"); // get pointer
  sv->set_value_pointer(&pretty_length); // set value to external
  // RANDOM_SEED
  assign_variable("random_seed", 0); // create variable
  sv = sq_var_ptr("random_seed"); // get pointer
  sv->set_value_pointer(&cmd_line.sm.random_seed()); // set value to external
  // RUNLENGTH
  assign_variable("runlength", 0); // create variable
  sv = sq_var_ptr("runlength"); // get pointer
  sv->set_value_pointer(&cmd_line.sm.runlength()); // set value to external
  // SUBOPTIMAL
  assign_variable("suboptimal", 0); // create variable
  sv = sq_var_ptr("suboptimal"); // get pointer
  sv->set_value_pointer(&cmd_line.sm.suboptimal()); // set value to external
  // USEANGLE
  assign_variable("useangle", 0); // create variable
  sv = sq_var_ptr("useangle"); // get pointer
  sv->set_value_pointer(&cmd_line.sm.useangle()); // set value to external

  return toplevel;
}


int cli::parse_and_run2(const char * buffline, ostream & os)
{
  is_help = FALSE; // TRUE when we want help rather than execution
  const char * newpos = discard_leading_whitespace(buffline);
  // blank lines should be ignored
  if (strlen(newpos) < 1) {return 0;}
  toplevel.rest_of_line = newpos;
  sequoia_variable * svp = NULL;
  if ((svp = sq_var_ptr("ECHO")) != NULL)
    if (svp->value()) os << newpos << endl;
  return toplevel.run(newpos, os);
}

// **** Operations on command_tokens follow ****

// "run" recursively plows through command tokens
// "line" has the current token removed
int command_token::run(const char * line, 
		       ostream & os) // run the command 
{
  const char * line2 = line;
  // first try to find a filename, if we want one
  if (takes_an_input_filename_argument)
    {
      Mystring fname = get_one_token(line2);
      line2 = advance_to_next_token(line2);
      if (!is_help)
	{
	  if (strlen(fname.c_str()) < 1)
	    {
	      os << "!!! SYNTAX ERROR !!! (not enough arguments)\n"; 
	      return 1;
	    }
	  else
	    {
	      // this may be a bit funny, keeping a pointer to this
	      // static variable
	      static ifstream inf;
	      inf.close();
	      inf.open(fname.c_str(), ios::in);
	      infile = & inf;
	    }
	}
    }
  if (takes_an_output_filename_argument)
    {
      Mystring fname = get_one_token(line2);
      line2 = advance_to_next_token(line2);
      if (!is_help)
	{
	  if (strlen(fname.c_str()) < 1)
	    {
	      outfile = & os;
	    }
	  else
	    {
	      // this may be a bit funny, keeping a pointer to this
	      // static variable
	      static ofstream of;
	      of.close();
	      of.open(fname.c_str(), ios::out);
	      outfile = & of;
	    }
	}
    }
  if (strlen(line2) < 1) // we are at end of command line
    {
      if (is_help) 
	{
	  return print_help(*this, os);
	  return 0;
	}
      // else if (myjob != NULL) // no subtokens are taken, does action
      else if (myjob != &command_token::has_no_function) // no subtokens are taken, does action
	return ((this->*myjob)(os)); // run action function
      else 
	{
	  os << "!!! SYNTAX ERROR !!! (not enough arguments)\n"; 
	  return 1;
	}
    }
  else if (subtoken.dim() > 0) // check for subtokens on rest of line
    {
      unsigned int i;
      for (i=0; i < subtoken.dim(); ++i)
	{
	  if (command_matches(subtoken[i], line2))
	    {
	      return subtoken[i].run(subtoken[i].rest_of_line, os);
	    }
	}
    }
  // if we reach this point, there are no matching subtokens
  if (is_help) 
    {
      return print_help(*this, os);
    }
  // else if (myjob != NULL) // is there an action function?
  else if (myjob != &command_token::has_no_function) // is there an action function?
    return (this->*myjob)(os); // run action function
  else os << "!!! SYNTAX ERROR !!! (too many (or wrong) arguments)\n";
  return 1; 
}


int command_matches(command_token & c, const char * s)
{
  // for now, just go for exact matches
  // (later, permit some tokens to have partial matches)
  Mystring s2 = get_one_token(s);

  // first check "official" command name
  if (c.prepend_token) 
    {
      if (!strncmp(c.name().c_str(),s2.upcase().c_str(),strlen(c.name().c_str()))) 
	{
	  c.rest_of_line = s + strlen(c.name().c_str());
	  return TRUE;
	}
    }
  else if (!strcmp(c.name().c_str(),s2.upcase().c_str())) 
    {
      c.rest_of_line = advance_to_next_token(s);
      return TRUE;
    }
  unsigned int i;
  // now check synonyms
  for (i=0; i < c.synonym.dim(); ++i)
    {
      if (c.prepend_token) 
	{
	  if (!strncmp(c.synonym[i].c_str(),s2.upcase().c_str(),strlen(c.synonym[i].c_str()))) 
	    {
	      c.rest_of_line = s + strlen(c.synonym[i].c_str());
	      return TRUE;
	    }
	}
      else if (!strcmp(c.synonym[i].c_str(),s2.upcase().c_str())) 
	{
	  c.rest_of_line = advance_to_next_token(s);
	  return TRUE;
	}
    }
  return FALSE;
}


command_token::command_token() // default constructor
{
  // myjob = NULL;
  myjob = &command_token::has_no_function;
  takes_an_input_filename_argument = FALSE;
  takes_an_output_filename_argument = FALSE;
  prepend_token = FALSE;
  infile = NULL;
  outfile = NULL;
  parent = NULL;
}


int command_token::add_synonym(const char * s)
{
  synonym += s;
  return 0;
}


command_token & command_token::add_subtoken( const char * s,
			  int (command_token::*fn)(ostream &))
{
  command_token newtoken;
  newtoken.name_string = s;
  newtoken.myjob = fn;
  newtoken.parent = this;
  subtoken += newtoken;
  return subtoken[subtoken.dim() - 1];
}


int print_help(const command_token & c, ostream & os)
{
  int result = 0;

  os << "\n";
  if (strlen(c.help_string.c_str()) < 1) 
    {
      os << "No help available\n";
      result = 1;
    }
  else
    {
      os << c.help_string << "\n";
    }
  if (c.synonym.dim() > 0)
    {
      os << "SYNONYMS:  ";
      unsigned int i;
      os << c.name();
      for (i=0; i < c.synonym.dim(); ++i)
	os << ",  " << c.synonym[i];
      os << "\n";
    }
  if (c.subtoken.dim() > 0)
    {
      os << "SUBCOMMANDS: \n";
      unsigned int i;
      for (i=0; i < c.subtoken.dim(); ++i)
	os << "  " << c.subtoken[i].name() << "\n";
    }
  return result;
}

int command_token::evaluate_expression(ostream & os) // SET something
{
  // at first, just accept "variable value" expressions
  Mystring varname = get_one_token(rest_of_line);
  if (varname.length() < 1) list_variables(os); // no arguments, list
						// everything  
  else
    {
      rest_of_line = advance_to_next_token(rest_of_line);
      double value = atof(get_one_token(rest_of_line).c_str());
      rest_of_line = advance_to_next_token(rest_of_line);
      assign_variable(varname.c_str(), value);
    }
  return 0;
}

// **** Parsing operations on strings follow: ****

// return a pointer to the first non-whitespace character of string
const char * discard_leading_whitespace(const char * s)
{
  const char * s2 = s;
  while (isspace(*s2)) ++s2;
  return s2;
}

// return a pointer to the next token of string
const char * advance_to_next_token(const char * s)
{
  const char * s2 = discard_leading_whitespace(s);
  while (isgraph(*s2)) ++s2;
  return discard_leading_whitespace(s2);
}

Mystring get_one_token(const char * s)
{
  Mystring answer = "";
  const char * s2 = discard_leading_whitespace(s);
  while (isgraph(*s2)) 
    {
      answer += (char) *s2;
      ++s2;
    }
  return answer;
}

bool is_a_register_name(const char * s)
{
  Mystring s2 = s;
  if (!strcmp(s2.upcase().c_str(), "SEQ1")) return TRUE;
  else if (!strcmp(s2.upcase().c_str(), "SEQ2")) return TRUE;
  else if (!strcmp(s2.upcase().c_str(), "ALIGN")) return TRUE;
  else if (!strcmp(s2.upcase().c_str(), "STRUCT1")) return TRUE;
  else if (!strcmp(s2.upcase().c_str(), "STRUCT2")) return TRUE;
  else if (!strcmp(s2.upcase().c_str(), "OVERLAY")) return TRUE;
  else return false;
}

bool is_a_number(const char * s)
{
  if (atof(s)) return true;
  // must be zero or not a number
  else if (s[0] == '0') return true;
  else if (s[0] == '.') return true;
  else if (s[0] == '-') return true;
  else return false;
}

// 
// Put action functions for commands here, at end of this file
// 

// Linux does not like setting 'myjob' to NULL, so this function
// is there to replace NULL
int command_token::has_no_function(ostream & os) // no assigned function
{
  os << "This function should not be called!!!" << endl;
  os << "(has_no_function)" << endl;
  return 1;
}
// 

int command_token::align_sequences(ostream & os) // ALIGN
{
      if (!(cmd_line.align[1].n_seq() && cmd_line.align[2].n_seq()))
        {
          os << "Need two sequences/alignments to align!!!\n";
          os << "(Use READ command)\n";
          return 1;
        }
      if (!cmd_line.sm.have_table) 
	{
	  os << "Running TABULATE command...\n";
          myoshook(os);
	  char buffer[20];
	  strcpy(buffer, "TABULATE");
	  if (cmd_line.parse_and_run2(buffer, os)) return 1;
	}
      os << "Aligning sequences...\n";
      myoshook(os);
      cmd_line.align[0] = cmd_line.sm.dyn_prog();
      os << cmd_line.align[0].newgaps() << " gaps added to make alignment\n";
      os << "Alignment score is " << cmd_line.align[0].score_val << "\n";
      cmd_line.sm.have_align = true;
      return 0;
}

int command_token::do_nothing(ostream & os) // COMMENT
{
  return 0;
}

int command_token::make_consensus(ostream & os) // CONSENSUS
{
  double cutoff = 1.2; // default sigma cutoff
  if (strlen(get_one_token(rest_of_line).c_str()))
    {
      cutoff = atof(get_one_token(rest_of_line).c_str());
      os << "Cutoff is " << cutoff << " sigma\n";
    }
  else os << "Using default cutoff of " << cutoff << " sigma\n";
  Sequence cons = consensus(cmd_line.align[0], cmd_line.cm, cutoff);
  cons.title() = "Consensus";
  cmd_line.align[2] = cons.remove_gaps();
  os << "Consensus of ALIGN has been placed in SEQ2\n";
  cmd_line.sm.have_seq2 = true;
  return 0;
}

// This is the new version that is not done yet (apparently)
int command_token::capitalize_equivalences(ostream & os) // EQUIVALENCE
{
  os << "Reassigning equivalences...\n";
  myoshook(os);
  uint i;
  i = cmd_line.align[0][0].n_cap_aligned(cmd_line.align[0][1]);
  os << "Total starting equivalences = " << i << "\n";
  cmd_line.align[0] = struct_equiv(cmd_line.align[0],
				   cmd_line.sm.dcutoff(),
				   cmd_line.sm.acutoff(),
				   (unsigned int) cmd_line.sm.runlength());
  i = cmd_line.align[0][0].n_cap_aligned(cmd_line.align[0][1]);
  os << "Total new equivalences = " << i << "\n";
  return 0;
}

int command_token::capitalize_equivalences2(ostream & os) // EQUIVALENCE2
{
  os << "Reassigning equivalences...\n";
  myoshook(os);
  uint i;
  i = cmd_line.align[0][0].n_cap_aligned(cmd_line.align[0][1]);
  os << "Total starting equivalences = " << i << "\n";
  // if (!cmd_line.sm.have_stable)
  //   if (cmd_line.parse_and_run2("STABULATE",os)) return 1;
  cmd_line.align[0] = struct_equiv(cmd_line.strct[1],cmd_line.strct[0],
			  cmd_line.align[0],cmd_line.sm.dcutoff(),
			  (unsigned int) cmd_line.sm.runlength());
  i = cmd_line.align[0][0].n_cap_aligned(cmd_line.align[0][1]);
  os << "Total new equivalences = " << i << "\n";
  return 0;
}

int command_token::start_help(ostream & os) // HELP
{
  is_help = TRUE; // TRUE when we want help rather than execution
  const char * newpos = discard_leading_whitespace(rest_of_line);
  toplevel.rest_of_line = newpos;
  return toplevel.run(newpos, os);
}

int command_token::optimize(ostream & os) // OPTIMIZE
{
      if (cmd_line.align[0].n_seq() < 2)
	{
	  os << "At least two sequences required to optimize!!\n";
	  return 1;
	}

      // randomize selection of residues
      int random_order[cmd_line.align[0].n_seq()];
      unsigned int i;
      // start with a definite order
      for (i = 0; i < cmd_line.align[0].n_seq(); ++i)
	random_order[i] = i;
      // swap each member of order list with another random member
      for (i = 0; i < cmd_line.align[0].n_seq(); ++i)
        {
	  int p, t;
	  // choose random partner
	  p = rand() % cmd_line.align[0].n_seq();
	  // swap
	  t = random_order[i];
	  random_order[i] = random_order[p];
	  random_order[p] = t;
	}

      // perform realignment
      Sequence tempseq;
      for (i = 0; i < cmd_line.align[0].n_seq(); ++i)
        {
	  char cmd[50];
	  int j;
	  int p = random_order[i];
	  sprintf(cmd, "SPLIT %d", p+1);
	  os << "SPLITing out sequence " << p+1 << ":" << endl;
          myoshook(os);
	  if (cmd_line.parse_and_run2(cmd, os)) return 1;
	  os << "reALIGNing:\n";
          myoshook(os);
	  if (cmd_line.parse_and_run2("ALIGN", os)) return 1;
	  
	  // replace original order
	  tempseq = cmd_line.align[0][0];
	  for (j=0; j < p; ++j)
	    {
	      cmd_line.align[0][j] = cmd_line.align[0][j+1];
	    }
	  cmd_line.align[0][p] = tempseq;
        }
      return 0;
}

int command_token::overlay_structures(ostream & os) // OVERLAY
{
  Alignment old_align;
  // overlay struct 2 onto struct 1
  if (cmd_line.parse_and_run2("SUPERPOSE",os)) return 1;
  myoshook(os);
  uint cycles = 0;
  while ((old_align != cmd_line.align[0]) && (cycles < 10))
    {
      old_align = cmd_line.align[0];
      if (cmd_line.parse_and_run2("EQUIVALENCE",os)) return 1;
      myoshook(os);
      if (cmd_line.parse_and_run2("SUPERPOSE",os)) return 1;
      myoshook(os);
      ++ cycles;
    }
  if (old_align != cmd_line.align[0]) 
    os << "Failed to converge in " << cycles << " cycles\n";
  else os << "Equivalences converged in " << cycles << " cycles\n";
  os << "Realigning based upon final overlay\n";
  if (cmd_line.parse_and_run2("SALIGN",os)) return 1;
  if (cmd_line.parse_and_run2("EQUIVALENCE",os)) return 1;
  return 0;
}

int command_token::print_sequence_file(ostream & os) // PRINT
{
  // Parse two arguments after READ command
  Mystring arg1 = get_one_token(rest_of_line);
  rest_of_line = advance_to_next_token(rest_of_line);
  Mystring arg2 = get_one_token(rest_of_line);
  rest_of_line = advance_to_next_token(rest_of_line);
  Mystring arg3 = get_one_token(rest_of_line);

  // Print contents of register
  if (is_a_register_name(arg1.c_str()))
    {
      Alignment * al_ptr = NULL;
      pdboverlay * str_ptr = NULL;
      if (arg1.upcase() == "SEQ1") al_ptr = & cmd_line.align[1];
      else if (arg1.upcase() == "SEQ2") al_ptr = & cmd_line.align[2];
      else if (arg1.upcase() == "ALIGN") al_ptr = & cmd_line.align[0];
      else if (arg1.upcase() == "STRUCT1") str_ptr = & cmd_line.strct[1];
      else if (arg1.upcase() == "STRUCT2") str_ptr = & cmd_line.strct[2];
      else if (arg1.upcase() == "OVERLAY") str_ptr = & cmd_line.strct[0];
      else 
	{
	  os << "Internal error!!!, unknown register!!!\n";
	  os << arg1 << "\n";
	  return 1;
	}	

      // Open file, if necessary
      ofstream of;
      ostream * ofp;
      if (strlen(arg2.c_str()) < 1)
	ofp = & os;
      else
	{
	  of.open(arg2.c_str(), ios::out);
	  ofp = 	& of;
	}
	
      // Is it a sequence file?
      if (al_ptr == NULL)
	  os << "Use WRITE for coordinates" << endl;
      // Is the file OK?
      else if (ofp->good()) {
        *ofp << "There are " << al_ptr->n_seq() << " sequences in the register." << endl;
        *ofp << endl;
	al_ptr->pretty(*ofp, (unsigned int) pretty_length);
      }
      else 
	{
	  os << "Could not open file\n";
	  of.close();
	  return 1;
	}
      of.close();
    }


  // Print superposition statistics
  else if (arg1.upcase() == "STRUCTSTATS")
    {
      if (!(cmd_line.align[0].dim() && cmd_line.strct[0].ca().dim()))
	{
	  os << "Need some residues to compare!\n";
	  return 1;
	}
      cmd_line.strct[0].update_refatoms();
      cmd_line.strct[0].update_orientations();
      cmd_line.strct[1].update_refatoms();
      cmd_line.strct[1].update_orientations();
      cmd_line.align[0].transfer_atoms(cmd_line.strct[1], cmd_line.strct[0]);
      uint i;
      // Just output first structure in each register
      uint off = cmd_line.strct[1].dim();
      // Make sure align has enough sequences
      if ((cmd_line.align[0].dim() < 1) || (cmd_line.align[0].dim() < off))
	{
	  os << "Not enough sequences in ALIGN\n";
	  return 1;
	}
      for (i=0; i < cmd_line.align[0][0].dim(); ++i)
	{
	  const SeqRes & r1 = cmd_line.align[0][0][i];
	  const SeqRes & r2 = cmd_line.align[0][off][i];
	  os << r1.letter() << " ";
	  os << r2.letter() << " ";
	  
	  // Print distance
	  os.width(8);
	  if ((r1.refatomp() == NULL) || (r2.refatomp() == NULL))
	    os << "    ----";
	  else os << r1.refatom().distance(r2.refatom());
	  os << "   ";
	  
	  // Print orientation change
	  os.width(8);
	  if ((r1.refatomp() == NULL) || (r2.refatomp() == NULL))
	    os << " No atom";
	  else if ((r1.orientation().dim() < 3) ||
		   (r2.orientation().dim() < 3))
	    {
	      os << " No matx " << r1.orientation().dim();
	      os << " " << r2.orientation().dim();
	    }
	  else 
	    {
	      float diffangle = r1.orientation().angle(r2.orientation());
	      os << setw(8);
	      os << setiosflags(ios::fixed);
	      os << diffangle * 180 / PI;
	    }
	  os << " \n";
	}
    }
  
  // Don't know what to print
  else
    {
      os << "SYNTAX ERROR: don't know how to print ``" << arg1 << "''"
	<< endl;
      return 1;
    }

  return 0;
}


  // Print table of identity
int command_token::print_id(ostream & os) // PRINT ID
{
  Mystring arg2 = get_one_token(rest_of_line);
  rest_of_line = advance_to_next_token(rest_of_line);
  Mystring arg3 = get_one_token(rest_of_line);
  rest_of_line = advance_to_next_token(rest_of_line);
  if (!is_a_register_name(arg2.c_str()))
    {
      os << "SYNTAX ERROR: Usage PRINT ID <register>\n";
      return 1;
    }
  Alignment * al_ptr = NULL;
  if (arg2.upcase() == "SEQ1") al_ptr = & cmd_line.align[1];
  else if (arg2.upcase() == "SEQ2") al_ptr = & cmd_line.align[2];
  else if (arg2.upcase() == "ALIGN") al_ptr = & cmd_line.align[0];
  else 
    {
      os << "SYNTAX ERROR: Usage PRINT ID <sequence register>\n";
      return 1;
    }
  // Open file, if necessary
  ofstream of;
  ostream * ofp;
  if (strlen(arg3.c_str()) < 1)
    ofp = & os;
  else
    {
      of.open(arg3.c_str(), ios::out);
      ofp = & of;
    }
  if (ofp->good())
    {
      *ofp << "\nPercent sequence identities:\n     ";
      uint i, j;
      for (i=0; i < al_ptr->n_seq(); ++i)
	{
	  ofp->width(4);
	  *ofp << i + 1;
	}
      *ofp << " Avg";
      *ofp << "\n";
      for (i=0; i < al_ptr->n_seq(); ++i) 
	{
          double average_numerator = 0.0000001;
          double average_denominator = 0.0000001;
	  ofp->width(4);
	  *ofp << i + 1 << ")";
	  for (j=0; j < al_ptr->n_seq(); ++j) 
	    {
	      ofp->width(4);
	      *ofp << (int) (100 * (*al_ptr)[i].identity((*al_ptr)[j]));
              if (i != j) {
                average_numerator += ((*al_ptr)[i].identity((*al_ptr)[j])) * (*al_ptr)[j].weight();
                average_denominator += (*al_ptr)[j].weight();
              }
	    }
	  ofp->width(4);
          *ofp << (int) (100 * average_numerator / average_denominator);
	  *ofp << "\n";
	}
    }
  else 
    {
      os << "Could not open file\n";
      of.close();
      return 1;
    }
  of.close();
  return 0;
}


// Print comparison matrix
int command_token::print_matrix(ostream & os) // PRINT MATRIX
{
  Mystring arg2 = get_one_token(rest_of_line);
  rest_of_line = advance_to_next_token(rest_of_line);
  Mystring arg3 = get_one_token(rest_of_line);
  rest_of_line = advance_to_next_token(rest_of_line);

  // Open file, if necessary
  ofstream of;
  ostream * ofp;
  if (strlen(arg2.c_str()) < 1)
    ofp = & os;
  else
    {
      of.open(arg2.c_str(), ios::out);
      ofp = & of;
    }
  if (ofp->good())
    {
      *ofp << cmd_line.cm;
    }
  else 
    {
      os << "Could not open file\n";
      of.close();
      return 1;
    }
  of.close();
  return 0;
}


int command_token::quit_program(ostream & os) // QUIT
{
  os << "Program terminated normally\n";
  myoshook(os);
  exit(0);
}


int command_token::random(ostream & os) // RANDOM
{
  static Array<MATRIX_TYPE> old_rscores;
  static double old_S;

  Mystring arg1 = get_one_token(rest_of_line);
  rest_of_line = advance_to_next_token(rest_of_line);

  // Run align so that current score is correct
  os << "Running ALIGN command...\n";
  if (cmd_line.parse_and_run2("ALIGN", os)) return 1;
  
  if (!cmd_line.sm.have_table) 
    {
      os << "Running TABULATE command...\n";
      if (cmd_line.parse_and_run2("TABULATE", os)) return 1;
    }

  uint nrand;
  uint report;
  double mean, stddev;
  uint i;

  if (!is_a_number(arg1.c_str())) nrand = 10;
  else nrand = atoi(arg1.c_str());
  if (nrand < 1) nrand = 10;
  
  if (nrand < 15) report = 1;
  else if (nrand < 150) report = 10;
  else if (nrand < 1500) report = 100;
  else if (nrand < 15000) report = 1000;
  else report = 10000;
  
  srand((int)cmd_line.sm.random_seed());
  mean = 0;
  // MATRIX_TYPE * rscores = new MATRIX_TYPE [nrand];
  Array<MATRIX_TYPE> rscores(nrand);
  for(i=0; i < nrand; ++i)
    {
      if ((!i) || (! ((i+1) % report)))
	os << "Calculating randomization number " << i+1 << "\n";
      myoshook(os);
      rscores[i] = cmd_line.sm.randomize();
      mean += rscores[i];
    }
  mean = mean/nrand;
  stddev = 0;
  double moment1, moment2, moment3, moment4;
  double m1;
  moment1 = moment2 = moment3 = moment4 = 0;
  for(i=0; i < nrand; ++i) {
    m1 = rscores[i] - mean;
    moment1 += m1;
    moment2 += m1 * m1;
    moment3 += m1 * m1 * m1;
    moment4 += m1 * m1 * m1 * m1;
  }
  moment2 = moment2 / (nrand - 1);
  // moment2 = moment2 / nrand;
  stddev = sqrt(moment2);
  moment1 = moment1 / nrand;
  moment3 = moment3 / nrand;
  moment3 = moment3 / (stddev * moment2); // skewness
  moment4 = moment4 / nrand;
  moment4 = moment4 / (moment2 * moment2) - 3.0; // kurtosis
  os << "Mean is " << mean;
  os << ", variance is " << moment2 << "\n";
  os << "Skewness is " << moment3;
  os << ", Kurtosis is " << moment4 << " (vs. 1.14,2.4 for EVD)" << endl;

  // Should we do delta statistics?
  double S = cmd_line.align[0].score_val;
  if (rscores.dim() == old_rscores.dim()) {
    Array<MATRIX_TYPE> delta_rscores = rscores;
    double delta_mean = 0;
    for (i = 0; i < nrand; ++i) {
      delta_rscores[i] = rscores[i] - old_rscores[i];
      delta_mean += delta_rscores[i];
    }
    delta_mean /= nrand;
    double delta_variance = 0;
    for (i = 0; i < nrand; ++i) {
      double diff = delta_rscores[i] - delta_mean;
      delta_variance += diff * diff;
    }
    delta_variance /= (nrand -1);
    os << "Delta statistics vs. previous run: " << endl;
    os << "mean = " << delta_mean << ", variance = ";
    os << delta_variance << ", zscore = ";
    double zscore;
    if (((S - old_S) - delta_mean) == 0) zscore = 0;
    else if (delta_variance == 0) zscore = 0;
    else zscore = ((S - old_S) - delta_mean) / sqrt(delta_variance);
    os << zscore;
    if (zscore > 0) os << " (BETTER) ";
    else os << " (WORSE) ";
    os << endl;
  }

  // Print out histogram
  // Figure out maximum and minimum score
  double max_score, min_score;
  max_score = min_score = rscores[0];
  for(i=0; i < nrand; ++i) {
    if (max_score < rscores[i]) max_score = rscores[i];
    if (min_score > rscores[i]) min_score = rscores[i];
  }
  uint n_bins = 15;
  int score_bins[n_bins];
  for (i=0; i<n_bins; ++i) score_bins[i] = 0;
  int largest_bin_size = 0;
  for(i=0; i < nrand; ++i) {
    uint my_bin = (uint) ((rscores[i] - min_score) * 
      (n_bins / (max_score - min_score)));
    if (my_bin >= n_bins) my_bin = n_bins - 1;
    score_bins[my_bin] = score_bins[my_bin] + 1;
    if (score_bins[my_bin] > largest_bin_size) 
      largest_bin_size = score_bins[my_bin];
  }
  double bin_scale = 60.0 / largest_bin_size;
  int j;
  cout << "Score histogram:" << endl;
  for (i=0; i < n_bins; ++i) {
    for (j=0; j < score_bins[i] * bin_scale; ++j)
      cout << "#"; // use hash marks for histogram
    cout << endl;
  }

  // Print final statistics
  os << "Alignment score of " << S << " is ";
  os << (S - mean) / stddev << " standard deviations";
  os << " above the mean random score\n";

  double euler_gamma = 0.57721566;
  double sqrt6 = 2.449489743;
  double lambda = PI / (sqrt6 * stddev);
  double mu = mean - (euler_gamma / lambda);
  cout << "E.V.D. lambda = " << lambda;
  cout << "     mu = " << mu << endl;
  double exponent = lambda * (S - mu);
  double p_value;
  os << "E.V.D. P-value = ";
  if (exponent < 10) // use precise expression
    p_value = 1.0 - exp(-exp(-exponent));
  else // use numerically stable approximation
    p_value = exp(-exponent);
  os << p_value;
  if (p_value > 0.05) os << " (not significant)";
  os << endl;
  
  old_rscores = rscores;
  old_S = S;
  // delete [] rscores;
  return 0;
}


int command_token::read_sequence_file(ostream & os) // READ
{
  // Parse two arguments after READ command
  Mystring arg1 = get_one_token(rest_of_line);
  rest_of_line = advance_to_next_token(rest_of_line);
  Mystring arg2 = get_one_token(rest_of_line);

  // open file
  ifstream inf;
  if (strlen(arg2.c_str()) < 1)
    {
      os << "!!! SYNTAX ERROR !!! (no filename)\n"; 
      return 1;
    }
  inf.open(arg2.c_str(), ios::in);
  if (!inf.good())
    {
      inf.close();
      os << "!!! Error: unable to open file ``" << arg2 << "''" << endl;
      return 1;
    }
	
  // Sequence/structure registers
  if (is_a_register_name(arg1.c_str()))
    {
      Alignment * al_ptr = NULL;
      pdboverlay * str_ptr = NULL;
      if (arg1.upcase() == "SEQ1") al_ptr = & cmd_line.align[1];
      else if (arg1.upcase() == "SEQ2") al_ptr = & cmd_line.align[2];
      else if (arg1.upcase() == "ALIGN") al_ptr = & cmd_line.align[0];
      else if (arg1.upcase() == "STRUCT1") 
	{str_ptr = & cmd_line.strct[1];al_ptr = & cmd_line.align[1];}
      else if (arg1.upcase() == "STRUCT2") 
	{str_ptr = & cmd_line.strct[2];al_ptr = & cmd_line.align[2];}
      else if (arg1.upcase() == "OVERLAY") 
	{str_ptr = & cmd_line.strct[0];al_ptr = & cmd_line.align[0];}
      else 
	{
	  os << "Internal error!!!, unknown register!!!\n";
	  os << name() << "\n";
	  inf.close();
	  return 1;
	}

      os << "Opening file ``" << arg2 << "''..." << endl;
      if (str_ptr == NULL) // sequence
	{
	  al_ptr->clear();
	  inf >> *al_ptr;
	  os << "Read in " << al_ptr->n_seq() << " sequence(s), ";
	  os << "each " << al_ptr->length() << " long." << endl;
	}
      else // coordinates
	{
	  pdboverlay empty;
	  *str_ptr = empty;
	  inf >> * str_ptr;
	  os << "Read in " << str_ptr->dim() << " model(s), ";
	  os << "with a total of " << str_ptr->n_atoms() << " atoms" <<endl;
	  *al_ptr = *str_ptr;
	  uint i;
          for (i=0; i<al_ptr->dim(); ++i)
	    if ((*al_ptr)[i].title().length() < 1)
	      (*al_ptr)[i].title() = arg2;
	  if (al_ptr == &cmd_line.align[0]) os << "Replaced ALIGN\n";
	  else if (al_ptr == &cmd_line.align[1])
	    {
	      os << "Replaced SEQ1\n";
	      cmd_line.sm.have_stable = false;
	      cmd_line.sm.have_table = false;
	      cmd_line.sm.have_seq1 = true;
	    }
	  else if (al_ptr == &cmd_line.align[2])
	    {
	      os << "Replaced SEQ2 and OVERLAY\n";
	      cmd_line.strct[0] = cmd_line.strct[2];
	      cmd_line.sm.have_stable = false;
	      cmd_line.sm.have_table = false;
	      cmd_line.sm.have_seq2 = true;
	      cmd_line.sm.have_overlay = true;
	    }
	}
      cmd_line.sm.have_table = false;
      if (al_ptr == & cmd_line.align[0]) cmd_line.sm.have_align = false;
      else if (al_ptr == &cmd_line.align[1]) 
	{
	  cmd_line.sm.have_seq1 = true;
	  cmd_line.sm.have_align = false;
	  cmd_line.sm.have_table = false;
	  cmd_line.sm.have_stable = false;
        }
      else if (al_ptr == &cmd_line.align[2]) 
	{
	  cmd_line.sm.have_seq2 = true;
	  cmd_line.sm.have_align = false;
	  cmd_line.sm.have_table = false;
	  cmd_line.sm.have_stable = false;
        }
    }

  // read a comparison matrix
  else if (arg1.upcase() == "MATRIX")
    {
      os << "Reading comparison matrix from file ``" << arg2 << "''"
	<< endl;
      inf >> cmd_line.cm; 
      cmd_line.sm.have_table = false;
    }

  // Don't know what the first argument is
  else
    {
      os << "SYNTAX ERROR: ``" << arg1 << "'' not understood" << endl;
      inf.close();
      return 1;
    }

  inf.close();
  return 0;
}

// **** NOTE, ROAD command has been skipped here ****

// this is for execution of SEQUOIA script files
int command_token::run_command(ostream & os) // RUN
{
  Mystring fname = get_one_token(rest_of_line);
  rest_of_line = advance_to_next_token(rest_of_line);
  if (strlen(fname.c_str()) < 1)
    {
      os << "!!! SYNTAX ERROR !!! (need a filename)\n"; 
      return 1;
    }
  else
    {
      //
      // in order for recursion to work here,
      // we need to have our own stream in this function
      // 
      // Therefore, we need to parse the filename and open it here too
      //
      ifstream fp(fname.c_str(), ios::in);
      os << "Opening file " << fname << "\n";
      myoshook(os);
      char buffer[2000];
      if (!fp.good())
	{os << "Could not open file " << fname << "\n"; return 1;}
      else while (fp.good())
	{
	  fp.getline(buffer, 1998);
	  if (cmd_line.parse_and_run2(buffer, os)) return 1;
	  myoshook(os);
	}
      return 0;
    }
  return 0;
}


int command_token::salign_sequences(ostream & os) // SALIGN
{
  if (!(cmd_line.align[1].n_seq() && cmd_line.align[2].n_seq()))
    {
      os << "Need two sequences to align!!!\n";
      os << "(Use READ command)\n";
      return 1;
    }
  if (!cmd_line.sm.have_stable)
    {
      os << "Running STABULATE command...\n";
      myoshook(os);
      if (cmd_line.parse_and_run2("STABULATE", os)) return 1;
    }
  os << "Aligning sequences...\n";
  myoshook(os);
  cmd_line.align[0] = cmd_line.sm.dyn_prog();
  os << cmd_line.align[0].newgaps() << " gaps added to make alignment\n";
  os << "Alignment score is " << cmd_line.align[0].score_val << "\n";
  cmd_line.sm.have_align = true;
  return 0;
}


int command_token::set_variable(ostream & os) // SET
{
  // Parse two arguments after SET command
  Mystring arg1 = get_one_token(rest_of_line);
  rest_of_line = advance_to_next_token(rest_of_line);
  Mystring arg2 = get_one_token(rest_of_line);

  // List all variables
  if (arg1.length() < 1) list_variables(os); // no arguments, list

  // Copy register contents
  else if (is_a_register_name(arg1.c_str()))
    {
      Alignment * al_ptr = NULL;
      pdboverlay * str_ptr = NULL;
      if (arg1.upcase() == "SEQ1") al_ptr = & cmd_line.align[1];
      else if (arg1.upcase() == "SEQ2") al_ptr = & cmd_line.align[2];
      else if (arg1.upcase() == "ALIGN") al_ptr = & cmd_line.align[0];
      else if (arg1.upcase() == "STRUCT1") str_ptr = & cmd_line.strct[1];
      else if (arg1.upcase() == "STRUCT2") str_ptr = & cmd_line.strct[2];
      else if (arg1.upcase() == "OVERLAY") str_ptr = & cmd_line.strct[0];
      else 
	{
	  os << "Internal error!!!, unknown register!!!\n";
	  os << name() << "\n";
	  return 1;
	}
      // only one register named?
      if (!is_a_register_name(arg2.c_str()))
	{
	  if (al_ptr != NULL) al_ptr->summarize(os);
	  else {} // need summary output for overlays too
	  return 0;
	}
      else // COPY register
	{
	  Alignment * al_ptr2 = NULL;
	  pdboverlay * str_ptr2 = NULL;
	  if (arg2.upcase() == "SEQ1") al_ptr2 = & cmd_line.align[1];
	  else if (arg2.upcase() == "SEQ2") al_ptr2 = & cmd_line.align[2];
	  else if (arg2.upcase() == "ALIGN") al_ptr2 = & cmd_line.align[0];
	  else if (arg2.upcase() == "STRUCT1") str_ptr2 = & cmd_line.strct[1];
	  else if (arg2.upcase() == "STRUCT2") str_ptr2 = & cmd_line.strct[2];
	  else if (arg2.upcase() == "OVERLAY") str_ptr2 = & cmd_line.strct[0];
	  else 
	    {
	      os << "Internal error!!!, unknown register!!!\n";
	      os << name() << "\n";
	      return 1;
	    }  
	  if (str_ptr != NULL)
	    {
	      // We can't read a sequence into a structure (with current
	      // technology!)
	      if (al_ptr2 != NULL)
		{
		  os << "Cannot make a structure from a sequence!!!\n";
		  os << "Sorry, maybe next version... :)\n";
		  return 1;
		}
	      else *str_ptr = *str_ptr2;
	    }
	  else
	    {
	      if (al_ptr2 != NULL) *al_ptr = *al_ptr2;
	      else *al_ptr = *str_ptr2;
	    }
	  os << "Copied " << arg2 << " to " << arg1 << endl;
	  return 0;
	}
    }

  // Assign a Real value to a variable
  else if (is_a_number(arg2.c_str()))
    {
      double value = atof(arg2.c_str());
      assign_variable(arg1.c_str(), value);
    }

  // Print out value of variable
  else if (variable_exists(arg1.c_str()))
    {
      sequoia_variable * var = sq_var_ptr(arg1.c_str());
      os << var->name() << " = " << var->value() << endl;
    }

  // could not understand first argument
  else 
    {
      os << "SYNTAX ERROR: usage SET <variable> <value>" << endl;
      return 1;
    }
  return 0;
}


int command_token::split_sequences(ostream & os) // SPLIT
{
  Mystring arg1 = get_one_token(rest_of_line);
  rest_of_line = advance_to_next_token(rest_of_line);

  if (cmd_line.align[0].n_seq() < 2)
    {
      os << "At least two sequences required to split!!\n";
      return 1;
    }
  uint extract_seq = 1;
  if (is_a_number(arg1.c_str()))
    {
      extract_seq = atoi(arg1.c_str());
    }
  if ((extract_seq < 1) || (extract_seq > cmd_line.align[0].n_seq()))
    {
      os << "There is no sequence " << extract_seq << " !!\n";
      return 1;
    }
  --extract_seq; // to make it index the C-like matrix
  uint j;
  Alignment empty;
  os << "Placing sequence " << extract_seq+1 << " in SEQ1...\n";
  myoshook(os);
  cmd_line.align[1] = empty;
  cmd_line.align[1] += cmd_line.align[0][extract_seq].remove_gaps();
  os << "Placing others in SEQ2...\n";
  cmd_line.align[2] = empty;
  for (j=0; j < cmd_line.align[0].n_seq(); ++j)
    if (j != extract_seq) cmd_line.align[2] += cmd_line.align[0][j];
  cmd_line.align[2] = cmd_line.align[2].remove_gaps();
  cmd_line.sm.have_seq1 = true;
  cmd_line.sm.have_seq2 = true;
  cmd_line.sm.have_table = false;
  os << "reTABULATEing scores...\n";
  if (cmd_line.parse_and_run2("TABULATE", os)) return 1;
  return 0;
}


int command_token::stabulate(ostream & os) // STABULATE
{
  if (!(cmd_line.strct[1].ca().dim() && cmd_line.strct[0].ca().dim()))
    {
      os << "Need some alpha carbons to compare!\n";
      return 1;
    }
  cmd_line.strct[1].update_refatoms();
  cmd_line.strct[0].update_refatoms();
  cmd_line.strct[1].update_orientations();
  cmd_line.strct[0].update_orientations();
  cmd_line.sm.stabulate(cmd_line.strct[1], cmd_line.strct[0]);
  os << "Mean residue alignment score = " << cmd_line.sm.mean();
  cmd_line.sm.sigma = cmd_line.sm.std_dev();
  os << " with a standard deviation of " << cmd_line.sm.sigma << "\n";
  cmd_line.sm.have_table = false;
  cmd_line.sm.have_stable = true;
  cmd_line.sm.have_align = false;
  return 0;
}


int command_token::superpose(ostream & os) // SUPERPOSE
{
  os << "Superposing STRUCT2 onto STRUCT1\n";
  if (cmd_line.strct[1].n_atoms() < 4)
    {
      os << "STRUCT1 needs at least 4 atoms!\n";
      return 1;
    }
  if (cmd_line.strct[2].n_atoms() < 4)
    {
      os << "STRUCT2 needs at least 4 atoms!\n";
      return 1;
    }
  if (!cmd_line.sm.have_align)
    {
      os << "ALIGN register must be current for overlay to work!\n";
      os << "Running ALIGN, to see if this helps...\n";
      if (cmd_line.parse_and_run2("ALIGN", os)) return 1;
      else if (cmd_line.align[0].dim() < 2) return 1;
    }
  if (cmd_line.align[0].dim() < 2)
    {
      os << "ALIGN must have at least 2 sequences to use for overlay!\n";
      os << "Running ALIGN, to see if this helps...\n";
      if (cmd_line.parse_and_run2("ALIGN", os)) return 1;
      else if (cmd_line.align[0].dim() < 2) return 1;
    }
  unsigned int n1 = cmd_line.align[1].dim();
  unsigned int n2 = cmd_line.align[2].dim();
  if (cmd_line.align[0].dim() < (n1 + n2))
    {
      os << "ALIGN must have at least as many sequences as SEQ1+SEQ2!\n";
      return 1;
    }
  if (cmd_line.align[0][0].n_cap_aligned(cmd_line.align[0][n1]) < 4)
    {
      os << "ALIGN must have at least 4 capitalized equivalences!\n";
      return 1;
    }
  RMat xform;
  os.setf(ios::showpoint);
  os.precision(4);
  os << "Starting RMS deviation is " 
    << cmd_line.align[0].rms(cmd_line.strct[1],cmd_line.strct[0])
      << " for " 
	<< cmd_line.align[0][0].n_cap_aligned(cmd_line.align[0][n1]) 
	  << " equivalent matches.\n";
  os << "Starting RMS deviation is " 
    << cmd_line.align[0].upcase().rms(cmd_line.strct[1],cmd_line.strct[0])
      << " for "
	<< cmd_line.align[0][0].n_aligned(cmd_line.align[0][n1]) 
	  << " total matches.\n";
  os << "Transforming STRUCT2...\n";
  os << "calculating rotation matrix...\n";
  myoshook(os);
  xform = overlay(cmd_line.strct[1],cmd_line.strct[2],cmd_line.align[0]); 
  os << "Transformation matrix is:\n" << xform << endl;
  os << "applying matrix...\n";
  cmd_line.strct[0] = cmd_line.strct[2].transform(xform);
  cmd_line.sm.have_stable = false;
  os << "New RMS deviation is " 
    << cmd_line.align[0].rms(cmd_line.strct[1],cmd_line.strct[0]) 
      << " for " 
	<< cmd_line.align[0][0].n_cap_aligned(cmd_line.align[0][n1]) 
	  << " equivalent matches.\n";
  os << "New RMS deviation is " 
    << cmd_line.align[0].upcase().rms(cmd_line.strct[1],cmd_line.strct[0]) 
      << " for " 
	<< cmd_line.align[0][0].n_aligned(cmd_line.align[0][1]) 
	  << " total matches.\n";
  return 0;
}


int command_token::system_command(ostream & os) // SYSTEM
{
  return system(rest_of_line);
}


int command_token::tabulate(ostream & os) // TABULATE
{
  if (!(cmd_line.align[1].n_seq() && cmd_line.align[2].n_seq()))
    {
      os << "Need two sequences/alignments to tabulate!!!\n";
      os << "(Use READ command)\n";
      return 1;
    }
  cmd_line.sm.tabulate(cmd_line.align[1], cmd_line.align[2], cmd_line.cm);
  os << "Mean residue alignment score = " << cmd_line.sm.mean();
  cmd_line.sm.sigma = cmd_line.sm.std_dev();
  os << " with a standard deviation of " << cmd_line.sm.sigma << "\n";
  cmd_line.sm.have_table = true;
  cmd_line.sm.have_stable = false;
  cmd_line.sm.have_align = false;
  return 0;
}


int command_token::weight(ostream & os) // WEIGHT
{
  Mystring arg1 = get_one_token(rest_of_line);
  rest_of_line = advance_to_next_token(rest_of_line);
  Mystring arg2 = get_one_token(rest_of_line);
  rest_of_line = advance_to_next_token(rest_of_line);
  Mystring arg3 = get_one_token(rest_of_line);
  rest_of_line = advance_to_next_token(rest_of_line);

  // If more arguments are given, manually weight one sequence
  // WEIGHT <register> <index> <value>
  // weight seq1 1 3.0
  if (is_a_register_name(arg1.c_str()))
  {
    int sequence_index = atoi(arg2.c_str());
    double new_weight = atof(arg3.c_str());
    Alignment * al_ptr = NULL;
    pdboverlay * str_ptr = NULL;
    if (arg1.upcase() == "SEQ1") al_ptr = & cmd_line.align[1];
    else if (arg1.upcase() == "SEQ2") al_ptr = & cmd_line.align[2];
    else if (arg1.upcase() == "ALIGN") al_ptr = & cmd_line.align[0];
    else if (arg1.upcase() == "STRUCT1") str_ptr = & cmd_line.strct[1];
    else if (arg1.upcase() == "STRUCT2") str_ptr = & cmd_line.strct[2];
    else if (arg1.upcase() == "OVERLAY") str_ptr = & cmd_line.strct[0];
    if (al_ptr != NULL) {
      Alignment & al = *al_ptr;
      al[sequence_index - 1].weight() = new_weight;
      os << "The weight of sequence " << sequence_index << " of register " <<
            arg1 << " has been set to " << new_weight << endl;
    }
  }
  else {
    cmd_line.align[0].set_weights();
    cmd_line.align[1].set_weights();
    cmd_line.align[2].set_weights();
    cmd_line.sm.have_table = false; // affects tabulation
    os << "SEQ1, SEQ2, and ALIGN have been Bruns-weighted\n";
    cmd_line.sm.have_table = false;
  }
  return 0;
}

int command_token::write_sequence_file(ostream & os) // WRITE
{
  Mystring arg1 = get_one_token(rest_of_line);
  rest_of_line = advance_to_next_token(rest_of_line);
  Mystring arg2 = get_one_token(rest_of_line);
  rest_of_line = advance_to_next_token(rest_of_line);

  // Print contents of register
  if (is_a_register_name(arg1.c_str()))
    {
      Alignment * al_ptr = NULL;
      pdboverlay * str_ptr = NULL;
      if (arg1.upcase() == "SEQ1") al_ptr = & cmd_line.align[1];
      else if (arg1.upcase() == "SEQ2") al_ptr = & cmd_line.align[2];
      else if (arg1.upcase() == "ALIGN") al_ptr = & cmd_line.align[0];
      else if (arg1.upcase() == "STRUCT1") str_ptr = & cmd_line.strct[1];
      else if (arg1.upcase() == "STRUCT2") str_ptr = & cmd_line.strct[2];
      else if (arg1.upcase() == "OVERLAY") str_ptr = & cmd_line.strct[0];
      else 
	{
	  os << "Internal error!!!, unknown register!!!\n";
	  os << arg1 << "\n";
	  return 1;
	}

      // Open file, if necessary
      ofstream of;
      ostream * ofp;
      if (strlen(arg2.c_str()) < 1)
	ofp = & os;
      else
	{
	  of.open(arg2.c_str(), ios::out);
	  ofp = & of;
	}
      if (ofp->good())
	{
	  os << "Opening file ...\n";
	  if (str_ptr != NULL)
	    {
	      if (str_ptr == &cmd_line.strct[0]) // Is the overlay, 
		// needs special handling
		{
		  pdboverlay special(0);
		  uint str;
		  for (str = 0; str < cmd_line.strct[1].dim(); ++str)
		    {
		      special += cmd_line.strct[1][str];
		      special[str].sequence() = cmd_line.align[0][str];
		    }
		  for (str = cmd_line.strct[1].dim(); 
		       str < (cmd_line.strct[0].dim() + cmd_line.strct[1].dim()); 
		       ++str)
		    {
		      special += cmd_line.strct[0]
			[str - cmd_line.strct[1].dim()]; 
		      special[str].sequence() = cmd_line.align[0][str];
		    }
		  *ofp << special;
		}
	      else
		{
		  *ofp << * str_ptr;
		  os << "Wrote " << str_ptr->n_atoms() << " atom(s)\n";
		}
	    }
	  else
	    {
	      *ofp << *al_ptr;
	      os << "Wrote " << al_ptr->n_seq() << " sequence(s)\n";
	    }
	}
      else 
	{
	  os << "Could not open file\n";
	  of.close();
	  return 1;
	}
      of.close();
    }

  return 0;
}

