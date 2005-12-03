#include "parse.hxx"
#include "oshook.hxx"
#include <iomanip>

// $Id: parse.cxx,v 1.4 2002/06/28 17:21:54 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/parse.cxx,v 1.4 2002/06/28 17:21:54 bruns Exp $
// $Log: parse.cxx,v $
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

cli::cli(const Mystring & s)
{
  prompt_string = s;
  Comparison_matrix temp(BLOSUM62);
  cm = temp;
  sm.have_table = false;
  sm.have_stable = false;
  sm.have_path = false;
  sm.have_align = false;
  sm.have_seq1 = false;
  sm.have_seq2 = false;
  sm.have_strct1 = false;
  sm.have_strct2 = false;
  sm.have_overlay = false;
  sm.useangle() = true; // this is where default value is set
}

const Mystring & cli::prompt(void) const
{
  return prompt_string;
}

// this is the meat of the command line parser
// 
// When adding things in, remember to
//  1) provide a return value, zero if no error
//  2) check the return value for any recursive calls to parse_and_run
//  3) update values of sm.have_align, have_stable, have_seq1 etc.
// 
int cli::parse_and_run(const char * buffline, ostream & os)
{
  char buffer[2000]; // limit on input line
  Array<Mystring> tokens;
  Mystring command;
  Mystring arg;
  uint next_tok;
  ifstream fp;
  ofstream ofp;
  Alignment * al_ptr;
  pdboverlay * str_ptr;
  static Comparison_matrix b62(BLOSUM62);
  static int pretty_length = 50;

  Mystring buf_string = buffline;
  tokens = buf_string.split();

  if (tokens.dim()) {command = tokens[0]; next_tok = 1;}
  else return 1; // command line is blank?

  // escape a system command
  if (command.upcase().at(0,1) == "!")
    {
      return system(buf_string.after('!').c_str());
    }

  // invoke a script file
  else if (command.upcase().at(0,1) == "@")
    {
      arg = command.after((uint)0);
      if (arg != "")
	{
	  strcpy(buffer, arg.c_str());
	  os << "Opening file " << arg << "\n";
          myoshook(os);
	  fp.close();
	  fp.open(buffer);
	  if (!fp.good())
	    {os << "Could not open file " << arg << "\n"; return 1;}
	  else while (fp.good())
	    {
	      fp.getline(buffer, 1998);
	      if (this->parse_and_run(buffer, os)) return 1;
              myoshook(os);
	    }
	}
      else {os << "Could not open file\n"; return 1;}
      return 0;
    }

  // Many ways to leave the program
  else if ( ( command.upcase() == "BYE") ||
	   (command.upcase() == "END") ||
	   (command.upcase() == "EXIT") ||
	   (command.upcase() == "KILL") ||
	   (command.upcase() == "Q") ||
	   (command.upcase() == "QUIT") ||
	   (command.upcase() == "STOP") )
    {
      os << "Program terminated by user\n";
      myoshook(os);
      exit(0);
    }

  // Comments
  else if ( (upcase(command.at(0,1)) == "!") ||
	   (upcase(command.at(0,1)) == "#") ||
	   (upcase(command.at(0,1)) == "{") ||
	   (upcase(command.at(0,3)) == "REM") ||
	   (upcase(command.at(0,1)) == "/") )
    return 0;

  // Angular orientation cutoff for structural alignments
  else if (upcase(command.at(0,5)) == "ACUTO")
    {
      if (tokens.dim() > next_tok)
	sm.acutoff() = atof(tokens[next_tok].c_str());
      os << "Orientation cutoff for calculation is now " ;
      os << sm.acutoff() << " degrees\n";
      sm.have_stable = false; // affects stabulate result
      return 0;
    }

  else if (upcase(command.at(0,5)) == "ALIGN")
    {
      if (!(align[1].n_seq() && align[2].n_seq()))
        {
          os << "Need two sequences/alignments to align!!!\n";
          os << "(Use READ command)\n";
          return 1;
        }
      if (!sm.have_table) 
	{
	  os << "Running TABULATE command...\n";
          myoshook(os);
	  strcpy(buffer, "TABULATE");
	  if (this->parse_and_run(buffer, os)) return 1;
	}
      os << "Aligning sequences...\n";
      myoshook(os);
      align[0] = sm.dyn_prog();
      os << align[0].newgaps() << " gaps added to make alignment\n";
      os << "Alignment score is " << align[0].score_val << "\n";
      sm.have_align = true;
      return 0;
    }

  // generate a consensus sequence
  else if (upcase(command.at(0,6)) == "CONSEN")
    {
      double cutoff = 1.2; // default sigma cutoff
      if (tokens.dim() > next_tok)
	{
	  cutoff = atof(tokens[next_tok].c_str());
	  os << "Cutoff is " << cutoff << " sigma\n";
	}
      else os << "Using default cutoff of " << cutoff << " sigma\n";
      Sequence cons = consensus(align[0], cm, cutoff);
      cons.title() = "Consensus";
      align[2] = cons.remove_gaps();
      os << "Consensus of ALIGN has been placed in SEQ2\n";
      sm.have_seq2 = true;
      return 0;
    }

  else if (upcase(command.at(0,4)) == "COPY")
    {
      al_ptr = NULL;
      str_ptr = NULL;
      // figure out first argument
      if (tokens.dim() <= next_tok)
	{al_ptr = & align[0]; --next_tok;}
      else if (upcase(tokens[next_tok]) == "SEQ1")
	al_ptr = & align[1];
      else if (upcase(tokens[next_tok]) == "SEQ2")
	al_ptr = & align[2];
      else if (upcase(tokens[next_tok]).at(0,5) == "ALIGN")
	al_ptr = & align[0];
      else if (upcase(tokens[next_tok]) == "STRUCT1")
	str_ptr = & strct[1];
      else if (upcase(tokens[next_tok]) == "STRUCT2")
	str_ptr = & strct[2];
      else if (upcase(tokens[next_tok]).at(0,4) == "OVER")
	str_ptr = & strct[0];
      else
	{al_ptr = & align[0]; --next_tok;}
      ++next_tok;

      // just one argument - put into seq1
      if (tokens.dim() <= next_tok)
	{
	  if (al_ptr == & align[0])
	    {align[1] = *al_ptr; os << "copied ALIGN to SEQ1\n";}
	  else if (al_ptr == & align[1])
	    {align[1] = *al_ptr; os << "copied SEQ1 to SEQ1\n";}
	  else if (al_ptr == & align[2])
	    {align[1] = *al_ptr; os << "copied SEQ2 to SEQ1\n";}
	  else if (str_ptr == & strct[0])
	    {align[1] = *str_ptr; os << "copied OVERLAY to SEQ1\n";}
	  else if (str_ptr == & strct[1])
	    {align[1] = *str_ptr; os << "copied STRUCT1 to SEQ1\n";}
	  else if (str_ptr == & strct[2])
	    {align[1] = *str_ptr; os << "copied STRUCT2 to SEQ1\n";}
	  sm.have_table = false;
	  sm.have_seq1 = true;
	}
      // argument 2 = seq2
      else if (upcase(tokens[next_tok]) == "SEQ2")
	{
	  if (al_ptr == & align[0])
	    {align[2] = *al_ptr; os << "copied ALIGN to SEQ2\n";}
	  else if (al_ptr == & align[1])
	    {align[2] = *al_ptr; os << "copied SEQ1 to SEQ2\n";}
	  else if (al_ptr == & align[2])
	    {align[2] = *al_ptr; os << "copied SEQ2 to SEQ2\n";}
	  else if (str_ptr == & strct[0])
	    {align[2] = *str_ptr; os << "copied OVERLAY to SEQ2\n";}
	  else if (str_ptr == & strct[1])
	    {align[2] = *str_ptr; os << "copied STRUCT1 to SEQ2\n";}
	  else if (str_ptr == & strct[2])
	    {align[2] = *str_ptr; os << "copied STRUCT2 to SEQ2\n";}
	  sm.have_table = false;
	  sm.have_seq2 = true;
	}
      // argument 2 = align
      else if (upcase(tokens[next_tok]).at(0,5) == "ALIGN")
	{
	  if (al_ptr == & align[0])
	    {align[0] = *al_ptr; os << "copied ALIGN to ALIGN\n";}
	  else if (al_ptr == & align[1])
	    {align[0] = *al_ptr; os << "copied SEQ1 to ALIGN\n";}
	  else if (al_ptr == & align[2])
	    {align[0] = *al_ptr; os << "copied SEQ2 to ALIGN\n";}
	  else if (str_ptr == & strct[0])
	    {align[0] = *str_ptr; os << "copied OVERLAY to ALIGN\n";}
	  else if (str_ptr == & strct[1])
	    {align[0] = *str_ptr; os << "copied STRUCT1 to ALIGN\n";}
	  else if (str_ptr == & strct[2])
	    {align[0] = *str_ptr; os << "copied STRUCT2 to ALIGN\n";}
	  sm.have_align = true;
	}
      // argument 2 should be seq1
      else
	{
	  if (al_ptr == & align[0])
	    {align[1] = *al_ptr; os << "copied ALIGN to SEQ1\n";}
	  else if (al_ptr == & align[1])
	    {align[1] = *al_ptr; os << "copied SEQ1 to SEQ1\n";}
	  else if (al_ptr == & align[2])
	    {align[1] = *al_ptr; os << "copied SEQ2 to SEQ1\n";}
	  else if (str_ptr == & strct[0])
	    {align[1] = *str_ptr; os << "copied OVERLAY to SEQ1\n";}
	  else if (str_ptr == & strct[1])
	    {align[1] = *str_ptr; os << "copied STRUCT1 to SEQ1\n";}
	  else if (str_ptr == & strct[2])
	    {align[1] = *str_ptr; os << "copied STRUCT2 to SEQ1\n";}
	  sm.have_table = false;
	  sm.have_seq1 = true;
	}
      return 0;
    }

  else if (upcase(command.at(0,4)) == "CUTO")
    {
      if (tokens.dim() > next_tok)
	sm.dcutoff() = atof(tokens[next_tok].c_str());
      os << "Atomic distance cutoff for calculation is now " ;
      os << sm.dcutoff() << " Angstroms\n";
      sm.have_stable = false; // affects stabulate result
      return 0;
    }

  else if (upcase(command.at(0,4)) == "EPEN")
    {
      if (tokens.dim() > next_tok)
	sm.epen() = atof(tokens[next_tok].c_str());
      os << "Gap extension penalty is set to " << sm.epen();
      os << " standard deviations.\n";
      sm.have_align = false; // affects alignment
      return 0;
    }

  // update capitalized residues based upon strct[0]:strct[1] alignment
  else if (upcase(command.at(0,5)) == "EQUIV")
    {
      os << "Reassigning equivalences...\n";
      myoshook(os);
      uint i;
      i = align[0][0].n_cap_aligned(align[0][1]);
      os << "Total starting equivalences = " << i << "\n";
      if (!sm.have_stable) if (parse_and_run("STABULATE",os)) return 1;
      align[0] = struct_equiv(strct[1],strct[0],align[0],sm.dcutoff(),(int)sm.runlength());
      i = align[0][0].n_cap_aligned(align[0][1]);
      os << "Total new equivalences = " << i << "\n";
      return 0;
    }


  else if (upcase(command.at(0,4)) == "GPEN")
    {
      if (tokens.dim() > next_tok)
	sm.gpen() = atof(tokens[next_tok].c_str());
      os << "Gap penalty is set to " << sm.gpen();
      os << " standard deviations.\n";
      sm.have_align = false; // affects alignment
      return 0;
    }


  else if (upcase(command.at(0,3)) == "LIS")
    {
      parse_and_run("GPEN",os);
      parse_and_run("EPEN",os);
      parse_and_run("CUTOFF",os);
      parse_and_run("ACUTOFF",os);
      parse_and_run("RUNLENGTH",os);      
      parse_and_run("ROAD",os); // OUCH
      parse_and_run("SUBOPTIMAL",os);
      if (sm.useangle()) os << "Orientation WILL be used in coordinate alignment\n";
      else os << "Orientation will NOT be used in coordinate alignment\n";
      os << "SEQ1: \n";
      align[1].summarize(os);
      os << "SEQ2: \n";
      align[2].summarize(os);
      return 0;
    }


  else if (upcase(command.at(0,3)) == "MAT")
    {
      if (tokens.dim() <= next_tok)
	{
	  os << "Current scoring matrix : \n";
	  os << cm;
	  return 0;
	}
      else if (upcase(tokens[next_tok]) == "BLOSUM62")
	{
	  cm = b62; 
	  sm.have_table = false;
	  return 0;
	}
      else
	{
	  fp.open(tokens[next_tok].c_str());
	  if (fp.good()) 
	    {
	      fp >> cm; 
	      sm.have_table = false;
	      return 0;
	    }
	  else 
	    {
	      os << "Could not open file " << tokens[next_tok] << "\n";
	      return 1;
	    }
	}
      return 0;
    }


  else if (upcase(command.at(0,4)) == "OPTI")
    {
      if (align[0].n_seq() < 2)
	{
	  os << "At least two sequences required to optimize!!\n";
	  return 1;
	}
      int i;
      for (i = align[0].n_seq() - 1; i >= 0; --i)
        {
	  char cmd[50];
	  sprintf(cmd, "SPLIT %d", align[0].n_seq());
	  os << "SPLITing out sequence " << i+1 << ":\n";
          myoshook(os);
	  if (this->parse_and_run(cmd, os)) return 1;
	  os << "reALIGNing:\n";
          myoshook(os);
	  if (this->parse_and_run("ALIGN", os)) return 1;
        }
      return 0;
    }


  // structure based overlay
  else if (upcase(command.at(0,4)) == "OVER")
    {
      Alignment old_align;
      // overlay struct 2 onto struct 1
      if (parse_and_run("SUPERPOSE",os)) return 1;
      myoshook(os);
      uint cycles = 0;
      while ((old_align != align[0]) && (cycles < 10))
	{
	  old_align = align[0];
	  if (parse_and_run("EQUIVALENT",os)) return 1;
          myoshook(os);
	  if (parse_and_run("SUPERPOSE",os)) return 1;
          myoshook(os);
	  ++ cycles;
	}
      if (old_align != align[0]) 
	os << "Failed to converge in " << cycles << " cycles\n";
      else os << "Equivalences converged in " << cycles << " cycles\n";
      os << "Realigning based upon final overlay\n";
      if (parse_and_run("SALIGN",os)) return 1;
      if (parse_and_run("EQUIVALENT",os)) return 1;
      return 0;
    }

  else if (upcase(command.at(0,4)) == "PLEN")
    {
      if (tokens.dim() > next_tok)
	{
	  pretty_length = atoi(tokens[next_tok].c_str());
	}
      os << "Line length for PRINT output is now set to ";
      os << pretty_length << " characters\n";
    }

  // print [seq1|seq2|align] [id] [filename]
  else if (upcase(command.at(0,5)) == "PRINT")
    {
      al_ptr = NULL;
      str_ptr = NULL;
      Mystring job = "mystery";
      ostream * destination = & os;
      // No arguments?, just print ALIGN then
      if (tokens.dim() <= next_tok)
	{al_ptr = & align[0]; job = "pretty"; destination = & os;}
      // Figure out first argument
      else{
	if (upcase(tokens[next_tok]) == "SEQ1")
	  {al_ptr = & align[1]; ++ next_tok;}
	else if (upcase(tokens[next_tok]) == "SEQ2")
	  {al_ptr = & align[2]; ++ next_tok;}
	else if (upcase(tokens[next_tok]) == "ALIGN")
	  {al_ptr = & align[0]; ++ next_tok;}
	else if (upcase(tokens[next_tok]) == "STRUCT1")
	  {str_ptr = & strct[1]; ++ next_tok;}
	else if (upcase(tokens[next_tok]) == "STRUCT2")
	  {str_ptr = & strct[2]; ++ next_tok;}
	else if (upcase(tokens[next_tok]) == "OVERLAY")
	  {str_ptr = & strct[0]; ++ next_tok;}
	else al_ptr = & align[0];
	if (tokens.dim() > next_tok)
	  {
	    if (al_ptr  && (upcase(tokens[next_tok]).at(0,2) == "ID"))
	      {
		job = "id";
		++ next_tok;
		// Is there another argument?
		if (tokens.dim() > next_tok)
		  {
		    ofp.open(tokens[next_tok].c_str());
		    if (ofp.good()) 
		      {
			destination = & ofp;
			os << "Opening file ``" << tokens[next_tok] << "''\n";
		      }
		    ++next_tok;
		  }
	      }
	    else
	      {
		job = "pretty";
		ofp.open(tokens[next_tok].c_str());
		if (ofp.good())
		  {
		    destination = & ofp;
		    os << "Opening file ``" << tokens[next_tok] << "''\n";
		  }
		++next_tok;
	      }
	  }
	else job = "pretty";
      }
      if (job == "pretty") 
	{
	  if (al_ptr) al_ptr->pretty(*destination, pretty_length);
	  // If it's a coordinate file
	  else 
	    {
	      if (str_ptr == &strct[0]) // Is the overlay, 
		                        // needs special handling
		{
		  pdboverlay special(0);
		  uint str;
		  for (str = 0; str < strct[1].dim(); ++str)
		    {
		      special += strct[1][str];
		      special[str].sequence() = align[0][str];
		    }
		  for (str = strct[1].dim(); 
		       str < (strct[1].dim() + strct[0].dim()); 
		       ++str)
		    {
		      special += strct[0][str - strct[1].dim()];
		      special[str].sequence() = align[0][str];
		    }
		  *destination << special;
		}
	      else
		{
		  *destination << * str_ptr;
		  os << "Wrote " << str_ptr->n_atoms() << " atom(s)\n";
		}
	    }
	}
      // (else) Print out identity table
      else
	{
	  *destination << "\nPercent sequence identities:\n     ";
	  uint i, j;
	  for (i=0; i < al_ptr->n_seq(); ++i)
	    {
	      destination->width(4);
	      *destination << i + 1;
	    }
	  *destination << "\n";
	  for (i=0; i < al_ptr->n_seq(); ++i) 
	    {
	      destination->width(4);
	      *destination << i + 1 << ")";
	      for (j=0; j < al_ptr->n_seq(); ++j) 
		{
		  destination->width(4);
		  *destination << (int) (100 * (*al_ptr)[i].identity((*al_ptr)[j]));
		}
	      *destination << "\n";
	    }
	}
      return 0;
    }
  
  else if (upcase(command.at(0,4)) == "RAND")
    {
      if (!sm.have_table) 
	{
	  os << "Running TABULATE command...\n";
	  strcpy(buffer, "TABULATE");
	  if (this->parse_and_run(buffer, os)) return 1;
        }
      uint nrand;
      uint report;
      double mean, stddev;
      uint i;
      if (tokens.dim() <= next_tok) nrand = 10;
      else nrand = atoi(tokens[next_tok].c_str());
      if (nrand < 1) nrand = 10;

      if (nrand < 15) report = 1;
      else if (nrand < 150) report = 10;
      else if (nrand < 1500) report = 100;
      else if (nrand < 15000) report = 1000;
      else report = 10000;

      srand((int)sm.random_seed());
      mean = 0;
      MATRIX_TYPE * rscores = new MATRIX_TYPE [nrand];
      for(i=0; i < nrand; ++i)
	{
	  if ((!i) || (! ((i+1) % report)))
	    os << "Calculating randomization number " << i+1 << "\n";
	  myoshook(os);
	  rscores[i] = sm.randomize();
	  mean += rscores[i];
	}
      mean = mean/nrand;
      stddev = 0;
      for(i=0; i < nrand; ++i)
	stddev += (mean - rscores[i]) * (mean - rscores[i]);
      stddev = stddev / (nrand - 1);
      stddev = sqrt(stddev);
      os << "Mean is " << mean << "\n";
      os << "Standard deviation is " << stddev << "\n";
      os << "Alignment score of " << align[0].score_val << " is ";
      os << (align[0].score_val - mean) / stddev << " standard deviations";
      os << " above the mean random score\n";
      delete [] rscores;
      return 0;
    }

  else if (upcase(command.at(0,4)) == "READ")
    {
      str_ptr = NULL;
      al_ptr = NULL;
      if (tokens.dim() <= next_tok)
	al_ptr = & align[1];
      else
	// figure out first argument
	{
	  if (upcase(tokens[next_tok]) == "STRUCT1")
	    {str_ptr = & strct[1]; al_ptr = & align[1];}
	  else if (upcase(tokens[next_tok]) == "STRUCT2")
	    {str_ptr = & strct[2]; al_ptr = & align[2];}
	  else if (upcase(tokens[next_tok]).at(0,4) == "OVER")
	    {str_ptr = & strct[0]; al_ptr = & align[0];}
	  else if (upcase(tokens[next_tok]) == "SEQ1")
	    al_ptr = & align[1];
	  else if (upcase(tokens[next_tok]) == "SEQ2")
	    al_ptr = & align[2];
	  else if (upcase(tokens[next_tok]).at(0,5) == "ALIGN")
	    al_ptr = & align[0];
	  else {al_ptr = & align[1]; --next_tok;}
	  ++ next_tok;
	}
      if (tokens.dim() > next_tok)
	{
	  fp.open(tokens[next_tok].c_str());
	  if (fp.good())
	    {
	      if (str_ptr)
		{
		  os << "Opening file ``" << tokens[next_tok] << "''\n";
		  pdboverlay empty;
		  *str_ptr = empty;
		  fp >> * str_ptr;
		  os << "Read in " << str_ptr->n_atoms() << " atoms\n";
		  *al_ptr = *str_ptr;
		  (*al_ptr)[0].title() = tokens[next_tok];
		  if (al_ptr == &align[0]) os << "Replaced ALIGN\n";
		  else if (al_ptr == &align[1])
		    {
		      os << "Replaced SEQ1\n";
		      sm.have_stable = false;
		      sm.have_table = false;
		      sm.have_seq1 = true;
		    }
		  else if (al_ptr == &align[2])
		    {
		      os << "Replaced SEQ2 and OVERLAY\n";
		      strct[0] = strct[2];
		      sm.have_stable = false;
		      sm.have_table = false;
		      sm.have_seq2 = true;
		      sm.have_overlay = true;
		    }
		}
	      else
		{
		  os << "Opening file ``" << tokens[next_tok] << "''\n";
		  al_ptr->clear();
		  fp >> * al_ptr;
		  os << "Read in " << al_ptr->n_seq() << " sequence(s)\n";
		  sm.have_table = false;
		  if (al_ptr == &align[0]) sm.have_align = true;
		  else if (al_ptr == &align[1]) sm.have_seq1 = true;
		  else if (al_ptr == &align[2]) sm.have_seq2 = true;
		}
	    }
	  else 
	    {
	      os << "Filename ``" << tokens[next_tok] << "'' not found\n";
	      return 1;
	    }
	}
      else 
	{
	  os << "You must specify a filename\n";
	  return 1;
	}
      return 0;
    }

  else if (upcase(command.at(0,4)) == "ROAD")
    {
      if (tokens.dim() > next_tok)
	{
	  sm.lowroad=false;
	  sm.highroad=false;
	  sm.midroad=false;
	  sm.randroad=false;
	  if (upcase(tokens[next_tok].at(0,2)) == "HI") 
	    sm.highroad=true;
	  else if (upcase(tokens[next_tok].at(0,3)) == "LOW") 
	    sm.lowroad = true;
	  else if (upcase(tokens[next_tok].at(0,3)) == "MID") 
	    sm.midroad = true;
	  else sm.randroad = true;
	}
      os << "Decision strategy is now set to ";
      if (sm.highroad) os << "HIGHROAD\n";
      else if (sm.lowroad) os << "LOWROAD\n";
      else if (sm.midroad) os << "MIDROAD\n";
      else os << "RANDOM\n";
      return 0;
    }


  else if (upcase(command.at(0,4)) == "RUNL")
    {
      if (tokens.dim() > next_tok)
	{
	  sm.runlength() = atoi(tokens[next_tok].c_str());
	}
      os << "Minimum equivalence RUNLENGTH is now set to ";
      os << sm.runlength() << " residues\n";
      return 0;
    }


  else if (upcase(command.at(0,6)) == "SALIGN")
    {
      if (!(align[1].n_seq() && align[2].n_seq()))
        {
          os << "Need two sequences to align!!!\n";
          os << "(Use READ command)\n";
          return 1;
        }
      if (!sm.have_stable)
	{
	  os << "Running STABULATE command...\n";
          myoshook(os);
	  strcpy(buffer, "STABULATE");
	  if (this->parse_and_run(buffer, os)) return 1;
	}
      os << "Aligning sequences...\n";
      myoshook(os);
      align[0] = sm.dyn_prog();
      os << align[0].newgaps() << " gaps added to make alignment\n";
      os << "Alignment score is " << align[0].score_val << "\n";
      sm.have_align = true;
      return 0;
    }


  else if (upcase(command.at(0,4)) == "SEED")
    {
      if (tokens.dim() > next_tok)
        sm.random_seed() = atoi(tokens[next_tok].c_str());
      os << "Randomization seed is now set to " << sm.random_seed() << "\n";
      return 0;
    }

  else if (upcase(command.at(0,5)) == "SPLIT")
    {
      if (align[0].n_seq() < 2)
        {
          os << "At least two sequences required to split!!\n";
          return 1;
        }
      uint extract_seq = 1;
      if (tokens.dim() > next_tok)
        {
          extract_seq = atoi(tokens[next_tok].c_str());
        }
      if ((extract_seq < 1) || (extract_seq > align[0].n_seq()))
        {
          os << "There is no sequence " << extract_seq << " !!\n";
          return 1;
        }
      --extract_seq; // to make it index the C-like matrix
      uint j;
      Alignment empty;
      os << "Placing sequence " << extract_seq+1 << " in SEQ1...\n";
      myoshook(os);
      align[1] = align[0][extract_seq].remove_gaps();
      os << "Placing others in SEQ2...\n";
      align[2] = empty;
      for (j=0; j < align[0].n_seq(); ++j)
        if (j != extract_seq) align[2] += align[0][j];
      align[2] = align[2].remove_gaps();
      sm.have_seq1 = true;
      sm.have_seq2 = true;
      os << "reTABULATEing scores...\n";
      if (this->parse_and_run("TABULATE", os)) return 1;
      return 0;
    }

  else if (upcase(command.at(0,5)) == "STABU")
    {
      if (!(strct[1].ca().dim() && strct[0].ca().dim()))
	{
	  os << "Need some alpha carbons to compare!\n";
	  return 1;
	}
      strct[1].update_refatoms();
      strct[0].update_refatoms();
      strct[1].update_orientations();
      strct[0].update_orientations();
      sm.stabulate(strct[1], strct[0]);
      os << "Mean residue alignment score = " << sm.mean();
      sm.sigma = sm.std_dev();
      os << " with a standard deviation of " << sm.sigma << "\n";
      sm.have_table = false;
      sm.have_stable = true;
      sm.have_align = false;
      return 0;
    }


  else if (upcase(command.at(0,10)) == "STRUCSTATS")
    {
      if (!(align[0].dim() && strct[0].ca().dim()))
	{
	  os << "Need some residues to compare!\n";
	  return 1;
	}
      strct[0].update_refatoms();
      strct[0].update_orientations();
      strct[1].update_refatoms();
      strct[1].update_orientations();
      align[0].transfer_atoms(strct[1], strct[0]);
      uint i;
      // Just output first structure in each register
      uint off = strct[1].dim();
      // Make sure align has enough sequences
      if ((align[0].dim() < 1) || (align[0].dim() < off))
	{
	  os << "Not enough sequences in ALIGN\n";
	  return 1;
	}
      for (i=0; i < align[0][0].dim(); ++i)
	{
	  const SeqRes & r1 = align[0][0][i];
	  const SeqRes & r2 = align[0][off][i];
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

      return 0;
    }


  else if (upcase(command.at(0,6)) == "SUBOPT")
    {
      if (tokens.dim() > next_tok)
	sm.suboptimal() = atof(tokens[next_tok].c_str());
      os << "The suboptimal alignment factor is now ";
      os << sm.suboptimal() << "\n";
      return 0;
    }


  // superpose struct2 onto struct1; put result in overlay
  else if (upcase(command.at(0,6)) == "SUPERP")
    {
      os << "Superposing STRUCT2 onto STRUCT1\n";
      if (strct[1].n_atoms() < 4)
	{
	  os << "STRUCT1 needs at least 4 atoms!\n";
	  return 1;
	}
      if (strct[2].n_atoms() < 4)
	{
	  os << "STRUCT2 needs at least 4 atoms!\n";
	  return 1;
	}
      if (align[0].dim() < 2)
	{
	  os << "ALIGN must have at least 2 sequences to use for overlay!\n";
	  os << "Running ALIGN, to see if this helps...\n";
	  if (this->parse_and_run("ALIGN", os)) return 1;
	  else if (align[0].dim() < 2) return 1;
	}
      if (align[0][0].n_cap_aligned(align[0][1]) < 4)
	{
	  os << "ALIGN must have at least 4 capitalized equivalences!\n";
	  return 1;
	}
      RMat xform;
      os.setf(ios::showpoint);
      os.precision(4);
      os << "Starting RMS deviation is " 
	<< align[0].rms(strct[1],strct[0])
	  << " for " 
	    << align[0][0].n_cap_aligned(align[0][1]) 
	      << " equivalent matches.\n";
      os << "Starting RMS deviation is " 
	<< align[0].upcase().rms(strct[1],strct[0])
	  << " for "
	    << align[0][0].n_aligned(align[0][1]) 
	      << " total matches.\n";
      os << "Transforming STRUCT2...\n";
      os << "calculating rotation matrix...\n";
      myoshook(os);
      xform = overlay(strct[1],strct[2],align[0]);
      os << "applying matrix...\n";
      strct[0] = strct[2].transform(xform);
      sm.have_stable = false;
      os << "New RMS deviation is " 
	<< align[0].rms(strct[1],strct[0]) 
	  << " for " 
	    << align[0][0].n_cap_aligned(align[0][1]) 
	      << " equivalent matches.\n";
      os << "New RMS deviation is " 
	<< align[0].upcase().rms(strct[1],strct[0]) 
	  << " for " 
	    << align[0][0].n_aligned(align[0][1]) 
	      << " total matches.\n";
      return 0;
    }


  else if (upcase(command.at(0,4)) == "TABU")
    {
      if (!(align[1].n_seq() && align[2].n_seq()))
        {
          os << "Need two sequences/alignments to tabulate!!!\n";
          os << "(Use READ command)\n";
          return 1;
        }
      sm.tabulate(align[1], align[2], cm);
      os << "Mean residue alignment score = " << sm.mean();
      sm.sigma = sm.std_dev();
      os << " with a standard deviation of " << sm.sigma << "\n";
      sm.have_table = true;
      sm.have_stable = false;
      sm.have_align = false;
      return 0;
    }


  // toggle setting to use orientation in structural alignment
  else if (upcase(command.at(0,5)) == "USEAN")
    {
      if (sm.useangle()) 
	{
	  sm.useangle() = false;
	  os << "Orientation will NOT be used in coordinate alignment\n";
	}
      else 
	{
	  sm.useangle() = true;
	  os << "Orientation WILL be used in coordinate alignment\n";
	}
      sm.have_stable = false; // affects tabulation
      return 0;
    }


  else if (upcase(command.at(0,5)) == "WEIGH")
    {
      align[0].set_weights();
      align[1].set_weights();
      align[2].set_weights();
      sm.have_table = false; // affects tabulation
      os << "SEQ1, SEQ2, and ALIGN have been Bruns-weighted\n";
      sm.have_table = false;
      return 0;
    }

  // write [seq1|seq2|align] [filename]
  else if (upcase(command.at(0,5)) == "WRITE")
    {
      str_ptr = NULL;
      al_ptr = NULL;
      ostream * destination = & os;
      if (tokens.dim() <= next_tok)
	al_ptr = & align[0];
      else
	{
	  if (upcase(tokens[next_tok]) == "SEQ1")
	    {al_ptr = & align[1]; ++ next_tok;}
	  else if (upcase(tokens[next_tok]) == "SEQ2")
	    {al_ptr = & align[2]; ++ next_tok;}
	  else if (upcase(tokens[next_tok]) == "ALIGN")
	    {al_ptr = & align[0]; ++ next_tok;}
	  else if (upcase(tokens[next_tok]) == "STRUCT2")
	    {str_ptr = & strct[2]; ++ next_tok;}
	  else if (upcase(tokens[next_tok]) == "STRUCT1")
	    {str_ptr = & strct[1]; ++ next_tok;}
	  else if (upcase(tokens[next_tok]) == "OVERLAY")
	    {str_ptr = & strct[0]; ++ next_tok;}
	  else al_ptr = & align[0];
	}
      if (tokens.dim() <= next_tok)
	destination = & os;
      else
	{
	  ofp.open(tokens[next_tok].c_str());
	  if (ofp.good()) 
	    {
	      os << "Opening file ``" << tokens[next_tok] << "''\n";
              myoshook(os);
	      destination = & ofp;
	    }
	  else 
	    {
	      os << "Failed to open file ``" << buffer << "''\n";
	      return 1;
	    }
	}
      if (al_ptr)
	{
	  *destination << * al_ptr;
	  os << "Wrote " << al_ptr->n_seq() << " sequence(s)\n";
	}
      else if (str_ptr)
	{
	  if (str_ptr == &strct[0]) // Is the overlay, 
	    // needs special handling
	    {
	      pdboverlay special(0);
	      uint str;
	      for (str = 0; str < strct[1].dim(); ++str)
		{
		  special += strct[1][str];
		  special[str].sequence() = align[0][str];
		}
	      for (str = strct[1].dim(); 
		   str < (strct[1].dim() + strct[0].dim()); 
		   ++str)
		{
		  special += strct[0][str - strct[1].dim()];
		  special[str].sequence() = align[0][str];
		}
	      *destination << special;
	    }
	  else
	    {
	      *destination << * str_ptr;
	      os << "Wrote " << str_ptr->n_atoms() << " atom(s)\n";
	    }
	}
      return 0;
    }

  else
    {
      os << "Unknown command: ``" << command << "''\n";
      return 1;
    }
  
  return 1;
}

