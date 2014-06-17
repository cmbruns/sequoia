#include "oshook.hxx"

// $Id: xoshook.cxx,v 1.2 2001/11/29 00:04:43 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/xoshook.cxx,v 1.2 2001/11/29 00:04:43 bruns Exp $
// $Log: xoshook.cxx,v $
// Revision 1.2  2001/11/29 00:04:43  bruns
// Removed ^M characters
// Added cvs tags
//

// extern void seqcpy(char * s, const Sequence & seq);

// This file attempts to bridge the C Motif functions with the
// C++ SEQUOIA functions

// must use statically allocated buffer to rewind and reuse (I think)
char * osbuf = new char[OSBUFSIZE];
ostrstream status_os(osbuf, OSBUFSIZE);

// minimal declarations to permit access to Motif Text output window
extern "C" {
#include <Xm/Text.h>
#include <Xm/Command.h>
#include <X11/cursorfont.h>
#include "strctyp.h"
extern Widget text_w, main_menubar, toplevel, command_w;
extern XmTextPosition textpos;
void exec_cmd(Widget widget, XtPointer client_data, XtPointer call_data);
void parseit(char * cmd);
/* sequence editor data */
void getseq(int seqnum);
int seqed_nres = 0;
int seqed_nseq = 0;
char ** seqed_array = NULL;
/* coordinates window data */
int getstrct(int strctnum);
unsigned int strctwin_nstrct = 0;
unsigned int strctwin_refstrct = 0;
coordarray *strctwin_array = NULL;
reslabelp *label_array = NULL;
int * strctwin_size;
}

#include "parse.hxx"
extern cli cmd_line;

// must use as "myoshook(os)", for "os << myoshook" dumps core
ostrstream & myoshook(ostrstream & os)
{
  // Add a null character to ensure the character array is terminated
  os << '\0';

  // output character array to X Window
  char * cp = os.str();
  XmTextReplace(text_w, textpos, textpos, cp);
  textpos += strlen(cp);

  // rewind character array to beginning, to accept further input
  os.seekp(0, ios::beg);
  return os;
}


/* Command callback function (from command window) */
void exec_cmd(Widget widget, XtPointer client_data, XtPointer call_data)
{
  char * cmd;
  XmCommandCallbackStruct *cbs = (XmCommandCallbackStruct *) call_data;

  XmStringGetLtoR(cbs->value, XmFONTLIST_DEFAULT_TAG, &cmd);

  if (!cmd || !*cmd) /* nothing typed? */
    {
      if (cmd)
	XtFree(cmd);
      return;
    }
  parseit(cmd);
  XtFree(cmd);
}

/* for direct execution of commands (from menu selections, say) */
void parseit(char * cmd)
{
  status_os << "\nXSEQUOIA> " << cmd << "\n";
  myoshook(status_os);

  // make clock cursor
  Display * dpy = XtDisplay(toplevel);
  Cursor cursor;
  cursor = XCreateFontCursor(dpy, XC_watch);
  XDefineCursor(dpy, XtWindow(toplevel), cursor);
  XFreeCursor(dpy, cursor);
  // Don't let the user select stuff at this time
  XtSetSensitive(main_menubar, False);
  XtSetSensitive(command_w, False);

  cmd_line.parse_and_run(cmd, status_os);

  // restore cursor
  XUndefineCursor(dpy, XtWindow(toplevel));
  XtSetSensitive(main_menubar, True);
  XtSetSensitive(command_w, True);
  // Give focus to the command input box
  XmProcessTraversal(XmCommandGetChild(command_w, XmDIALOG_COMMAND_TEXT), 
		     XmTRAVERSE_CURRENT);
  
  myoshook(status_os);
}

/* fetch sequence alignment in a form readable by C code */
void getseq(int seqnum)
{
  int i;
  /* delete old structure */
  for (i=0; i < (seqed_nseq); ++i)
    {
      delete [] seqed_array[i];
      seqed_array[i] = NULL;
    }
  delete [] seqed_array;
  seqed_nres = 0;
  seqed_nseq = 0;

  seqed_nres = cmd_line.align[seqnum].length();
  seqed_nseq = cmd_line.align[seqnum].dim();
  seqed_array = new char*[seqed_nseq];
  for (i=0; i < (seqed_nseq); ++i)
    {
      seqed_array[i] = new char[seqed_nres];
      seqcpy(seqed_array[i], cmd_line.align[seqnum][i]);
    }
}


/* Same thing, but reads from a file and does not access a register */
/* fetch sequence alignment in a form readable by C code */
int getseqfile(char * fname)
{
  ifstream infile(fname);
  if (infile.good())
    {
      int i;
      /* delete old structure */
      for (i=0; i < (seqed_nseq); ++i)
	{
	  delete [] seqed_array[i];
	  seqed_array[i] = NULL;
	}
      delete [] seqed_array;
      seqed_nres = 0;
      seqed_nseq = 0;

      Alignment align;
      infile >> align;
      seqed_nres = align.length();
      seqed_nseq = align.dim();
      seqed_array = new char*[seqed_nseq];
      for (i=0; i < (seqed_nseq); ++i)
	{
	  seqed_array[i] = new char[seqed_nres];
	  seqcpy(seqed_array[i], align[i]);
	}
      return 0;
    }
  else return 1;
}


/* fetch PDB structure in a form readable by C code */
int getstrct(int strctnum)
{
  int answer = 1;
  uint i;
  /* delete old structure */
  for (i=0; i < strctwin_nstrct; ++i)
    {
      delete [] strctwin_array[i];
      delete [] label_array[i];
      strctwin_array[i] = NULL;
      label_array[i] = NULL;
    }
  delete [] strctwin_array;
  delete [] label_array;
  delete [] strctwin_size;
  strctwin_nstrct = 0;
  int refnum = 1;

  // overlay should show both struct1 and overlay
  if (strctnum == 0) // OVERLAY
    {
      strctwin_nstrct = cmd_line.strct[refnum].dim() + cmd_line.strct[0].dim();
      strctwin_refstrct = cmd_line.strct[refnum].dim();
    }  
  else 
    {
      strctwin_nstrct = cmd_line.strct[strctnum].dim();
      strctwin_refstrct = 0;
    }
  strctwin_array = new coordarray[strctwin_nstrct];
  label_array = new reslabelp[strctwin_nstrct];
  strctwin_size = new int[strctwin_nstrct];
  int j;
  int nres = 0;
  pdbprotein prot;
  for (i=0; i < (strctwin_nstrct); ++i)
    {
      if (i < strctwin_refstrct)
	prot = cmd_line.strct[refnum][i].ca();
      else
	prot = cmd_line.strct[strctnum][i-strctwin_refstrct].ca();	
      nres = prot.n_atoms();
      if (nres > 0) answer = 0;
      strctwin_array[i] = new vector[nres];
      label_array[i] = new reslabel[nres];
      strctwin_size[i] = nres;
      for (j=0; j < nres; ++ j)
	{
	  strctwin_array[i][j].x = prot[j].coord().x();
	  strctwin_array[i][j].y = prot[j].coord().y();
	  strctwin_array[i][j].z = prot[j].coord().z();
	  // Text label for coordinate
	  sprintf(label_array[i][j], "%s %d %s\0", 
		  prot[j].res_name(), prot[j].res_num(), (char *)prot.title());
	}
    }
  return answer;
}


