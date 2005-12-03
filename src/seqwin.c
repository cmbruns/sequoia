#include <Xm/MainW.h>
#include <Xm/DrawingA.h>
#include <Xm/FileSB.h>
#include <Xm/RowColumn.h>
#include <stdio.h>

/* $Id: seqwin.c,v 1.2 2001/11/29 00:04:43 bruns Exp $ */
/* $Header: /usr/data/cvs/sequoia1/seqwin.c,v 1.2 2001/11/29 00:04:43 bruns Exp $ */
/* $Log: seqwin.c,v $
/* Revision 1.2  2001/11/29 00:04:43  bruns
/* Removed ^M characters
/* Added cvs tags
/* */

extern Widget toplevel;
extern void getseq(int seqnum);

/* these are defined in xoshook.cxx */
extern int seqed_nres;
extern int seqed_nseq;
extern char ** seqed_array;

void seqwin_activate();

static Widget seq_mainw, seq_shell, seq_drawa;
static GC textgc;
static char read_reg[100] = "SEQ1";
static Display * dpy;
/* sequence layout variables */
static int top_margin = 20;
static int bottom_margin = 20;
static int left_margin = 10;
static int right_margin = 10;
static Font font;
static XFontStruct * font_info;

static void expose_seqdraw(Widget w, XtPointer client_data, XtPointer
			   call_data);

void resize_seqdraw(Widget w, XtPointer client_data, XEvent * event,
                 Boolean *continue_to_dispatch); /* */

static void file_cb(Widget widget, XtPointer client_data, XtPointer
		    call_data); 
static void option_cb(Widget widget, XtPointer client_data, XtPointer call_data);
static void loadreg_cb(Widget widget, XtPointer client_data, XtPointer call_data);
static void resize_drawa();
static int alignheight();
static int linelength();

/* Eventually resize events should change the horizontal size of all
   of the window objects, while the vertical size is only changed in
   the main window, with the internal stuff accessed via a vertical
   scrollbar
   */

void seqwin_activate()
{
  /* If the window has not yet been created, create it */
  if (!seq_shell)
    {
      int mheight = 300;
      int dheight;
      int mwidth = 300;
      int dwidth;
      int vbar_w = 0;
      int barspace = 0;
      short vbar_margin_w = 0;
      short margin_h = 0; /* =0 to avoid leftover bits upon assignment */
      short margin_w = 0;
      XmString file, edit, help; /* Menus */
      XmString loadfile, loadreg, savefile, savereg, close; /* File Menu */
      XmString seqonelab, seqtwolab, alignlab; /* File->Load Register Menu */
      Widget menu, menubar, vbar, pullright;

      /* Need a background shell for the window manager */
      dpy = XtDisplay(toplevel);
      seq_shell = XtVaAppCreateShell(NULL, 
		    "XSequoia", topLevelShellWidgetClass, dpy,
		    XtNtitle, "XSequoia Sequence Viewer",
		    /* Don't destroy window on deletion */
		    XmNdeleteResponse, XmUNMAP,
		    NULL);
      
      /* A main window on this will permit a menubar */
      seq_mainw = XtVaCreateManagedWidget
	("seq_mainw",
	 xmMainWindowWidgetClass, seq_shell,
	 /* only show bars when needed */
	 XmNscrollBarDisplayPolicy, XmAS_NEEDED,
	 /* you move the children for me */
	 XmNscrollingPolicy, XmAUTOMATIC,
	 XmNheight, mheight,
	 XmNwidth, mwidth,
	 NULL);
      XtManageChild(seq_mainw);

      /* Menubar */
      file = XmStringCreateLocalized("File");
      edit = XmStringCreateLocalized("Edit");
      help = XmStringCreateLocalized("Help");
      menubar = XmVaCreateSimpleMenuBar(seq_mainw, "menubar",
					XmVaCASCADEBUTTON, file, 'F',
					XmVaCASCADEBUTTON, edit, 'E',
					XmVaCASCADEBUTTON, help, 'H',
					NULL);
      /* Tell which is the help button */
      if ((menu = XtNameToWidget(menubar, "button_2")))
	XtVaSetValues (menubar, XmNmenuHelpWidget, menu, NULL);
      XmStringFree(file);
      XmStringFree(edit);
      XmStringFree(help);
      /* Disable option for now */
      /* XtSetSensitive(XtNameToWidget(menubar, "button_0"), False); /* File */
      XtSetSensitive(XtNameToWidget(menubar, "button_1"), False); /* Edit */
      XtSetSensitive(XtNameToWidget(menubar, "button_2"), False); /* Help */
      XtManageChild(menubar);

      /* Now the actual menu buttons */
      /* File Menu */
      loadfile = XmStringCreateLocalized("Load from File...");
      loadreg = XmStringCreateLocalized("Load from Register");
      savefile = XmStringCreateLocalized("Save to File...");
      savereg = XmStringCreateLocalized("Save to Register");
      close = XmStringCreateLocalized("Close Window");
      menu = XmVaCreateSimplePulldownMenu(menubar, "file_menu", 0, file_cb,
                     XmVaPUSHBUTTON, loadfile, 'L', NULL, NULL,
                     XmVaCASCADEBUTTON, loadreg, 'R',
                     XmVaPUSHBUTTON, savefile, 'S', NULL, NULL,
                     XmVaCASCADEBUTTON, savereg, 'G',
                     XmVaSEPARATOR,
                     XmVaPUSHBUTTON, close, 'C', NULL, NULL,
                     NULL);
      XmStringFree(loadfile);
      XmStringFree(loadreg);
      XmStringFree(savefile);
      XmStringFree(savereg);
      XmStringFree(close);
      /* Disable option for now */
      XtSetSensitive(XtNameToWidget(menu, "button_0"), False); /* Load file*/
      /* XtSetSensitive(XtNameToWidget(menu, "button_1"), False); /* Register */
      XtSetSensitive(XtNameToWidget(menu, "button_2"), False); /* Save */
      XtSetSensitive(XtNameToWidget(menu, "button_3"), False); /* Save */
      /* XtSetSensitive(XtNameToWidget(menu, "button_4"), False); /* Close */

      /* File->Load Register menu */
      seqonelab = XmStringCreateLocalized("SEQ1");
      seqtwolab = XmStringCreateLocalized("SEQ2");
      alignlab = XmStringCreateLocalized("ALIGN");
      pullright = XmVaCreateSimplePulldownMenu (menu,
		      "projection", 1 /* menu item offset */, loadreg_cb,
                      XmVaPUSHBUTTON, seqonelab, 'S', NULL, NULL,
                      XmVaPUSHBUTTON, seqtwolab, 'Q', NULL, NULL,
                      XmVaPUSHBUTTON, alignlab, 'A', NULL, NULL,
                      NULL);
      XmStringFree(seqonelab);
      XmStringFree(seqtwolab);
      XmStringFree(alignlab);

      /* The drawing area permits general graphics operations */
      dheight = mheight - 2*margin_h;
      /* must be narrower than main window to avoid automatic */
      /* horizontal scrollbar */
      /* need margin sizes to determine exact drawing area size */
      XtVaGetValues(seq_mainw,
		    XmNmainWindowMarginHeight, &margin_h,
		    XmNmainWindowMarginWidth, &margin_w,
		    XmNverticalScrollBar, &vbar,
		    XmNspacing, &barspace,
		    NULL);
      if (vbar == NULL) vbar_w = vbar_margin_w = 0;
      else XtVaGetValues(vbar,
		    XmNwidth, &vbar_w,
		    XmNmarginWidth, &vbar_margin_w,
		    NULL);
      /* 8 is enough, 7 is not (Linux 2, Moo-Tif) */
      dwidth = mwidth - 2*margin_w - vbar_w - 2*vbar_margin_w -
	barspace - 8;
      if (dheight < 1) dheight = 1;
      if (dwidth < 1) dwidth = 1;
      seq_drawa = XtVaCreateManagedWidget("seq_drawa",
		    xmDrawingAreaWidgetClass, seq_mainw,
		    XmNheight, dheight,
		    XmNwidth, dwidth,
		    XmNresizePolicy, XmRESIZE_ANY,
		    NULL);
      /* resize does not get sent with scrolled children */
      /* XtAddCallback(seq_drawa, XmNresizeCallback, resize_seqdraw, NULL); */
      XtAddCallback(seq_drawa, XmNexposeCallback, expose_seqdraw, NULL);
      XtManageChild(seq_drawa);

      /* Yet another attempt to simultaneously have automatic */
      /* scrollbars while getting resize events.  Further it appears */
      /* that the SHELL is the only window for which I can capture */
      /* these events (at least with ResizeRedirectMask). */
      XtAddEventHandler(seq_shell, StructureNotifyMask,
                        False, resize_seqdraw, NULL);

      XtVaSetValues(seq_mainw, 
		    XmNmenuBar, menubar,
		    XmNworkWindow, seq_drawa, 
		    NULL);
      XtRealizeWidget(seq_shell);

      textgc = XCreateGC(dpy, RootWindowOfScreen(XtScreen(seq_drawa)), 0,
			 NULL);
      /* Font selection */
      font = XLoadFont(dpy, "fixed");
      XSetFont(dpy, textgc, font);
      font_info = XQueryFont(dpy, font);
      XSetForeground(dpy, textgc, BlackPixelOfScreen(XtScreen(seq_mainw)));
      XSetBackground(dpy, textgc, WhitePixelOfScreen(XtScreen(seq_mainw)));
    }

  /* Pop up sequence window, if not already popped up */
  XtPopup(seq_shell, XtGrabNone);
  /* If it was already up, bring it to the front */
  XMapRaised(XtDisplay(seq_shell), XtWindow(seq_shell));
}


/* how many characters can we put on one line of sequence? */
int linelength()
{
  int cwidth, maxline;
  Dimension dwidth = 0;
  XtVaGetValues(seq_drawa,
		XmNwidth, &dwidth,
		NULL);
  /* width of characters */
  cwidth = font_info->max_bounds.width;
  if (cwidth < 1) cwidth = 1;
  /* how many characters can be put on a line? */
  maxline = (dwidth - left_margin - right_margin)/cwidth;
  if (maxline < 1) maxline = 0;
  return maxline;
}


/* how many vertical pixels does the alignment require? */
int alignheight()
{
  int rows;
  int columns;
  int pixels;
  int cheight;

  columns = linelength();
  if (columns < 1) columns = 1;
  rows = (seqed_nres/columns) * (seqed_nseq + 2);
  if (seqed_nres%columns) rows += (seqed_nseq + 2);
  cheight = font_info->max_bounds.ascent + font_info->max_bounds.descent;
  pixels = cheight * rows + top_margin + bottom_margin;
  return pixels;
}


/* resize drawing area to fit main window and alignment. */
/* Actual XConfigure events are sent to resize_seqdraw(), which then */
/* calls this function */
void resize_drawa()
{
  Dimension width = 0;
  Dimension height = 0;
  /* menubar */
  Dimension mheight = 50;
  /* drawing area */
  Dimension dwidth, dheight, pheight;
  Dimension barspace = 0;
  Dimension margin_h = 0;
  Dimension margin_w = 0;
  Widget menubar, vbar;
  Dimension vbar_w = 0;
  short vbar_margin_w = 0;

  XtVaGetValues(seq_mainw,
		XmNheight, &height,
		XmNwidth, &width,
		XmNmarginHeight, &margin_h,
		XmNmarginWidth, &margin_w,
		XmNmenuBar, &menubar,
		XmNverticalScrollBar, &vbar,
		XmNspacing, &barspace,
		NULL);
  XtVaGetValues(menubar,
		XmNheight, &mheight,
		NULL);
  /* Don't crash with really small windows */
  if ((mheight + margin_h) >= height) return;
  if (margin_w >= width) return;

  if (vbar == NULL) vbar_w = vbar_margin_w = 0;
  else XtVaGetValues(vbar,
		     XmNwidth, &vbar_w,
		     XmNmarginWidth, &vbar_margin_w,
		     NULL);
  /* 8 is enough, 7 is not (Linux 2, Moo-Tif) */
  dwidth = width - 2*margin_w - vbar_w - 2*vbar_margin_w -
    barspace - 8;
  if (dwidth < 1) return;
  /* Height should be the maximum of the space provided by the */
  /* main window, and the space required by the alignment */
  dheight = height - mheight;
  pheight = alignheight();
  if (pheight > dheight) dheight = pheight;
  if (dheight < 1) return;
  XtVaSetValues(seq_drawa,
		XmNwidth, dwidth,
		XmNheight, dheight,
		NULL); /* */
}


void expose_seqdraw(Widget w, XtPointer client_data, XtPointer call_data)
{
  int cheight, maxline, thisline, i;
  int voffset;
  int j;

  cheight = font_info->max_bounds.ascent + font_info->max_bounds.descent;
  /* how many characters can be put on a line? */
  maxline = linelength();
  /* none?, then do nothing */
  if (maxline < 1) return;
  voffset = top_margin;
  /* First remove old junk from display */
  XClearWindow(dpy, XtWindow(w));
  /* Now output lines of sequence */
  if (seqed_nseq > 0)
    {
      for (i=0; i < seqed_nres; i += maxline)
	{
	  /* avoid junk on last line of sequence */
	  if ((i + maxline) > seqed_nres) thisline = (seqed_nres - i);
	  else thisline = maxline;
	  for (j=0; j < seqed_nseq; ++j)
	    {
	      XDrawImageString(dpy, XtWindow(w), textgc, left_margin,
			  voffset, seqed_array[j] + i, thisline); 
	      voffset += cheight;
	    }
	  /* put some space between each group */
	  voffset += 2*cheight;
	}
    }
  else
    XDrawString(dpy, XtWindow(w), textgc, left_margin, top_margin, 
		"No sequence ", 11);  
}


/* Call this for actual possibly resizing X-events; call resize_drawa */
/* to directly do the resizing */
/* when you change the size of the window, you must change the size of */
/* the drawing area, hence the following function */
void resize_seqdraw(Widget w, XtPointer client_data, XEvent * event,
                 Boolean *continue_to_dispatch)
{
  /* we really only want resize events */
  if (event->type == ConfigureNotify)
    {  
      static Dimension oldheight;
      static Dimension oldwidth;
      Dimension width = event->xconfigure.width;
      Dimension height = event->xconfigure.height;
      if ((height != oldheight) && (width != oldwidth))
	{
	  /* now we know there was a resize event */
	  /* for some reason the scrollbar does not always update the */
	  /* first time, so do it twice (Linux Mootif) */
	  resize_drawa();
	  resize_drawa();
	  oldheight = height;
	  oldwidth = width;
	}
    }
}


/* File menu callback function */
void file_cb(Widget widget, XtPointer client_data, XtPointer call_data)
{
  int item_no = (int) client_data;
  static Widget fdialog;

  switch (item_no)
    {
    case 0: /* Load from file */
      if (!fdialog)
        {
	  XmString seqone, seqtwo, align;
	  XmString reg, okbutton, cancelbutton;
	  Widget rc, option_menu;

          fdialog = XmCreateFileSelectionDialog (seq_shell, "file_sel",
						 NULL, 0);
          /* Change the default button names */
 	  okbutton = XmStringCreateLocalized("Read File");
 	  cancelbutton = XmStringCreateLocalized("Dismiss");
	  XtVaSetValues(fdialog,
			XmNokLabelString, okbutton,
			XmNcancelLabelString, cancelbutton,
			NULL); /* */
	  /* make option menu for what to do with the file */
	  rc = XtVaCreateManagedWidget("rowcol",
				       xmRowColumnWidgetClass, fdialog, 
				       NULL);
	  seqone = XmStringCreateLocalized("SEQ1");
	  seqtwo = XmStringCreateLocalized("SEQ2");
	  align = XmStringCreateLocalized("ALIGN");
	  reg = XmStringCreateLocalized("Register:");
	  option_menu = XmVaCreateSimpleOptionMenu (rc, "option_menu",
			  reg, 'R', 0 /* initial selection */, option_cb,
                          XmVaPUSHBUTTON, seqone, 'S', NULL, NULL,
                          XmVaPUSHBUTTON, seqtwo, 'Q', NULL, NULL,
                          XmVaPUSHBUTTON, align, 'A', NULL, NULL,
                          NULL);
	  XtManageChild(option_menu);
          /* XtAddCallback(fdialog, XmNokCallback, load_seq, NULL); */
          XtAddCallback(fdialog, XmNcancelCallback,
                        (void *) XtUnmanageChild, NULL);
          /* No help exists, so disable */
	  XtUnmanageChild(XmFileSelectionBoxGetChild(fdialog,
						     XmDIALOG_HELP_BUTTON)); 
	  XmStringFree(seqone);
	  XmStringFree(seqtwo);
	  XmStringFree(align);
        }
      /* make the dialog close automatically when done */
      XtVaSetValues (fdialog,
                     XmNautoUnmanage, True,
                     NULL); /* */
      XtManageChild(fdialog);
      XtPopup(XtParent(fdialog), XtGrabNone);
      break;
    case 1: /* Load from register */
      /* Do nothing, the pullright pushbuttons have their own */
      /* callback, loadreg_cb */
      break;
    case 2: /* save to file */
      break;
    case 3: /* save to register */
      break;
    case 4: /* close */
      XtPopdown(seq_shell);
      break;
    default:
      break;
    }
}


/* File input chooser option menu callback function */
void option_cb(Widget widget, XtPointer client_data, XtPointer call_data)
{
  int item_no = (int) client_data;
  if (item_no == 0) strcpy(read_reg, "SEQ1");
  else if (item_no == 1) strcpy(read_reg, "SEQ2");
  else if (item_no == 2) strcpy(read_reg, "ALIGN");
  else strcpy(read_reg, "SEQ1");
  return;
}

void loadreg_cb(Widget widget, XtPointer client_data, XtPointer call_data)
{
  int item_no = (int) client_data;

  switch (item_no)
    {
    case 0: /* SEQ1 */
      getseq(1);
      break;
    case 1: /* SEQ2 */
      getseq(2);
      break;
    case 2: /* OVERLAY */
      getseq(0);
      break;
    default:
      break;
    }
  /* redraw window */
  resize_drawa();
  expose_seqdraw(seq_drawa, client_data, call_data);
}
