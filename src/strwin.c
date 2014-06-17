/* This file is for the XSequoia Coordinate Viewer Window */
#include <Xm/MainW.h>
#include <Xm/DrawingA.h>
#include <Xm/FileSB.h>
#include <Xm/RowColumn.h>
#include <Xm/Label.h>
#include <Xm/ScrolledW.h>
#include <Xm/PanedW.h>
#include <stdio.h>
#include "cvector.h"
#include "strctyp.h"

/* $Id: strwin.c,v 1.2 2001/11/29 00:04:43 bruns Exp $ */
/* $Header: /usr/data/cvs/sequoia1/strwin.c,v 1.2 2001/11/29 00:04:43 bruns Exp $ */
/* $Log: strwin.c,v $
/* Revision 1.2  2001/11/29 00:04:43  bruns
/* Removed ^M characters
/* Added cvs tags
/* */

#define STARTWIDTH 270
#define STARTHEIGHT 200

#define debug 0

extern Widget toplevel;
/* these are defined in xoshook.c */
extern int getstrct(int strctnum);
extern unsigned int strctwin_nstrct;
extern unsigned int strctwin_refstrct;
extern coordarray * strctwin_array;
extern int * strctwin_size;
extern reslabelp *label_array;

void loadreg_cb(Widget widget, XtPointer client_data, XtPointer call_data);

static Widget str_mainw, str_shell, str_drawa, menubar, str_messwin;
static Widget strseq_scrollw, strseq_drawa;
static Display * dpy;
static GC molgc, textgc;
static char read_reg[100] = "STRUCT1";
/* coordinate orientation variables */
static matrix rotmat;
static vector transvec;
static double scale;
static int is_first_structure = 1;
static int cenx = 0;
static int ceny = 0; /* center of drawing area */
static Pixmap scratch;
static int use_slab = 0;
static float slab_size = 10;
static int mouse_moved = 0;

static void file_cb(Widget widget, XtPointer client_data, XtPointer call_data);
static void view_cb(Widget widget, XtPointer client_data, XtPointer call_data);
static void center_cb(Widget widget, XtPointer client_data, XtPointer call_data);
static void option_cb(Widget widget, XtPointer client_data, XtPointer call_data);
static int center_mol();
static int scale_mol();
static void draw_mol();
static void expose_strdraw(Widget w, XtPointer client_data, XtPointer
			   call_data);
static void expose_strseqdraw(Widget w, XtPointer client_data, XtPointer
			   call_data);
static void buttonlayer(Widget w, XtPointer client_data, XEvent * event, 
                 Boolean *continue_to_dispatch);
static void xyrotate(Widget w, XtPointer client_data, XEvent * event, 
             Boolean *continue_to_dispatch);
static void xytranslate(Widget w, XtPointer client_data, XEvent * event, 
             Boolean *continue_to_dispatch);
static void zoom(Widget w, XtPointer client_data, XEvent * event, 
             Boolean *continue_to_dispatch);
static void zrotate(Widget w, XtPointer client_data, XEvent * event, 
             Boolean *continue_to_dispatch);
static void resize_strdraw(Widget w, XtPointer client_data, XtPointer call_data);
static void update_pixmap();
static void make_selection_string(int x, int y, char * message);


void strwin_activate()
{
  dpy = XtDisplay(toplevel);
  if (!str_shell)
    {
      XmString file, view, help; /* Menus */
      XmString loadfile, loadreg, savefile, savereg, close; /* File Menu */
      XmString stronelab, strtwolab, alignlab; /* File->Load Register
						  Menu */
      XmString center; /* View menu */
      XmString onmolecule, onselection; /* View->Center menu */
      XmString msg; /* message window */
      Widget menu, pullright, panedw;

      rotmat = EYE;

      str_shell = XtVaAppCreateShell
	(NULL, 
	 "XSequoia", topLevelShellWidgetClass, dpy,
	 XtNtitle, "XSequoia Molecule Viewer",
	 /* Don't destroy window on deletion */
	 XmNdeleteResponse, XmUNMAP,
	 NULL);

      str_mainw = XtVaCreateManagedWidget
	("matrix_window",
	 xmMainWindowWidgetClass, str_shell,
	 XmNvisualPolicy, XmVARIABLE,
	 NULL);
      XtManageChild(str_mainw);

      /* Menubar */
      file = XmStringCreateLocalized("File");
      view = XmStringCreateLocalized("View");
      help = XmStringCreateLocalized("Help");
      menubar = XmVaCreateSimpleMenuBar(str_mainw, "menubar",
					XmVaCASCADEBUTTON, file, 'F',
					XmVaCASCADEBUTTON, view, 'V',
					XmVaCASCADEBUTTON, help, 'H',
					NULL);
      /* Tell which is the help button */
      if ((menu = XtNameToWidget(menubar, "button_2")))
	XtVaSetValues (menubar, XmNmenuHelpWidget, menu, NULL);
      XmStringFree(file);
      XmStringFree(view);
      XmStringFree(help);
      /* Disable option for now */
      /* XtSetSensitive(XtNameToWidget(menubar, "button_0"), False); /* File */
      /* XtSetSensitive(XtNameToWidget(menubar, "button_1"), False); /* View */
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
      stronelab = XmStringCreateLocalized("STRUCT1");
      strtwolab = XmStringCreateLocalized("STRUCT2");
      alignlab = XmStringCreateLocalized("OVERLAY");
      pullright = XmVaCreateSimplePulldownMenu (menu,
		      "projection", 1 /* menu item offset */, loadreg_cb,
                      XmVaPUSHBUTTON, stronelab, 'S', NULL, NULL,
                      XmVaPUSHBUTTON, strtwolab, 'T', NULL, NULL,
                      XmVaPUSHBUTTON, alignlab, 'O', NULL, NULL,
                      NULL);
      XmStringFree(stronelab);
      XmStringFree(strtwolab);
      XmStringFree(alignlab);

      /* View Menu */
      center = XmStringCreateLocalized("Center");
      menu = XmVaCreateSimplePulldownMenu(menubar, "view_menu", 1, view_cb,
                     XmVaCASCADEBUTTON, center, 'C',
                     NULL);
      XmStringFree(center);
      /* Disable option for now */
      /* XtSetSensitive(XtNameToWidget(menu, "button_0"), False); /* Center */

      /* View->Center menu */
      onmolecule = XmStringCreateLocalized("on Molecule");
      onselection = XmStringCreateLocalized("on Selection");
      pullright = XmVaCreateSimplePulldownMenu (menu,
		      "center", 0 /* menu item offset */, center_cb,
                      XmVaPUSHBUTTON, onmolecule, 'M', NULL, NULL,
                      XmVaPUSHBUTTON, onselection, 'S', NULL, NULL,
                      NULL);
      XmStringFree(onmolecule);
      XmStringFree(onselection);
      XtSetSensitive(XtNameToWidget(pullright, "button_1"), False); /* on Selection */

      /* Paned window, parent of structure and sequence viewers */
      panedw = XtVaCreateManagedWidget
	("panedw",
	 xmPanedWindowWidgetClass, str_mainw,
	 NULL);
      XtManageChild(panedw);

      /* *** Drawing Area *** */
      str_drawa = XtVaCreateManagedWidget
	("str_drawa",
	 xmDrawingAreaWidgetClass, panedw,
	 XmNwidth, STARTWIDTH,
	 XmNheight, STARTHEIGHT,
	 NULL);
      cenx = STARTWIDTH/2;
      ceny = STARTHEIGHT/2;
      /* Capture button events */
      XtAddEventHandler(str_drawa, ButtonPressMask | ButtonReleaseMask, 
                        False, buttonlayer, NULL);
      XtAddEventHandler(str_drawa, Button1MotionMask, 
                        False, xyrotate, NULL);
      XtAddEventHandler(str_drawa, Button2MotionMask, 
                        False, xytranslate, NULL);
      XtAddEventHandler(str_drawa, Button3MotionMask, 
                        False, zoom, NULL);
      XtAddCallback(str_drawa, XmNexposeCallback, expose_strdraw, NULL);
      XtAddCallback(str_drawa, XmNresizeCallback, resize_strdraw, NULL);
      XtManageChild(str_drawa);

      /* Graphics context for molecular drawing */
      molgc = XCreateGC(dpy, RootWindowOfScreen(XtScreen(str_mainw)), 0,
			 NULL);
      XSetForeground(dpy, molgc, BlackPixelOfScreen(XtScreen(str_mainw)));
      XSetBackground(dpy, molgc, WhitePixelOfScreen(XtScreen(str_mainw)));

      /* Pixmap for double buffering of drawing area */
      scratch = XCreatePixmap(dpy,
			      RootWindowOfScreen(XtScreen(str_mainw)),
			      STARTWIDTH, STARTHEIGHT,
			      DefaultDepth(dpy, DefaultScreen(dpy)));
      if (scratch == XmUNSPECIFIED_PIXMAP)
	{
	  printf("can't create initial pixmap\n");
	  exit(1);
	}

      /* Sequence viewing window */
      strseq_scrollw = XtVaCreateManagedWidget
	("strseq_scrollw",
	 xmScrolledWindowWidgetClass, panedw,
	 XmNwidth, STARTWIDTH,
	 XmNheight, 50,
	 XmNscrollBarDisplayPolicy, XmAS_NEEDED,
	 XmNvisualPolicy, XmCONSTANT,
	 XmNscrollingPolicy, XmAUTOMATIC,
	 XmNskipAdjust, True, /* Don't automatically resize from paned window */
	 NULL);

      /* Structure Sequence Drawing Area */
      strseq_drawa = XtVaCreateManagedWidget
	("strseq_drawa",
	 xmDrawingAreaWidgetClass, strseq_scrollw,
	 XmNwidth, STARTWIDTH + 20,
	 XmNheight, 40,
	 NULL);
      XtAddCallback(strseq_drawa, XmNexposeCallback, expose_strseqdraw, NULL);
      textgc = XCreateGC(dpy, RootWindowOfScreen(XtScreen(strseq_drawa)), 0,
			 NULL);
      XtManageChild(strseq_drawa);

      /* message window */
      msg = XmStringCreateLocalized("Your message here");
      str_messwin = XtVaCreateManagedWidget
	("str_messwin",
	 xmLabelWidgetClass, str_mainw,
	 XmNlabelString, msg,
	 NULL);
      XmStringFree(msg);

      XtVaSetValues(str_mainw,
		    XmNmenuBar, menubar,
		    XmNworkWindow, panedw,
		    XmNmessageWindow, str_messwin,
		    NULL);
    }
  XtPopup(str_shell, XtGrabNone);
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
	  XmString strone, strtwo, align;
	  XmString reg, okbutton, cancelbutton;
	  Widget rc, option_menu;

          fdialog = XmCreateFileSelectionDialog (str_shell, "file_sel",
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
	  strone = XmStringCreateLocalized("STRUCT1");
	  strtwo = XmStringCreateLocalized("STRUCT2");
	  align = XmStringCreateLocalized("OVERLAY");
	  reg = XmStringCreateLocalized("Register:");
	  option_menu = XmVaCreateSimpleOptionMenu (rc, "option_menu",
			  reg, 'R', 0 /* initial selection */, option_cb,
                          XmVaPUSHBUTTON, strone, 'S', NULL, NULL,
                          XmVaPUSHBUTTON, strtwo, 'T', NULL, NULL,
                          XmVaPUSHBUTTON, align, 'O', NULL, NULL,
                          NULL);
	  XtManageChild(option_menu);
          /* XtAddCallback(fdialog, XmNokCallback, load_str, NULL); */
          XtAddCallback(fdialog, XmNcancelCallback,
                        (void *) XtUnmanageChild, NULL);
          /* No help exists, so disable */
	  XtUnmanageChild(XmFileSelectionBoxGetChild(fdialog,
						     XmDIALOG_HELP_BUTTON)); 
	  XmStringFree(strone);
	  XmStringFree(strtwo);
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
      XtPopdown(str_shell);
      break;
    default:
      break;
    }
}


/* View menu callback function */
void view_cb(Widget widget, XtPointer client_data, XtPointer call_data)
{
}


/* View->Center */
void center_cb(Widget widget, XtPointer client_data, XtPointer call_data)
{
  int item_no = (int) client_data;

  switch (item_no)
    {
    case 0: /* on Molecule */
      center_mol();
      draw_mol();
      break;
    case 1: /* on Selection */
      break;
    default:
      break;
    }
}

/* File input chooser option menu callback function */
void option_cb(Widget widget, XtPointer client_data, XtPointer call_data)
{
  int item_no = (int) client_data;
  if (item_no == 0) strcpy(read_reg, "STRUCT1");
  else if (item_no == 1) strcpy(read_reg, "STRUCT2");
  else if (item_no == 2) strcpy(read_reg, "OVERLAY");
  else strcpy(read_reg, "STRUCT1");
  return;
}


void loadreg_cb(Widget widget, XtPointer client_data, XtPointer call_data)
{
  int item_no = (int) client_data;
  int result = 1;

  switch (item_no)
    {
    case 0: /* STRUCT1 */
      result = getstrct(1);
      break;
    case 1: /* STRUCT2 */
      result = getstrct(2);
      break;
    case 2: /* OVERLAY */
      result = getstrct(0);
      break;
    default:
      break;
    }
  /* printf("nstructs, size[0]: %d, %d\n",strctwin_nstrct, strctwin_size[0]); */
  if (!result)
    {
      if (is_first_structure)
	{
	  /* center on molecule */
	  result = center_mol();
	  /* initialize scale */
	  result = scale_mol();
	  is_first_structure = 0;
	}
      /* redraw window */
      draw_mol();
    }
}

/* adjust translation vector to center on molecule */
/* ( T = -(c.o.m.)) */
int center_mol()
{
  int answer = 1;
  int i, j;
  int atomcount = 0;
  double xcum = 0;
  double ycum = 0;
  double zcum = 0;

  for (i = 0; i < strctwin_nstrct; ++i)
    {
      for (j = 0; j < strctwin_size[i]; ++j)
	{
	  xcum -= strctwin_array[i][j].x;
	  ycum -= strctwin_array[i][j].y;
	  zcum -= strctwin_array[i][j].z;
	  ++atomcount;
	}
    }
  if (atomcount)
    {
      transvec.x = xcum/atomcount;
      transvec.y = ycum/atomcount;
      transvec.z = zcum/atomcount;
      answer = 0;
    }
  return answer;
}

/* scale model to window */
int scale_mol()
{
  int answer = 1;
  float winsize = 0;
  float protsize = 0;
  Dimension wwidth = 0;
  Dimension wheight = 0;
  float xmax, xmin, ymax, ymin, zmax, zmin;
  float xdiff, ydiff, zdiff;
  int i, j;

  /* 1) figure out scale of window, in pixels */
  XtVaGetValues(str_drawa,
		XmNwidth, &wwidth,
		XmNheight, &wheight,
		NULL);
  /* is the window too small? */
  if (wwidth < 2)
    {
      if (wheight < 2)
	winsize = 20; /* arbitrary small size */
      else
	winsize = wheight;
    }
  else if (wheight < 2) winsize = wwidth;
  else if (wheight > wwidth) winsize = wwidth;
  else winsize = wheight; /* adjust to SMALLER dimension */

  /* 2) figure out scale of coordinates, in Angstroms */
  /* initialize extents to first coordinate */
  if (strctwin_nstrct)
    {    
    if (strctwin_size[0])
      {
	xmin = xmax = strctwin_array[0][0].x;
	ymin = ymax = strctwin_array[0][0].y;
	zmin = zmax = strctwin_array[0][0].z;
      }
    }
  else return 1; /* no structures */
  /* cycle through the coordinates */
  for (i = 0; i < strctwin_nstrct; ++i)
    {
      for (j = 0; j < strctwin_size[i]; ++j)
	{
	  if (xmin > strctwin_array[i][j].x) xmin = strctwin_array[i][j].x;
	  if (ymin > strctwin_array[i][j].y) ymin = strctwin_array[i][j].y;
	  if (zmin > strctwin_array[i][j].z) zmin = strctwin_array[i][j].z;
	  if (xmax < strctwin_array[i][j].x) xmax = strctwin_array[i][j].x;
	  if (ymax < strctwin_array[i][j].y) ymax = strctwin_array[i][j].y;
	  if (zmax < strctwin_array[i][j].z) zmax = strctwin_array[i][j].z;
	}
    }
  xdiff = xmax - xmin;
  ydiff = ymax - ymin;
  zdiff = zmax - zmin;
  /* choose the largest extent */
  if (xdiff > ydiff)
    {
      if (xdiff > zdiff) protsize = xdiff;
      else protsize = zdiff;
    }
  else if (ydiff > zdiff) protsize = ydiff;
  else protsize = zdiff;

  /* fit largest protein dimension into smallest window dimension,
     using a scale in pixels per Angstrom */
  scale = winsize/protsize;
  answer = 0;

  return answer;
}


void draw_mol()
{
  Dimension width = 0;
  Dimension height = 0;

  update_pixmap();
  XtVaGetValues(str_drawa,
                XmNheight, &height,
                XmNwidth, &width,
                NULL);
  XCopyArea(dpy, scratch, XtWindow(str_drawa), molgc, 0, 0, width,
	    height, 0, 0);
}


void update_pixmap()
{
  int i, j;
  vector transform;
  matrix smat;
  int lastplot = 0;
  float oldx, oldy;
  Drawable dwin;
  Dimension width = 0;
  Dimension height = 0;

  Colormap cmap = DefaultColormapOfScreen(XtScreen(str_drawa));
  XColor black, white, red, unused;

  XAllocNamedColor(dpy, cmap, "Black", &black, &unused);
  XAllocNamedColor(dpy, cmap, "White", &white, &unused);
  XAllocNamedColor(dpy, cmap, "Red", &red, &unused);

  XtVaGetValues(str_drawa,
                XmNheight, &height,
                XmNwidth, &width,
                NULL);
  dwin = scratch;
  /* scale matrix */
  smat = matscale(scale, rotmat);
  /* erase old stuff */
  XSetForeground(dpy, molgc, white.pixel);
  XFillRectangle(dpy, dwin, molgc, 0, 0, width, height);
  XSetForeground(dpy, molgc, black.pixel);
  /* cycle through the coordinates */
  for (i = 0; i < strctwin_nstrct; ++i)
    {
      if (i < strctwin_refstrct) XSetForeground(dpy, molgc, red.pixel);
      else XSetForeground(dpy, molgc, black.pixel);
      for (j = 0; j < strctwin_size[i]; ++j)
	{
	  /* translate, rotate, scale coordinate */
	  transform = vectoradd(strctwin_array[i][j], transvec); /* translate */
	  transform = matvmult(smat, transform); /* rotate, scale */
	  transform.x += cenx; /* center in window */
	  transform.y += ceny;
	  if (lastplot) XDrawLine(dpy, dwin, molgc, oldx, oldy,
				  transform.x, transform.y);
	  else XDrawLine(dpy, dwin, molgc, transform.x, transform.y,
			 transform.x, transform.y);
	  oldx = transform.x;
	  oldy = transform.y;
	  lastplot = 1;
	}
      /* don't connect consecutive molecules */
      lastplot = 0;
    }
}


void expose_strdraw(Widget w, XtPointer client_data, XtPointer call_data)
{
  draw_mol();
}


/* Draw sequence in sequence window */
void expose_strseqdraw(Widget w, XtPointer client_data, XtPointer call_data)
{
  XDrawString(dpy, XtWindow(w), textgc, 5, 15, 
	      "Title: (001) ---YOUR--SEQUENCE----HERE--* (030)", 47);
}


void buttonlayer(Widget w, XtPointer client_data, XEvent * event, 
                 Boolean *continue_to_dispatch)
{
      switch (event->xbutton.button)
        {
        case Button1:
          if (debug) printf("Button 1\n");
          xyrotate(w, client_data, event, continue_to_dispatch);
          break;
        case Button2:
          if (debug) printf("Button 2\n");
          xytranslate(w, client_data, event, continue_to_dispatch);
          break;
        case Button3:
          if (debug) printf("Button 3\n");
          zoom(w, client_data, event, continue_to_dispatch);
          break;
        }
}


/* Button 1 Events in the pixmap are for rotation */
void xyrotate(Widget w, XtPointer client_data, XEvent * event, 
             Boolean *continue_to_dispatch) 
{
  static int lastx, lasty;
  static int mouse_moved = 1;
  int x, y;

  switch (event->type)
    {
    case ButtonPress:
      lastx = event->xbutton.x;
      lasty = event->xbutton.y;
      mouse_moved = 0; /* distinguish clicks from drags */
      break;
    case ButtonRelease:
      if (!mouse_moved) /* we have a selection click */
	{
	  static char message[500];
	  XmString msg;

	  make_selection_string(event->xbutton.x, event->xbutton.y, message);
	  // sprintf(message, "%d %d\0", event->xbutton.x, event->xbutton.y);
	  msg = XmStringCreateLocalized(message);
	  XtVaSetValues(str_messwin,
			XmNlabelString, msg,
			NULL);
	  XmStringFree(msg);
	}
      draw_mol();
      mouse_moved = 1;
      break;
    case MotionNotify:
      {
	mouse_moved = 1;
	/* Don't let these events pile up */
	if (XEventsQueued(dpy, QueuedAfterReading) > 0) break;
	else
	  {
	    vector axis;
	    float radius; /* radius of virtual trackball, in pixels */
	    double len;

	    x = event->xbutton.x;
	    y = event->xbutton.y;

	    /* update rotation matrix */

	    axis.x = y - lasty;
	    axis.y = lastx - x;
	    axis.z = 0;
	    len = length(axis);
	    /* must test length, to make sure that we have not
		   ended up at the last spot (by skipping over some
		   different intermediate position).  I think this is
		   not a bug causing problem for the other controls,
		   but keep it in mind...  */
	    if (len >= 0.9)
	      {
		radius = (cenx + ceny)/2; /* virtual trackball size */
		rotmat = matmult
		  (axis_angle(unit(axis), len/radius),
		   rotmat);
	      }
	    /* if (debug) printf("length = %f\n", length(axis)); */
	    lastx = x;
	    lasty = y;
	    draw_mol();
	  }
      }
      break;
    default:
      break;
    }
}


/* Button 2 Events in the pixmap are for rotation */
void xytranslate(Widget w, XtPointer client_data, XEvent * event, 
             Boolean *continue_to_dispatch) 
{
  static int lastx, lasty;
  int x, y;

  switch (event->type)
    {
    case ButtonPress:
      lastx = event->xbutton.x;
      lasty = event->xbutton.y;
      break;
    case ButtonRelease:
      draw_mol();
      break;
    case MotionNotify:
      {
	/* Don't let these events pile up */
	if (XEventsQueued(dpy, QueuedAfterReading) > 0) break;
	else
	  {
	    vector axis;

	    x = event->xbutton.x;
	    y = event->xbutton.y;

	    /* update translation vector */
	    axis.x = x - lastx;
	    axis.y = y - lasty;
	    axis.z = 0;
	    /* convert to Angstroms */
	    axis = vecscale(1/scale, axis);
	    /* convert to pre-rotation coordinates */
	    axis = matvmult(transpose(rotmat), axis);
	    transvec = vectoradd(transvec, axis);

	    lastx = x;
	    lasty = y;
	    draw_mol();
	  }
      }
      break;
    default:
      break;
    }
}


/* Button 2 Events in the pixmap are NO LONGER for Z-axis rotation */
void zrotate(Widget w, XtPointer client_data, XEvent * event,
              Boolean *continue_to_dispatch)
{
  static int lastx, lasty;
  int x, y;
  
  if (debug) printf("Button2Motion\n");
  switch (event->type)
    {
    case ButtonPress:
      lastx = event->xbutton.x;
      lasty = event->xbutton.y;
      break;
    case ButtonRelease:
      draw_mol();
      break;
    case MotionNotify:
      {
	/* Don't let these events pile up */
	if (XEventsQueued(dpy, QueuedAfterReading) > 0) break;
	else
	  {
	    vector axis, start, orth, move;
	    double len, rad, angle;

	    x = event->xbutton.x;
	    y = event->xbutton.y;

	    /* update rotation matrix */

	    /* rotation axis (Z, orthogonal to screen) */
	    axis.x = 0;
	    axis.y = 0;
	    axis.z = 1;

	    /* vector from center of screen to end mouse position */
	    /* Using the start position (lastx) seems to cause a brief
	       "backtrack" when the rotation direction is reversed */
	    start.x = x - cenx;
	    start.y = y - ceny;
	    start.z = 0;
	    rad = length(start);

	    /* unit vector orthogonal to "start" */
	    orth.x = start.y;
	    orth.y = -(start.x);
	    orth.z = 0;
	    orth = unit(orth);
	    
	    /* Pointer motion vector */
	    move.x = x - lastx;
	    move.y = y - lasty;
	    move.z = 0;

	    len = dot(move, orth);
	    /* len = PI * (y - lasty + lastx - x)/(double)width; /* */

	    /* don't let it blow up right at the center */
	    if (rad < 3) angle = 0;
	    else angle = len/rad;

	    rotmat = matmult
	      (axis_angle(axis, angle), rotmat); 
	  }
	lastx = x;
	lasty = y;
	draw_mol();
      }
      break;
    default:
      break;
    }
}


/* Button 3 Events in the pixmap are for zooming */
void zoom(Widget w, XtPointer client_data, XEvent * event,
              Boolean *continue_to_dispatch)
{
  static int lastx, lasty;
  int x, y;

  if (debug) printf("Btn3\n");
  switch (event->type)
    {
    case ButtonPress:
      lastx = event->xbutton.x;
      lasty = event->xbutton.y;
      break;
    case ButtonRelease:
      draw_mol();
      break;
    case MotionNotify:
      {
	/* Don't let these events pile up */
	if (XEventsQueued(dpy, QueuedAfterReading) > 0) break;
	else
	  {
	    x = event->xbutton.x;
	    y = event->xbutton.y;
	    scale *= 1 + ((float)(x - lastx + lasty - y))/200 ;

	    lastx = x;
	    lasty = y;
	    draw_mol();
	  }
      }
      break;
    }
}


void resize_strdraw(Widget w, XtPointer client_data, XtPointer call_data)
{
  /* these must be initialized, or they end up with funny
     values after the XtVaGetValues call */
  Dimension width = 0;
  Dimension height = 0;
  Dimension mheight = 50;
  Dimension margin_h = 0;
  Dimension margin_w = 0;
  Pixmap old = scratch;

  XtVaGetValues(str_mainw,
		XmNheight, &height,
		XmNwidth, &width,
		XmNmarginHeight, &margin_h,
		XmNmarginWidth, &margin_w,
		NULL);

  XtVaGetValues(menubar,
		XmNheight, &mheight,
		NULL);

  /* Don't crash with really small windows */
  if ((mheight + margin_h) >= height) return;
  if (margin_w >= width) return;

  height -= (mheight + margin_h);
  width -= margin_w;
  if ((height < 1) || (width < 1)) return;

  scratch = XCreatePixmap(dpy,
                         RootWindowOfScreen(XtScreen(str_mainw)),
                         width, height,
                         DefaultDepth(dpy, DefaultScreen(dpy)));


  if (scratch == XmUNSPECIFIED_PIXMAP)
    {
      printf("can't create pixmap\n");
      exit(1);
    }
  else
    {
      XmDestroyPixmap(XtScreen(str_mainw), old);
      XtVaSetValues(str_drawa,
		    XmNheight, height,
		    XmNwidth, width,
		    NULL);
    }
  cenx = width/2;
  ceny = height/2;
  draw_mol();
}


static void make_selection_string(int x, int y, char * message)
{
  uint i, j;
  vector transform;
  vector closest;
  float mindist, dist;
  matrix smat;

  /* cycle through all coordinates */
  smat = matscale(scale, rotmat);
  if (strctwin_nstrct < 1) sprintf(message, "No model\0");
  else
    {
      mindist = 5000;
      for (i = 0; i < strctwin_nstrct; ++i)
	{
	  for (j = 0; j < strctwin_size[i]; ++j)
	    {
	      /* translate, rotate, scale coordinate */
	      transform = vectoradd(strctwin_array[i][j], transvec); /* translate */
	      transform = matvmult(smat, transform); /* rotate, scale */
	      transform.x += cenx; /* center in window */
	      transform.y += ceny;
	      dist = (transform.x-x)*(transform.x-x) + 
		(transform.y-y)*(transform.y-y);
	      if (dist < mindist)
		{
		  closest = transform;
		  mindist = dist;
		  sprintf(message, "%s: %.1f %.1f %.1f (%d,%d)\0", 
			  label_array[i][j], closest.x, closest.y, closest.z, x, y);
		}
	    }
	}
    }
}
