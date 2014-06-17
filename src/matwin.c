#include <Xm/MainW.h>
#include <Xm/DrawingA.h>

/* $Id: */
/* $Header: */
/* $Log: */

extern Widget toplevel;

#define STARTWIDTH 200
#define STARTHEIGHT 150

static Widget mat_mainw, mat_shell, mat_drawa;
static Display * dpy;

static void file_cb(Widget widget, XtPointer client_data, XtPointer call_data);

void matwin_activate()
{
  dpy = XtDisplay(toplevel);
  if (!mat_shell)
    {
      XmString file, help; /* Menus */
      XmString loadalig, loadpath, loadddst, loaddrot, close; /* File Menu */
      Widget menu, menubar;

      mat_shell = XtVaAppCreateShell
	(NULL, 
	 "XSequoia", topLevelShellWidgetClass, dpy,
	 XtNtitle, "XSequoia Matrix Viewer",
	 /* Don't destroy window on deletion */
	 XmNdeleteResponse, XmUNMAP,
	 NULL);

      mat_mainw = XtVaCreateManagedWidget
	("matrix_window",
	 xmMainWindowWidgetClass, mat_shell,
	 XmNvisualPolicy, XmVARIABLE,
	 NULL);
      XtManageChild(mat_mainw);

      /* Menubar */
      file = XmStringCreateLocalized("Matrix");
      help = XmStringCreateLocalized("Help");
      menubar = XmVaCreateSimpleMenuBar(mat_mainw, "menubar",
					XmVaCASCADEBUTTON, file, 'F',
					XmVaCASCADEBUTTON, help, 'H',
					NULL);
      /* Tell which is the help button */
      if ((menu = XtNameToWidget(menubar, "button_1")))
	XtVaSetValues (menubar, XmNmenuHelpWidget, menu, NULL);
      XmStringFree(file);
      XmStringFree(help);
      /* Disable option for now */
      /* XtSetSensitive(XtNameToWidget(menubar, "button_0"), False); /* File */
      XtSetSensitive(XtNameToWidget(menubar, "button_1"), False); /* Help */
      XtManageChild(menubar);

      /* File menu */
      loadalig = XmStringCreateLocalized("Alignment");
      loadpath = XmStringCreateLocalized("Path Score");
      loadddst = XmStringCreateLocalized("Difference Distance");
      loaddrot = XmStringCreateLocalized("Difference Rotation");
      close = XmStringCreateLocalized("Close Window");
      menu = XmVaCreateSimplePulldownMenu(menubar, "file_menu", 0, file_cb,
                     XmVaPUSHBUTTON, loadalig, 'A', NULL, NULL,
                     XmVaPUSHBUTTON, loadpath, 'P', NULL, NULL, 
                     XmVaPUSHBUTTON, loadddst, 'D', NULL, NULL,
                     XmVaPUSHBUTTON, loaddrot, 'R', NULL, NULL,
                     XmVaSEPARATOR,
                     XmVaPUSHBUTTON, close, 'C', NULL, NULL,
                     NULL);
      XmStringFree(loadalig);
      XmStringFree(loadpath);
      XmStringFree(loadddst);
      XmStringFree(loaddrot);
      XmStringFree(close);
      /* Disable option for now */
      XtSetSensitive(XtNameToWidget(menu, "button_0"), False); 
      XtSetSensitive(XtNameToWidget(menu, "button_1"), False); 
      XtSetSensitive(XtNameToWidget(menu, "button_2"), False); 
      XtSetSensitive(XtNameToWidget(menu, "button_3"), False); 
      /* XtSetSensitive(XtNameToWidget(menu, "button_4"), False); /* Close */

      mat_drawa = XtVaCreateManagedWidget
	("mat_drawa",
	 xmDrawingAreaWidgetClass, mat_mainw,
	 XmNwidth, STARTWIDTH,
	 XmNheight, STARTHEIGHT,
	 NULL);
      XtManageChild(mat_drawa);
    }
  XtPopup(mat_shell, XtGrabNone);
}


/* File menu callback function */
void file_cb(Widget widget, XtPointer client_data, XtPointer call_data)
{
  int item_no = (int) client_data;
  static Widget fdialog;

  switch (item_no)
    {
    case 0:
      break;
    case 1: 
      break;
    case 2:
      break;
    case 3: 
      break;
    case 4:
      XtPopdown(mat_shell);
      break;
    default:
      break;
    }
}
