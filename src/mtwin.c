#include "mtwin.h"

/* $Id: mtwin.c,v 1.2 2001/11/28 23:40:22 bruns Exp $ */
/* $Header: /usr/data/cvs/sequoia1/mtwin.c,v 1.2 2001/11/28 23:40:22 bruns Exp $ */
/* $Log: mtwin.c,v $
/* Revision 1.2  2001/11/28 23:40:22  bruns
/* Added cvs header tags
/* Removed ^M characters
/* */

extern void exec_cmd(Widget widget, XtPointer client_data, XtPointer call_data);
extern void parseit(char * cmd);

/* variable used by other source files */
XtAppContext app;
Widget main_menubar, text_w, command_w, toplevel;
XmTextPosition textpos = 0;

void init_windows(int argc, char *argv[]);

/* global variables private to this file */
static Widget main_w;
static char read_reg[100] = "SEQ1";
static GC gc;

static void init_menus(void);
static void file_cb(Widget widget, XtPointer client_data, XtPointer call_data);
static void edit_cb(Widget widget, XtPointer client_data, XtPointer call_data);
static void action_cb(Widget widget, XtPointer client_data, XtPointer call_data);
static void view_cb(Widget widget, XtPointer client_data, XtPointer call_data);
static void settings_cb(Widget widget, XtPointer client_data, XtPointer call_data);
static void option_cb(Widget widget, XtPointer client_data, XtPointer call_data);
static void load_seq(Widget dialog, XtPointer client_data, XtPointer call_data);
static void parse_and_history(char * cmd);

void init_windows(int argc, char *argv[])
{
  Display * dpy;
  Dimension width = 300;
  Dimension height = 350;
  short margin_w = 0;
  short margin_h = 0;
  Dimension mheight = 0; /* menu height */
  XmString file;

  /* scrolled text */
  Arg args[5];
  int n = 0;
  /* end of variables */

  XtSetLanguageProc(NULL,NULL,NULL);

  /* Main window stuff */
  toplevel = XtVaAppInitialize(&app, "XSequoia", NULL, 0, &argc,
                               argv, NULL, NULL);

  XtVaSetValues(toplevel, 
		XtNtitle, "XSequoia by Chris Bruns",
		NULL);
  
  main_w = XtVaCreateManagedWidget ("main_window",
                     xmMainWindowWidgetClass, toplevel,
                     XmNvisualPolicy, XmVARIABLE,
                     XmNcommandWindowLocation, XmCOMMAND_BELOW_WORKSPACE,
                     NULL);

  dpy = XtDisplay(main_w);
  gc = DefaultGC(dpy, DefaultScreen(dpy));

  XtVaGetValues(main_w,
                XmNheight, &height,
                XmNwidth, &width,
                XmNmarginHeight, &margin_h,
                XmNmarginWidth, &margin_w,
                NULL);

  /* Command area */
  file = XmStringCreateLocalized("Command:");
  command_w = XtVaCreateWidget("command_w", xmCommandWidgetClass, main_w,
			       XmNpromptString, file,
			       NULL);
  XmStringFree(file);
  XtAddCallback(command_w, XmNcommandEnteredCallback, exec_cmd, text_w);
  XtManageChild(command_w);

  init_menus();

  XtVaGetValues(main_menubar,
                XmNheight, &mheight,
                NULL);

  height -= (mheight + margin_h);
  width -= margin_w;

  /* Create ScrolledText in main window work area */
  XtSetArg (args[n], XmNrows, 24); n++;
  XtSetArg (args[n], XmNcolumns, 80); n++;
  XtSetArg (args[n], XmNeditable, False); n++;
  XtSetArg (args[n], XmNeditMode, XmMULTI_LINE_EDIT); n++;
  text_w = XmCreateScrolledText (main_w, "text_w", args, n);
  XtManageChild(text_w);

  XmMainWindowSetAreas(main_w, main_menubar, command_w,
                       NULL, NULL, XtParent(text_w));

  XtVaSetValues(main_menubar, XmNuserData, text_w, NULL);

  XtVaSetValues(main_w,
                XmNmenuBar, main_menubar,
		XmNinitialFocus, command_w,
                NULL);

  XtRealizeWidget(toplevel);
}

void init_menus(void)
{
  /* Menu Bar variables */
  XmString file, edit, action, view, settings, help; /* Menus */
  XmString open, save, quit, quit_acc; /* File */
  XmString cut, copy, paste; /* Edit */
  XmString align, overlay, runscript, consensus, statistics; /* Action */
  XmString sequences, coordinates, matrices; /* View */
  XmString alignment; /* Settings */
  Widget menu;

  /* Main window variables */
  Widget widget;

  /* Menu bar stuff */
  file = XmStringCreateLocalized("File");
  edit = XmStringCreateLocalized("Edit");
  action = XmStringCreateLocalized("Action");
  view = XmStringCreateLocalized("View");
  settings = XmStringCreateLocalized("Settings");
  help = XmStringCreateLocalized("Help");
  main_menubar = XmVaCreateSimpleMenuBar(main_w, "main_menubar",
                                    XmVaCASCADEBUTTON, file, 'F',
                                    XmVaCASCADEBUTTON, edit, 'E',
                                    XmVaCASCADEBUTTON, action, 'A',
                                    XmVaCASCADEBUTTON, view, 'V',
                                    XmVaCASCADEBUTTON, settings, 'S',
                                    XmVaCASCADEBUTTON, help, 'H',
                                    NULL);
  XmStringFree(file);
  XmStringFree(edit);
  XmStringFree(action);
  XmStringFree(view);
  XmStringFree(settings);
  XmStringFree(help);

  /* Tell which is the help button */
  if ((widget = XtNameToWidget(main_menubar, "button_5")))
    XtVaSetValues (main_menubar, XmNmenuHelpWidget, widget, NULL);

  /* File Menu */
  open = XmStringCreateLocalized("Open File...");
  save = XmStringCreateLocalized("Save File...");
  quit = XmStringCreateLocalized("Quit");
  quit_acc = XmStringCreateLocalized("Ctrl-C");
  menu = XmVaCreateSimplePulldownMenu(main_menubar, "file_menu", 0, file_cb,
                     XmVaPUSHBUTTON, open, 'O', NULL, NULL,
                     XmVaPUSHBUTTON, save, 'S', NULL, NULL,
                     XmVaSEPARATOR,
                     XmVaPUSHBUTTON, quit, 'Q', "Ctrl<Key>c", quit_acc,
                     NULL);
  XmStringFree(open);
  XmStringFree(save);
  XmStringFree(quit);
  XmStringFree(quit_acc);
  /* Disable option for now */
  /* XtSetSensitive(XtNameToWidget(menu, "button_0"), False); /* Open */
  XtSetSensitive(XtNameToWidget(menu, "button_1"), False); /* Save */

  /* Edit Menu */
  cut = XmStringCreateLocalized("Cut");
  copy = XmStringCreateLocalized("Copy");
  paste = XmStringCreateLocalized("Paste");
  menu = XmVaCreateSimplePulldownMenu(main_menubar, "edit_menu", 1, edit_cb,
                     XmVaPUSHBUTTON, cut, 'X', NULL, NULL,
                     XmVaPUSHBUTTON, copy, 'C', NULL, NULL,
                     XmVaPUSHBUTTON, paste, 'P', NULL, NULL,
                     NULL);
  XmStringFree(cut);
  XmStringFree(copy);
  XmStringFree(paste);
  XtSetSensitive(XtNameToWidget(menu, "button_0"), False); /* Cut */
  XtSetSensitive(XtNameToWidget(menu, "button_1"), False); /* Copy */
  XtSetSensitive(XtNameToWidget(menu, "button_2"), False); /* Paste */

  /* Action Menu */
  /* XmString align, overlay, runscript, statistics; /* Action */
  align = XmStringCreateLocalized("Align");
  overlay = XmStringCreateLocalized("Overlay...");
  runscript = XmStringCreateLocalized("Run Script...");
  consensus = XmStringCreateLocalized("Consensus");
  statistics = XmStringCreateLocalized("Statistics...");
  menu = XmVaCreateSimplePulldownMenu(main_menubar, "action_menu", 2, action_cb,
                     XmVaPUSHBUTTON, align, 'A', NULL, NULL,
                     XmVaPUSHBUTTON, overlay, 'O', NULL, NULL,
                     XmVaPUSHBUTTON, runscript, 'R', NULL, NULL,
                     XmVaPUSHBUTTON, consensus, 'C', NULL, NULL,
                     XmVaPUSHBUTTON, statistics, 'S', NULL, NULL,
                     NULL);
  XmStringFree(align);
  XmStringFree(overlay);
  XmStringFree(runscript);
  XmStringFree(consensus);
  XmStringFree(statistics);
  /* XtSetSensitive(XtNameToWidget(menu, "button_0"), False); /* Align */
  XtSetSensitive(XtNameToWidget(menu, "button_1"), False); /* Overlay */
  XtSetSensitive(XtNameToWidget(menu, "button_2"), False); /* Run Script */
  XtSetSensitive(XtNameToWidget(menu, "button_3"), False); /* Consensus */
  XtSetSensitive(XtNameToWidget(menu, "button_4"), False); /* Statistics */

  /* View Menu */
  sequences = XmStringCreateLocalized("Sequences...");
  coordinates = XmStringCreateLocalized("Coordinates...");
  matrices = XmStringCreateLocalized("Matrix...");
  menu = XmVaCreateSimplePulldownMenu(main_menubar, "view_menu", 3, view_cb,
                     XmVaPUSHBUTTON, sequences, 'S', NULL, NULL,
                     XmVaPUSHBUTTON, coordinates, 'C', NULL, NULL,
                     XmVaPUSHBUTTON, matrices, 'M', NULL, NULL,
                     NULL);
  XmStringFree(sequences);
  XmStringFree(coordinates);
  XmStringFree(matrices);
  /* XtSetSensitive(XtNameToWidget(menu, "button_0"), False); /* Sequences */
  /* XtSetSensitive(XtNameToWidget(menu, "button_1"), False); /* Coordinates */
  /* XtSetSensitive(XtNameToWidget(menu, "button_2"), False); /* Matrices */

  /* Settings Menu */
  alignment = XmStringCreateLocalized("Alignment...");
  menu = XmVaCreateSimplePulldownMenu(main_menubar, "settings_menu", 4, settings_cb,
                     XmVaPUSHBUTTON, alignment, 'A', NULL, NULL,
                     NULL);
  XmStringFree(alignment);
  XtSetSensitive(XtNameToWidget(menu, "button_0"), False); /* Alignment */

  XtManageChild(main_menubar);
}

/* File menu callback function */
void file_cb(Widget widget, XtPointer client_data, XtPointer call_data)
{
  int item_no = (int) client_data;
  static Widget fdialog;

  switch (item_no)
    {
    case 0: /* Open file */
      if (!fdialog)
        {
	  XmString seqone, seqtwo, structone, structtwo, align, overlay;
	  XmString reg, okbutton, cancelbutton;
	  Widget rc, option_menu;

          fdialog = XmCreateFileSelectionDialog (toplevel, "file_sel",
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
	  structone = XmStringCreateLocalized("STRUCT1");
	  structtwo = XmStringCreateLocalized("STRUCT2");
	  align = XmStringCreateLocalized("ALIGN");
	  overlay = XmStringCreateLocalized("OVERLAY");
	  reg = XmStringCreateLocalized("Register:");
	  option_menu = XmVaCreateSimpleOptionMenu (rc, "option_menu",
			  reg, 'R', 0 /* initial selection */, option_cb,
                          XmVaPUSHBUTTON, seqone, 'S', NULL, NULL,
                          XmVaPUSHBUTTON, seqtwo, 'Q', NULL, NULL,
                          XmVaPUSHBUTTON, align, 'A', NULL, NULL,
                          XmVaPUSHBUTTON, structone, 'T', NULL, NULL,
                          XmVaPUSHBUTTON, structtwo, 'R', NULL, NULL,
                          XmVaPUSHBUTTON, overlay, 'O', NULL, NULL,
                          NULL);
	  XtManageChild(option_menu);
          XtAddCallback(fdialog, XmNokCallback, load_seq, NULL);
          XtAddCallback(fdialog, XmNcancelCallback,
                        (void *) XtUnmanageChild, NULL);
          /* No help exists, so disable */
          /* XtSetSensitive((Widget) XmSelectionBoxGetChild
                         (fdialog, XmDIALOG_HELP_BUTTON), False); /* */
	  XtUnmanageChild(XmFileSelectionBoxGetChild(fdialog, XmDIALOG_HELP_BUTTON));
	  XmStringFree(seqone);
	  XmStringFree(seqtwo);
	  XmStringFree(structone);
	  XmStringFree(structtwo);
	  XmStringFree(align);
	  XmStringFree(overlay);
        }
      /* make the dialog close automatically when done */
      /* XtVaSetValues (fdialog,
                     XmNautoUnmanage, True,
                     NULL); /* */
      XtManageChild(fdialog);
      XtPopup(XtParent(fdialog), XtGrabNone);
      break;
    case 2:
      exit(0); /* Quit */
      break;
    }
}

/* Edit menu callback function */
void edit_cb(Widget widget, XtPointer client_data, XtPointer call_data)
{
  return;
}

/* Action menu callback function */
void action_cb(Widget widget, XtPointer client_data, XtPointer call_data)
{
  int item_no = (int) client_data;
  char cmd[100];
  
  switch (item_no)
    {
    case 0: /* Align */
      strcpy(cmd, "ALIGN");
      parse_and_history(cmd);
      break;
    case 1:
      break;
    default:
      break;
    }
  return;
}

/* View menu callback function */
void view_cb(Widget widget, XtPointer client_data, XtPointer call_data)
{
  int item_no = (int) client_data;
  extern void seqwin_activate();
  extern void strwin_activate();
  extern void matwin_activate();

  switch (item_no)
    {
    case 0: /* Sequence editor window */
      seqwin_activate();
      break;
    case 1: /* PDB file viewer */
      strwin_activate();
      break;
    case 2: /* Matrix */
      matwin_activate();
      break;
    default:
      break;
    }
  return;
}

/* Settings menu callback function */
void settings_cb(Widget widget, XtPointer client_data, XtPointer call_data)
{
  return;
}

/* File input chooser option menu callback function */
void option_cb(Widget widget, XtPointer client_data, XtPointer call_data)
{
  int item_no = (int) client_data;
  if (item_no == 0) strcpy(read_reg, "SEQ1");
  else if (item_no == 1) strcpy(read_reg, "SEQ2");
  else if (item_no == 2) strcpy(read_reg, "ALIGN");
  else if (item_no == 3) strcpy(read_reg, "STRUCT1");
  else if (item_no == 4) strcpy(read_reg, "STRUCT2");
  else if (item_no == 5) strcpy(read_reg, "OVERLAY");
  else strcpy(read_reg, "SEQ1");
  return;
}

void load_seq(Widget widget, XtPointer client_data, XtPointer call_data)
{
  char cmd[2000];
  char fname[1000];
  char * file = NULL;
  XmFileSelectionBoxCallbackStruct *cbs =
    (XmFileSelectionBoxCallbackStruct *) call_data;
  XmString command;
  Widget history_list;
  int pos; /* n+1 position of history list */

  if (cbs)
    {
      if (!XmStringGetLtoR (cbs->value, XmFONTLIST_DEFAULT_TAG, &file))
        return; /* internal error */
      (void) strcpy(fname, file);
      XtFree(file);
    }
  else return;

  sprintf(cmd, "READ %s %s", read_reg, fname);

  /* append this command to the history list */
  command = XmStringCreateLocalized(cmd);
  history_list = XmCommandGetChild(command_w, XmDIALOG_HISTORY_LIST);
  XtVaGetValues(history_list,
	XmNitemCount, &pos,
        NULL);
  ++pos; /* add it AFTER the last item */
  XmListAddItemUnselected(history_list, command, pos);
  XmStringFree(command);

  /* actually execute the command */
  parseit(cmd);
}

void parse_and_history(char * cmd)
{
  int pos;
  XmString command;
  Widget history_list;

  /* append this command to the history list */
  command = XmStringCreateLocalized(cmd);
  history_list = XmCommandGetChild(command_w, XmDIALOG_HISTORY_LIST);
  XtVaGetValues(history_list,
		XmNitemCount, &pos,
		NULL);
  ++pos; /* add it AFTER the last item */
  XmListAddItemUnselected(history_list, command, pos);
  XmStringFree(command);
  /* actually execute the command */
  parseit(cmd);
}
