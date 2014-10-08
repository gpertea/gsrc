#ifndef CLVIEWMAINWIN_H
#define CLVIEWMAINWIN_H
#define COL_RAMP_SIZE 16
#define NUM_GROUPS 24

#include "fx.h"

class FXSplashDlg;
// Main Window
class MainWindow : public FXMainWindow {
  FXDECLARE(MainWindow)
protected:
  // Member data
  FXSplashDlg*  splashmsg;
  FXStatusBar*  statusbar;
  FXMenuBar*    menubar;
  FXMenuPane*   filemenu;
  FXMenuPane*   helpmenu;
  FXMenuPane*   windowmenu;
  FXMDIClient*  mdiclient;               // MDI Client area
  FXMDIMenu*    mdimenu;                 // MDI Window Menu
  FXGIFIcon*    mdiicon;                 // MDI Icon
  FXComboBox* filter;
  FXDirList*  dirlist;
  bool constructing;
//color tables:
  FXColor baseColors[6]; //base colors
  FXColor seqColors[12];
      /*
              buffer for color array to pass to
               the Seq Paint method
           0 = seq color
           1 = seq highlight
           2 = seq shadow
           3 = seqtext color
           4 = seqtrim color
           5 = seq range color
           6 = seq select color
           7 = alternate seq color (ET)
           8 = inter-segment connector color
           9 = internal trim color
          10 = major consensus splice site
          11 = minor consensus splice site
       */

  FXColor matchColors[COL_RAMP_SIZE+1]; //read luminance ramp
  FXColor ctgQuals[COL_RAMP_SIZE+1]; //contig luminance ramp
  FXColor grpColors[2][NUM_GROUPS+1]; //group colors
protected:
  void initColors();
  MainWindow(){}
public:
  // Message handlers
  long onCmdAbout(FXObject*,FXSelector,void*);
  long onRandomLayout(FXObject*,FXSelector,void*);
  long onCmdFilter(FXObject*,FXSelector,void* ptr);
  long onCmdFileDlg(FXObject*,FXSelector,void* ptr);
  long onCmdTreeList(FXObject*,FXSelector,void* ptr);
  long onMaxRestore(FXObject* obj,FXSelector sel,void* ptr);
  void splash(const char* message=NULL);
public:
  void setCurrentPattern(FXint n);
  void newMDIChild(const char* filename=NULL);
  FXString getPatterns() const;
  FXint getCurrentPattern() {return filter->getCurrentItem();}
  void setPatterns(const FXString& patterns);
  FXString getPatterns();
  // Messages
 enum {
    ID_ABOUT=FXMainWindow::ID_LAST,
    ID_FILEFILTER,
    ID_TREELIST,
    ID_FILEOPEN,
    ID_RANDOM,
    ID_MAXIMIZE,
    ID_CLOSEALL,
    ID_LAST
    };
 public:
  MainWindow(FXApp* a);
  long onCmdCloseAll(FXObject*,FXSelector,void*);
  virtual void create();
  virtual ~MainWindow();
 };

/*******************************************
 splash/dialog to use for various purposes
********************************************/
class FXSplashDlg : public FXDialogBox {
  FXDECLARE(FXSplashDlg)
protected:
  FXVerticalFrame* contents;
  FXLabel *msg;
private:
  FXSplashDlg(){}
public:
  FXSplashDlg(FXWindow* owner, const char* message);
  virtual ~FXSplashDlg();
  };

#endif
