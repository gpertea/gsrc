/********************************************************************************
*                                                                               *
*                                                                               *
*********************************************************************************/
#include "mainwin.h"
#include "mdichild.h"
#include "clrutils.h"
/*******************************************************************************/

// Map
FXDEFMAP(MainWindow) MainWindowMap[]={
  //________Message_Type____________ID___________________Message_Handler________
  FXMAPFUNC(SEL_COMMAND,  MainWindow::ID_ABOUT,            MainWindow::onCmdAbout),
  FXMAPFUNC(SEL_COMMAND,  MainWindow::ID_RANDOM,           MainWindow::onRandomLayout),
  FXMAPFUNC(SEL_COMMAND,  MainWindow::ID_FILEFILTER,       MainWindow::onCmdFilter),
  FXMAPFUNC(SEL_COMMAND,  MainWindow::ID_FILEOPEN,         MainWindow::onCmdFileDlg),
  //FXMAPFUNC(SEL_SIGNAL,   MainWindow::ID_CLOSEALL,         MainWindow::onCmdCloseAll),
  FXMAPFUNC(SEL_COMMAND,  MainWindow::ID_CLOSEALL,         MainWindow::onCmdCloseAll),
  FXMAPFUNC(SEL_COMMAND,  MDIChild::ID_MAXRESTORE,         MainWindow::onMaxRestore),
  FXMAPFUNC(SEL_DOUBLECLICKED,  MainWindow::ID_TREELIST, MainWindow::onCmdTreeList),
  };

// Object implementation
FXIMPLEMENT(MainWindow,FXMainWindow,MainWindowMap,ARRAYNUMBER(MainWindowMap))


// FXSplashDlg implementation
FXIMPLEMENT(FXSplashDlg,FXDialogBox,NULL,0)

/*******************************************************************************/

const unsigned char penguin[]={
  0x47,0x49,0x46,0x38,0x37,0x61,0x10,0x00,0x12,0x00,0xf2,0x00,0x00,0xb2,0xc0,0xdc,
  0x80,0x80,0x80,0x00,0x00,0x00,0xc0,0xc0,0xc0,0x10,0x10,0x10,0xff,0xff,0xff,0xe0,
  0xa0,0x08,0x00,0x00,0x00,0x2c,0x00,0x00,0x00,0x00,0x10,0x00,0x12,0x00,0x00,0x03,
  0x53,0x08,0xba,0x21,0x12,0x2b,0xc6,0xe6,0x9e,0x94,0x62,0x64,0x77,0xa3,0x20,0x4e,
  0x21,0x74,0x8b,0x60,0x9c,0x1a,0xa9,0x98,0xa8,0x45,0xb2,0x85,0x38,0x76,0x4f,0x6c,
  0xbb,0x93,0x60,0xdb,0x0d,0xe4,0xd9,0x83,0x1d,0xe7,0x57,0x18,0x04,0x6f,0xb8,0x4c,
  0xec,0x88,0x9c,0x01,0x0c,0x47,0x66,0xac,0xa2,0x38,0x19,0x76,0x36,0x83,0xc3,0xf0,
  0xb4,0x5e,0x77,0x03,0xaf,0xf8,0x7b,0x13,0x77,0xad,0xd3,0xad,0x75,0x61,0xa5,0x54,
  0x02,0x27,0x45,0x02,0x00,0x3b
  };


static const FXchar fltfiles[]="ACE files (*.ace,*ACE)\nLayout files (*.lyt)\nAny file (*)";

FXString fileNameOnly(const FXString& file){
  register FXint f,n;
  if(!file.empty()){
    n=0;
#ifdef WIN32
    if(isalpha((FXuchar)file[0]) && file[1]==':') n=2;
#endif
    f=n;
    while(file[n]){
      if(ISPATHSEP(file[n])) f=n+1;
      n++;
      }
    return FXString(file.text()+f,n-f);
    }
  return FXString::null;
  }

void MainWindow::initColors() {
  FXRegistry& ini=getApp()->reg();
  static const char* section="Colors";
  static const char* bsection="Base Color Mode";
  static const char* grpsect="Group Colors";

  matchColors[16]=ini.readColorEntry(section, "Mismatch",
          makeHiliteColor(makeHiliteColor(colorFromName("Red"))));
  ini.writeColorEntry(section,"Mismatch",matchColors[16]);
  matchColors[15]=ini.readColorEntry(section, "MaxCoverage",
               colorFromName("darkBlue"));
  ini.writeColorEntry(section,"MaxCoverage",matchColors[15]);
  ctgQuals[16]=ini.readColorEntry(section, "ContigNoCoverage",
                                                 getApp()->getBaseColor());
  ini.writeColorEntry(section,"ContigNoCoverage",ctgQuals[16]);
  ctgQuals[15]=ini.readColorEntry(section,"ContigMaxCoverage",0xFF00FFFF);
  ini.writeColorEntry(section,"ContigMaxCoverage",ctgQuals[15]);

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
  //-- base sequence color:

  seqColors[0]=ini.readColorEntry(section,"SeqBg",colorFromName("darkBlue"));
  // Hex Colors are BGR style !
  //-- highlight border:
  seqColors[1]=ini.readColorEntry(section,"SeqBorderHilight",
                                 makeHiliteColor(colorFromName("Gray50")));
  //--shadow border:
  seqColors[2]=ini.readColorEntry(section,"SeqBorderShadow",
                                 makeShadowColor(colorFromName("Gray20")));
  //nucleotide font color:
  seqColors[3]=ini.readColorEntry(section,"SeqText",colorFromName("White"));
  //-- clipping color:
  seqColors[4]=ini.readColorEntry(section,"SeqClipping",
                       colorFromName("MistyRose3")); //trimmed ends
  //seq range color:
  seqColors[5]=ini.readColorEntry(section,"SeqRange",colorFromName("Green"));
  //seq select color:
  seqColors[6]=ini.readColorEntry(section,"SeqSelect", colorFromName("Yellow"));
  //-- NP/ET color:
  seqColors[7]=ini.readColorEntry(section,"SeqSpecial", colorFromName("darkRed"));
  //-- inter-segment gap line:
  seqColors[8]=ini.readColorEntry(section, "Intron", colorFromName("LightSkyBlue3"));
  //-- internal clipping color:
  seqColors[9]=ini.readColorEntry(section, "InternalClipping", colorFromName("MistyRose1"));
  //-- splice site consensus
  seqColors[10]=ini.readColorEntry(section, "MajorSplice", colorFromName("Red")); //hard
  seqColors[11]=ini.readColorEntry(section, "MinorSplice", colorFromName("DarkPink")); //soft
  ini.writeColorEntry(section, "SeqBg",seqColors[0]);
  ini.writeColorEntry(section, "SeqBorderHilight",seqColors[1]);
  ini.writeColorEntry(section, "SeqBorderShadow",seqColors[2]);
  ini.writeColorEntry(section, "SeqText",seqColors[3]);
  ini.writeColorEntry(section, "SeqClipping",seqColors[4]);
  ini.writeColorEntry(section, "SeqRange",seqColors[5]);
  ini.writeColorEntry(section, "SeqSelect",seqColors[6]);
  ini.writeColorEntry(section, "SeqSpecial",seqColors[7]);
  ini.writeColorEntry(section, "Intron",seqColors[8]);
  ini.writeColorEntry(section, "InternalClipping",seqColors[9]);
  ini.writeColorEntry(section, "MajorSplice",seqColors[10]);
  ini.writeColorEntry(section, "MinorSplice",seqColors[11]);

  FXString s;
  FXColor ggapc;
  for (int i=0;i<NUM_GROUPS;i++) {
    s.format("ColorGrp%d",i);
    grpColors[0][i]=ini.readColorEntry(grpsect, s.text(),
                              getRColor(seqColors[0],NUM_GROUPS,i));
    ini.writeColorEntry(grpsect,s.text(),grpColors[0][i]);
    //intron color: increase luminance, decrease saturation of the group color
    ggapc=makeHiliteColor(makeHiliteColor(grpColors[0][i]));
    ggapc=modhls(ggapc,0,45,-200);
    if (ggapc==backColor) { //white, I assume..
        ggapc=modhls(ggapc, 0,-10,10);
        }
    grpColors[1][i]=ggapc;
    }
  for (int i=0;i<15;i++) {
    matchColors[14-i]=modhls(matchColors[15], 0, 10*i, -5*i);
    ctgQuals[14-i]=modhls(ctgQuals[15], 0, -4*i, 0);
    }
  baseColors[0]=ini.readColorEntry(bsection, "BaseText", colorFromName("White")); //font color
  baseColors[1]=ini.readColorEntry(bsection, "BaseA", colorFromName("DarkGreen"));  //A
  baseColors[2]=ini.readColorEntry(bsection, "BaseC", colorFromName("Blue"));   //C
  baseColors[3]=ini.readColorEntry(bsection, "BaseG", colorFromName("Orange")); //G
  baseColors[4]=ini.readColorEntry(bsection, "BaseT", colorFromName("Red"));    //T
  baseColors[5]=ini.readColorEntry(bsection, "BaseX", colorFromName("Gray20"));  // N,*, other codes
  ini.writeColorEntry(bsection, "BaseText",baseColors[0]);
  ini.writeColorEntry(bsection, "BaseA",baseColors[1]);
  ini.writeColorEntry(bsection, "BaseC",baseColors[2]);
  ini.writeColorEntry(bsection, "BaseG",baseColors[3]);
  ini.writeColorEntry(bsection, "BaseT",baseColors[4]);
  ini.writeColorEntry(bsection, "BaseX",baseColors[5]);
  ini.write();
}

// Make some windows
MainWindow::MainWindow(FXApp *a):FXMainWindow(a,"Cluster viewer",NULL,NULL,
   DECOR_ALL|PLACEMENT_SCREEN,40,40,760,560) {

  //MDIChild::ifont=NULL;
  constructing=true;
  initColors();
  new FXToolTip(a);
  menubar=new FXMenuBar(this,LAYOUT_SIDE_TOP|LAYOUT_FILL_X);
  // Status bar
  statusbar=new FXStatusBar(this,LAYOUT_SIDE_BOTTOM|LAYOUT_FILL_X|STATUSBAR_WITH_DRAGCORNER);
  FXSplitter* splitter=new FXSplitter(this,LAYOUT_SIDE_TOP|LAYOUT_FILL_X|LAYOUT_FILL_Y|SPLITTER_TRACKING);
  /*
  // Contents
  FXHorizontalFrame *contents = new FXHorizontalFrame(this,
          LAYOUT_SIDE_TOP|LAYOUT_FILL_X|LAYOUT_FILL_Y, 0,0,0,0, 0,0,0,0, 0,0);
  */

  // Sunken border for tree
   FXVerticalFrame* treebox=new FXVerticalFrame(splitter,LAYOUT_FILL_X|LAYOUT_FILL_Y,0,0,160,0,0,0,0,0);
   // Make tree
   FXHorizontalFrame* treeframe=new FXHorizontalFrame(treebox,FRAME_SUNKEN|FRAME_THICK|LAYOUT_FILL_X|LAYOUT_FILL_Y,0,0,0,0, 0,0,0,0);
   dirlist=new FXDirList(treeframe,this,ID_TREELIST,
      DIRLIST_SHOWFILES|TREELIST_BROWSESELECT|
      TREELIST_SHOWS_LINES|TREELIST_SHOWS_BOXES|LAYOUT_FILL_X|LAYOUT_FILL_Y);
   FXHorizontalFrame* filterframe=new FXHorizontalFrame(treebox,LAYOUT_FILL_X);
   new FXLabel(filterframe,"Filter:");
   filter=new FXComboBox(filterframe,25,this,ID_FILEFILTER,COMBOBOX_STATIC|LAYOUT_FILL_X|FRAME_SUNKEN|FRAME_THICK);
   filter->appendItem("ACE files (*.ace)");
   filter->appendItem("Layout files (*.lyt)");
   filter->appendItem("All Files (*)");
   filter->setNumVisible(3);
   setCurrentPattern(0);
   dirlist->setDirectory(FXSystem::getCurrentDirectory());
   treebox->hide();

  // --------- MDI Stuff  ------------------------------------------
  // Nice sunken box around the MDIClient
   FXVerticalFrame *mdibox=new FXVerticalFrame(splitter,
          FRAME_SUNKEN|LAYOUT_FILL_X|LAYOUT_FILL_Y,0,0,0,0, 0,0,0,0);
    // MDI Client
   mdiclient=new FXMDIClient(mdibox,LAYOUT_FILL_X|LAYOUT_FILL_Y);

   // Make MDI Menu
   mdimenu=new FXMDIMenu(this,mdiclient);
   mdiicon=new FXGIFIcon(a,penguin);

   // MDI buttons in menu:- note the message ID's!!!!!
   // Normally, MDI commands are simply sensitized or desensitized;
   // Under the menubar, however, they're hidden if the MDI Client is
   // not maximized.  To do this, they must have different ID's.
   new FXMDIWindowButton(menubar,mdimenu,mdiclient,FXMDIClient::ID_MDI_MENUWINDOW,LAYOUT_LEFT);
   new FXMDIDeleteButton(menubar,mdiclient,FXMDIClient::ID_MDI_MENUCLOSE,FRAME_RAISED|LAYOUT_RIGHT);
   new FXMDIRestoreButton(menubar,mdiclient,FXMDIClient::ID_MDI_MENURESTORE,FRAME_RAISED|LAYOUT_RIGHT);
   new FXMDIMinimizeButton(menubar,mdiclient,FXMDIClient::ID_MDI_MENUMINIMIZE,FRAME_RAISED|LAYOUT_RIGHT);

   //newMDIChild();
  //--------------------------------------------------------------------------------------------

  // File Menu
  filemenu=new FXMenuPane(this);
  new FXMenuCommand(filemenu,"File Browser\t\tDisplay file list.",NULL,treebox,
                                FXWindow::ID_TOGGLESHOWN);
  //new FXMenuCommand(filemenu,"&New\tCtl-N\tCreate new empty view",NULL,this,ID_NEW);
  new FXMenuCommand(filemenu,"&Open\tCtl-O\tOpen ace/layout file ", NULL, this, MainWindow::ID_FILEOPEN);
  new FXMenuSeparator(filemenu);
  /*
  new FXMenuCommand(filemenu,"&Random layout\t\tGenerate a random test layout.",NULL,this,
               MainWindow::ID_RANDOM,0);
  new FXMenuSeparator(filemenu);*/
  new FXMenuCommand(filemenu,"&Quit\tCtl-Q\tQuit application.",NULL,a,FXApp::ID_QUIT,0);
  new FXMenuTitle(menubar,"&File",NULL,filemenu);
  // Window menu
  windowmenu=new FXMenuPane(this);
  new FXMenuCommand(windowmenu,"Tile &Horizontally",NULL,mdiclient,FXMDIClient::ID_MDI_TILEHORIZONTAL);
  new FXMenuCommand(windowmenu,"Tile &Vertically",NULL,mdiclient,FXMDIClient::ID_MDI_TILEVERTICAL);
  new FXMenuCommand(windowmenu,"C&ascade",NULL,mdiclient,FXMDIClient::ID_MDI_CASCADE);
  new FXMenuCommand(windowmenu,"&Close",NULL,mdiclient,FXMDIClient::ID_MDI_CLOSE);
  new FXMenuCommand(windowmenu,"Close &All",NULL,mdiclient,MainWindow::ID_CLOSEALL);
  FXMenuSeparator* sep1=new FXMenuSeparator(windowmenu);
  sep1->setTarget(mdiclient);
  sep1->setSelector(FXMDIClient::ID_MDI_ANY);
  new FXMenuCommand(windowmenu,FXString::null,NULL,mdiclient,FXMDIClient::ID_MDI_1);
  new FXMenuCommand(windowmenu,FXString::null,NULL,mdiclient,FXMDIClient::ID_MDI_2);
  new FXMenuCommand(windowmenu,FXString::null,NULL,mdiclient,FXMDIClient::ID_MDI_3);
  new FXMenuCommand(windowmenu,FXString::null,NULL,mdiclient,FXMDIClient::ID_MDI_4);
  new FXMenuTitle(menubar,"&Window",NULL,windowmenu);
  // Help menu
  helpmenu=new FXMenuPane(this);
  new FXMenuCommand(helpmenu,"&About clView...",NULL,this,ID_ABOUT,0);
  new FXMenuTitle(menubar,"&Help",NULL,helpmenu,LAYOUT_RIGHT);
  constructing=false;
  a->setAnimSpeed(0);
  }


// Set current pattern
// Strip pattern from text if present
static FXString patternFromText(const FXString& pattern){
  FXint beg,end;
  end=pattern.rfind(')');
  beg=pattern.rfind('(',end-1);
  if(0<=beg && beg<end) return pattern.mid(beg+1,end-beg-1);
  return pattern;
  }

// Close all windows
long MainWindow::onCmdCloseAll(FXObject* sender, FXSelector, void*){
  mdiclient->forallWindows(sender, FXSEL(SEL_CLOSE, 0), NULL);
  return 1;
  }



void MainWindow::setCurrentPattern(FXint n){
  n=FXCLAMP(0,n,filter->getNumItems()-1);
  FXString cdir=dirlist->getDirectory();
  if (cdir.empty())
      cdir=FXSystem::getCurrentDirectory();
  filter->setCurrentItem(n);
  dirlist->setPattern(patternFromText(filter->getItemText(n)));
  dirlist->setDirectory(cdir);
  }

// Change the pattern
long MainWindow::onCmdFilter(FXObject*,FXSelector,void* ptr){
  dirlist->setPattern(FXFileSelector::patternFromText((FXchar*)ptr));

  FXTreeItem* curitem=dirlist->getCurrentItem();
  if (dirlist->isItemFile(curitem))
          curitem=curitem->getParent();
  dirlist->collapseTree(curitem,1);
  dirlist->expandTree(curitem,1);
  dirlist->selectItem(curitem,1);
  dirlist->makeItemVisible(curitem);
  return 1;
  }


// Clean up
MainWindow::~MainWindow(){
  delete filemenu;
  delete windowmenu;
  delete helpmenu;
  delete mdiicon;
  }

// Show up
void MainWindow::create(){
  FXMainWindow::create();
  show();
  FXTreeItem* curitem=dirlist->getCurrentItem();
  dirlist->expandTree(curitem,1);
  dirlist->makeItemVisible(curitem);
  }

// About
long MainWindow::onCmdAbout(FXObject*,FXSelector,void*){
  FXMessageBox::information(this,MBOX_OK,"About clView",
    "Assembly/Cluster layout viewer \n Written using FOX toolkit");
  return 1;
  }

// Change patterns, each pattern separated by newline
void MainWindow::setPatterns(const FXString& patterns){
  FXString pat; FXint i;
  filter->clearItems();
  for(i=0; !(pat=patterns.section('\n',i)).empty(); i++){
    filter->appendItem(pat);
    }
  if(!filter->getNumItems()) filter->appendItem("All Files (*)");
  setCurrentPattern(0);
  }


// Return list of patterns
FXString MainWindow::getPatterns() const {
  FXString pat; FXint i;
  for(i=0; i<filter->getNumItems(); i++){
    if(!pat.empty()) pat+='\n';
    pat+=filter->getItemText(i);
    }
  return pat;
  }

// Open existing file
long MainWindow::onCmdFileDlg(FXObject*,FXSelector,void*){
  FXFileDialog opendialog(this,"Open file");
  opendialog.setDirectory(FXSystem::getCurrentDirectory());
  opendialog.setSelectMode(SELECTFILE_EXISTING);
  opendialog.setPatternList(fltfiles);
  opendialog.setCurrentPattern(getCurrentPattern());
  //opendialog.setDirectory(FXFile::directory(filename));
  if(opendialog.execute()){
   FXString file=opendialog.getFilename();
   if (!file.empty()) newMDIChild(file.text());
   }
  //FXString file=FXFileDialog::getOpenFilename(this,"Open file",NULL,fltfiles,0);
  return 1;
  }


// Command from the tree list
long MainWindow::onCmdTreeList(FXObject*,FXSelector,void* ptr){
  FXTreeItem *item=(FXTreeItem*)ptr;
  FXint x,y;
  FXuint state;
  getCursorPosition(x,y,state);

  FXString file;
  if(!item || !dirlist->isItemFile(item)) return 1;
  file=dirlist->getItemPathname(item);
  //loadFile(file);         :TO DO
  MDIChild* mdichild=(MDIChild*)mdiclient->getActiveChild();
  if ((state & CONTROLMASK) && mdichild!=NULL) {
     mdichild->loadGroups(file);
     }
     else newMDIChild(file.text());
  return 1;
  }

void MainWindow::newMDIChild(const char* openfile) {
  if (openfile!=NULL) {
   FXString s;
   s.format(" Parsing file %s, please wait... \n", openfile);
   splash(s.text());
   s=openfile;
   //s.format("Cluster viewer: %s", fileNameOnly(s).text());
   s=fileNameOnly(s);
   MDIChild* mdichild=new MDIChild(mdiclient,s.text(),mdiicon,mdimenu,
     MDI_MAXIMIZED|MDI_TRACKING,10,10,600,400);
   if (!constructing)
        mdichild->create();
   mdichild->alignview->initColors(seqColors, baseColors, matchColors, 
                                        ctgQuals, grpColors[0], grpColors[1]);
   mdichild->openFile(openfile);
   splash();
   }
}

long MainWindow::onMaxRestore(FXObject* obj,FXSelector sel,void* ptr) {
 FXMDIChild* active=mdiclient->getActiveChild();
 FXString appName(getApp()->getAppName());
 if (active!=NULL && active->isMaximized()) {
       //FXString newtitle(appName);
       //newtitle.append(" - ");
       //newtitle.append(active->getTitle());
       FXString title=active->getTitle();
       //fprintf(stderr, "%s\n",title.text());
       obj->handle(this,FXSEL(SEL_COMMAND,ID_SETSTRINGVALUE),(void*)&title);
       //this->setTitle(active->getTitle());
       }
   else
       this->setTitle(appName);
//pass message to FXMainWindow
 return FXMainWindow::handle(obj,sel,ptr);
}


long MainWindow::onRandomLayout (FXObject*,FXSelector,void* ptr) {
  MDIChild* mdichild=new MDIChild(mdiclient,"Random layout",mdiicon,mdimenu,
     MDI_MAXIMIZED,10,10,600,400);
  if (!constructing)
        mdichild->create();
  mdichild->alignview->initColors(seqColors, baseColors, matchColors,
           ctgQuals, grpColors[0], grpColors[1]);
  FXString s;
  for (int i=0;i<1000;i++) {
                   s.format("TEST%04d",i);
                   mdichild->alignview->addSeq(s.text(), rand()%640+120, rand()%1024-120, rand()%2,
                         rand()%20, rand()%20);
                   }
  mdichild->alignview->buildLayout(); //only after this it is displayed correctly
  mdichild->alignview->setColorStyle(FXClView::csDensity);
  if (mdichild->alignview->HasSeqs())
      mdichild->clropt3->enable();
    else
      mdichild->clropt3->disable();
  return 1;
}

void MainWindow::splash(const char* mesg) {
 if (mesg==NULL) {
    splashmsg->hide();
    delete splashmsg;
    return;
    }
 splashmsg=new FXSplashDlg(this, mesg);
 splashmsg->create();
 splashmsg->show();
 getApp()->flush();
}

/**************************** splash/dialog */
// Construct a dialog box
FXSplashDlg::FXSplashDlg(FXWindow* ownr, const char* mesg):
  FXDialogBox(ownr,"splash",DECOR_NONE|PLACEMENT_OWNER,
                ownr->getX()+ownr->getWidth()/2-60,
                ownr->getY()+ownr->getHeight()/2-30,0,0,0,0,0,0) {
  // Contents
  contents=new FXVerticalFrame(this, FRAME_RAISED|FRAME_THICK|LAYOUT_FILL_X|LAYOUT_FILL_Y, 0,0,0,0,0,0,0,0);
  FXHorizontalFrame* fr=new FXHorizontalFrame(contents,
         FRAME_NONE|LAYOUT_FILL_X|LAYOUT_FILL_Y|LAYOUT_CENTER_Y|LAYOUT_CENTER_X,
                         0,0,0,0,20,20,20,20);
  msg=new FXLabel(fr, mesg, NULL, FRAME_NONE|LAYOUT_CENTER_Y|LAYOUT_FILL_X);
  }


// Must delete the menus
FXSplashDlg::~FXSplashDlg(){
  //anybody?
  }

