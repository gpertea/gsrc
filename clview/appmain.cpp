#include "fx.h"
#include "mainwin.h"
//---------------------------------------------------------------------------
// Start the whole thing
int main(int argc,char *argv[]){

  // Make application
  FXApp application("Cluster viewer","Fox app");

  // Open display; this reads registry
  application.init(argc,argv); //removes the fox-specific parameters, if any
  //FXFont* appfont=new FXFont(&application,"helvetica",9,FONTWEIGHT_NORMAL);
  FXFont* appfont=new FXFont(&application,"helvetica",9, FXFont::Normal);
  application.setNormalFont(appfont);

  // Build a window
  //new MainWindow(&application);
  // Make window
  MainWindow* mainwindow=new MainWindow(&application);
  application.setMenuPause(30);
  // Create app
  application.create();

  mainwindow->show();
  for (int i=1;i<argc;i++) {
    mainwindow->newMDIChild(argv[i]);
    }

  // Run
  int r=application.run();
  delete appfont;
  return r;
}
//---------------------------------------------------------------------------
