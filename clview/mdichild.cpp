/********************************************************************************
*                                                                               *
*                                                                               *
*********************************************************************************/

#include <stdio.h>
#include "mdichild.h"
#include "clv_icons.h"

/*******************************************************************************/


// ========= static members initialization: ================
short MDIChild::numInst=0; //number of instances
FXFont* MDIChild::ifont=NULL;
FXIcon* MDIChild::icnDefCS=NULL;
FXIcon* MDIChild::icnBaseCS=NULL;
FXIcon* MDIChild::icnDensCS=NULL;
//============================================================

// Map
FXDEFMAP(MDIChild) MDIChildMap[]={
  //________Message_Type____________________ID___________________Message_Handler________
  FXMAPFUNC(SEL_CHANGED, MDIChild::ID_CLVIEW, MDIChild::onViewUpdate),
  FXMAPFUNCS(SEL_CHANGED,           MDIChild::ID_CLXZOOM, 
                                    MDIChild::ID_CLYZOOM,   MDIChild::onCmdClZoom),
  FXMAPFUNC(SEL_COMMAND,            MDIChild::ID_CLNOZOOM,  MDIChild::onNoZoom),
  FXMAPFUNC(SEL_COMMAND,            MDIChild::ID_CLGRPS,    MDIChild::onCmdGrps),
  FXMAPFUNC(SEL_COMMAND,            MDIChild::ID_CL1PXZOOM,  MDIChild::onMaxZoom),
  FXMAPFUNC(SEL_CLICKED,            MDIChild::ID_CLVIEW,    MDIChild::onSeqSelect),
  FXMAPFUNC(SEL_RIGHTBUTTONPRESS,   MDIChild::ID_CLVIEW,    MDIChild::onRMouseDown),
  FXMAPFUNC(SEL_RIGHTBUTTONRELEASE, MDIChild::ID_CLVIEW,    MDIChild::onRMouseUp),
  FXMAPFUNC(SEL_MIDDLEBUTTONPRESS,   MDIChild::ID_CLVIEW,    MDIChild::onMMouseDown),
  FXMAPFUNC(SEL_MIDDLEBUTTONRELEASE, MDIChild::ID_CLVIEW,    MDIChild::onMMouseUp),
  FXMAPFUNC(SEL_MOTION,             MDIChild::ID_CLVIEW,    MDIChild::onMouseMove),
  FXMAPFUNCS(SEL_COMMAND,           MDIChild::ID_CSTYLE1,
                                    MDIChild::ID_CSTYLE3,   MDIChild::onColorOption),
  FXMAPFUNC(SEL_COMMAND,            MDIChild::ID_CONTIG,    MDIChild::onSelContig),
  FXMAPFUNC(SEL_COMMAND,            MDIChild::ID_REGION,    MDIChild::onSelRegion),
  FXMAPFUNC(SEL_COMMAND,            MDIChild::ID_HIDEGAPS,  MDIChild::onHideGaps),
  FXMAPFUNC(SEL_COMMAND,            MDIChild::ID_SELSEQ,    MDIChild::onSelSeq),
  // FXMAPFUNC(SEL_CONFIGURE,          0, MDIChild::onMaxRestore), //breaks on fox 1.7
  };


// Object implementation
FXIMPLEMENT(MDIChild,FXMDIChild,MDIChildMap,ARRAYNUMBER(MDIChildMap))


/*******************************************************************************/
// Make some windows
MDIChild::MDIChild(FXMDIClient* p,const FXString& name,FXIcon* ic,
    FXMenuPane* mn,FXuint opts,FXint x,FXint y,FXint w,FXint h)
     :FXMDIChild(p,name, ic, mn, opts, x, y, w, h) {
  FXString s;
  // Tooltip
  wasMaximized=false;
  new FXToolTip(getApp());
  layoutparser=NULL;
  groupHolder=NULL;
  zooming=false;
  panning=false;
  prevColumn=-1;
  contents=new FXVerticalFrame(this,FRAME_NONE|LAYOUT_FILL_X|LAYOUT_FILL_Y,0,0,0,0,0,0,0,0);

    vframe=new FXVerticalFrame(contents,FRAME_NONE|LAYOUT_FILL_Y|LAYOUT_FILL_X,
                 0,0,0,0,0,0,0,0);
      FXHorizontalFrame* hf=new FXHorizontalFrame(vframe,FRAME_RAISED|LAYOUT_FILL_X|LAYOUT_TOP);
         //contig and option buttons:
      //new FXLabel(hf, "Cluster:", NULL, LAYOUT_CENTER_Y);

      //selcontig=new FXComboBox(hf,32,this,ID_CONTIG,
      //                 COMBOBOX_NORMAL|FRAME_SUNKEN|FRAME_THICK|LAYOUT_CENTER_Y);
      selcontig=new FXListBox(hf,this,ID_CONTIG,
                             FRAME_SUNKEN|FRAME_THICK|LAYOUT_CENTER_Y);
      selcontig->setNumVisible(12); //TODO: should be min(count, 12)
      //selcontig->setEditable(FALSE);
      pop=new FXPopup(this);
      if (numInst<=0) {
          icnDefCS=new FXXPMIcon(getApp(), cs_default, FXRGB(192,192,192),IMAGE_ALPHACOLOR);
          icnBaseCS=new FXXPMIcon(getApp(), cs_base, FXRGB(192,192,192),IMAGE_ALPHACOLOR);
          icnDensCS=new FXXPMIcon(getApp(), cs_density, FXRGB(192,192,192),IMAGE_ALPHACOLOR);
          ifont=new FXFont(getApp(),"helvetica", 10, FXFont::Bold); //info font
          }
      numInst++;

        clropt1=new FXOption(pop,"\t\tCoverage/quality color style",
                    icnDensCS,
                     this,ID_CSTYLE2,JUSTIFY_HZ_APART);
        clropt2=new FXOption(pop,"\t\tSimple/Group color style", icnDefCS,
                     this, ID_CSTYLE1,JUSTIFY_HZ_APART);
        clropt3=new FXOption(pop,"\t\tNucleotide color style",
                    icnBaseCS,
                     this,ID_CSTYLE3,JUSTIFY_HZ_APART);
      new FXOptionMenu(hf, pop,
           LAYOUT_TOP|FRAME_RAISED|FRAME_THICK|JUSTIFY_HZ_APART);
      new FXButton(hf, "G\tLoad groups\tShow groups loaded from another layout file",
              NULL,this, ID_CLGRPS,
              LAYOUT_TOP|FRAME_RAISED|FRAME_THICK|JUSTIFY_HZ_APART);
      cbgaps=new FXCheckButton(hf,"hide contig gaps",this, MDIChild::ID_HIDEGAPS);
      new FXLabel(hf, "    Region:", NULL, LAYOUT_CENTER_Y);
      selregion=new FXTextField(hf, 32, this, ID_REGION, FRAME_SUNKEN|FRAME_THICK|LAYOUT_CENTER_Y);

      hframe=new FXHorizontalFrame(vframe,FRAME_NONE|LAYOUT_FILL_X|LAYOUT_FILL_Y|ICON_BEFORE_TEXT,
                 0,0,0,0,0,0,0,0);
      vframe= new FXVerticalFrame(hframe,FRAME_NONE|LAYOUT_FILL_Y|LAYOUT_FILL_X,
                 0,0,0,0,0,0,0,0);
           clframe=new FXHorizontalFrame(vframe,FRAME_SUNKEN|LAYOUT_FILL_X|LAYOUT_FILL_Y|LAYOUT_BOTTOM,
                     0,0,0,0,0,0,0,0);
               alignview=new FXClView(clframe, this,MDIChild::ID_CLVIEW, LAYOUT_FILL_X|LAYOUT_FILL_Y);

            hf = new FXHorizontalFrame(vframe,FRAME_SUNKEN|LAYOUT_FILL_X,
                     0,0,0,0,0,0,0,0);
             //=== seq_info box:
            FXVerticalFrame* vf = new FXVerticalFrame(hf, FRAME_NONE|LAYOUT_CENTER_Y|LAYOUT_FILL_X);
             FXHorizontalFrame* hhf=new FXHorizontalFrame(vf,
                       FRAME_NONE|LAYOUT_CENTER_Y|LAYOUT_FILL_X, 0,0,0,0,0,0,0,0);
             selseq=new FXComboBox(hhf,28,this,ID_SELSEQ,
                       COMBOBOX_NORMAL|FRAME_SUNKEN|FRAME_THICK|LAYOUT_CENTER_Y);
             selseq->setNumVisible(12);
             selseq->setEditable(FALSE);
             seqData = new FXLabel(hhf, "                      ", NULL,
                   FRAME_SUNKEN|JUSTIFY_LEFT|LAYOUT_FILL_X);
                   //seqComment=new FXTextField(vf,120,NULL,0,JUSTIFY_LEFT|FRAME_SUNKEN|FRAME_THICK|LAYOUT_FILL_X);
                seqData->setFont(ifont);
                seqAsmInfo=new FXLabel(vf, "                                      ", NULL,
                   FRAME_SUNKEN|JUSTIFY_LEFT|LAYOUT_FILL_X);
                seqAsmInfo->setFont(ifont);
             //===slider box:
             vf = new FXVerticalFrame(hf, FRAME_RAISED);
               sliderZX=new FXSlider(vf,this,MDIChild::ID_CLXZOOM,
                  LAYOUT_FIX_WIDTH|LAYOUT_FIX_HEIGHT|FRAME_RAISED,
                   0,0,160,22,0,0,0,0);

               hf=new FXHorizontalFrame(vf,FRAME_NONE|LAYOUT_FILL_X,0,0,0,0,0,0,0,0);
                 nozoomBtn=new FXButton(hf, "&1x1",NULL, this, MDIChild::ID_CLNOZOOM,
                    FRAME_RAISED|FRAME_THICK|LAYOUT_CENTER_Y);
                 new FXButton(hf, "1p&x",NULL, this, MDIChild::ID_CL1PXZOOM,
                    FRAME_RAISED|FRAME_THICK|LAYOUT_CENTER_X|LAYOUT_CENTER_Y);
                 zxLabel=new FXLabel(hf, "1.000", NULL,
                  FRAME_SUNKEN|JUSTIFY_RIGHT|LAYOUT_FIX_WIDTH|LAYOUT_FIX_HEIGHT|LAYOUT_RIGHT|LAYOUT_CENTER_Y,
                  0,0,6*6,14);

             s.format("%02.3f",alignview->getZoomX());
             zxLabel->setText(s);
             columnInfo=new FXLabel(vframe, "                                ",
                                NULL, JUSTIFY_LEFT|FRAME_SUNKEN|LAYOUT_FILL_X);
 //=== vframe filled ===

    toolframe=new FXVerticalFrame(hframe,FRAME_RAISED|LAYOUT_FILL_Y);
    //new FXLabel(toolframe, "Y Scale:", NULL, LAYOUT_CENTER_X);
    sliderZY=new FXSlider(toolframe,this,MDIChild::ID_CLYZOOM,
             SLIDER_VERTICAL|LAYOUT_FIX_WIDTH|LAYOUT_FIX_HEIGHT|LAYOUT_CENTER_X|FRAME_RAISED,
             0,0,22,140,0,0,0,0);

    sliderZY->setRange(10,100); //
    sliderZY->setValue((FXint)(alignview->getZoomY()*100));

    sliderZX->setRange(10,1000);//
    sliderZX->setValue((FXint)(alignview->getZoomX()*1000));
    s.format("%02.2f",alignview->getZoomY());
    zyLabel=new FXLabel(toolframe, s, NULL,
         FRAME_SUNKEN|LAYOUT_CENTER_X|LAYOUT_FIX_WIDTH|LAYOUT_FIX_HEIGHT,
         0,0,6*6,14);
         
  alignview->setFocus();
  }

void MDIChild::showSeqInfo(ClSeq* seq) {
 FXString s;

 if (seq!=NULL && seq==alignview->selSeq) { //selection
     s.format("%s (%d nt)", seq->name.text(), seq->len);
     selseq->setText(s);
     if (seq->ins>0 || seq->dels>0)
        s.format(" %s   range: %d to %d | clipL:%d clipR:%d | gaps: %d  dels: %d",
            seq->reversed?"-":"+",
            seq->clipL+1, seq->len-seq->clipR,
            seq->clipL, seq->clipR, 
            seq->ins, seq->dels);
       else 
        s.format(" %s   range: %d to %d | clipL:%d clipR:%d ",
            seq->reversed?"-":"+",
            seq->clipL+1, seq->len-seq->clipR,
            seq->clipL, seq->clipR);
    seqData->setText(s);
    //seqComment->setText(seq->comment);
    if (seq->group>=0)
       s.format("Layout location: %d to %d  [group %d]",
        seq->cl_xpos+seq->clipL, seq->cl_xpos+seq->len-seq->clipR-1, seq->group);
     else
       s.format("Layout location: %d to %d",
        seq->cl_xpos+seq->clipL, seq->cl_xpos+seq->len-seq->clipR-1);
    seqAsmInfo->setText(s);
    }
  else { //deselect
    selseq->setText("  ");
    seqData->setText("              ");
    seqAsmInfo->setText("             ");
    }
 int cidx=alignview->getSelCol();
 if (prevColumn==cidx) return;
 prevColumn=cidx;
 FXString cinfo;
 if (cidx>=0) {
    ColumnData* c=alignview->getColumnData(cidx);
    if (c->letter!=0) {
       cinfo.format("Column %d  ['%c']: %d layers, %d mismatches ",
          cidx+alignview->XLeft, c->letter, c->thickness, c->mismatches);
       if (c->mismatches>1 && c->thickness>3) {
         cinfo+=" [ ";
         for (int i=0;i<3;i++) {
           if (c->ntdata[i].count>0) {
            s.format(" %d%c ", c->ntdata[i].count, c->ntdata[i].letter);
            cinfo+=s;
           }
          }
         cinfo+=" ] "; 
         }
        }
      else
        cinfo.format("Column %d (%d layers)", cidx+alignview->XLeft, c->thickness);
    columnInfo->setText(cinfo);
    }
  else {
    columnInfo->setText(" ");
    }
}

long MDIChild::onSeqSelect(FXObject*,FXSelector sel,void* ptr){
 showSeqInfo((ClSeq*)ptr);
 return 1;
}

long MDIChild::onCmdClZoom(FXObject*,FXSelector sel,void*){
 FXuint sid=FXSELID(sel);
 FXString s;
 switch(sid) {
    case ID_CLXZOOM:
      alignview->ZoomX(((double)sliderZX->getValue())/1000);
      s.format("%2.3f",alignview->getZoomX());
      zxLabel->setText(s);
      break;
    case ID_CLYZOOM:
      alignview->ZoomY(((double)sliderZY->getValue())/100);
      s.format("%2.2f",alignview->getZoomY());
      zyLabel->setText(s);
      break;
    }
 return 1;
}

long MDIChild::onColorOption(FXObject*,FXSelector sel,void*){
 FXuint sid=FXSELID(sel);
 switch(sid) {
    case ID_CSTYLE1:
      alignview->setColorStyle(FXClView::csDefault);
      break;
    case ID_CSTYLE2:
      alignview->setColorStyle(FXClView::csDensity);
      break;
    case ID_CSTYLE3:
      alignview->setColorStyle(FXClView::csBaseColor);
      break;
    }
 this->ungrab(); //? bug?
 return 1;
}

// Show up
void MDIChild::create(){
  FXMDIChild::create();
  show();
  }


long MDIChild::onRMouseDown(FXObject*,FXSelector,void* ptr) {
  FXEvent *ev = (FXEvent*)ptr;
  //prevCursor=getApp();
  alignview->setDragCursor(getApp()->getDefaultCursor(DEF_MOVE_CURSOR));
  alignview->grab();
  // While the right-button mouse is down, do zooming
  zooming=true;
  startZX=alignview->getZoomX();
  startZY=alignview->getZoomY();
  zoomPt.x=ev->win_x;
  zoomPt.y=ev->win_y;
  return 1;
  }

long MDIChild::onMMouseDown(FXObject*,FXSelector,void* ptr) {
  //FXEvent *ev = (FXEvent*)ptr;
  //prevCursor=getApp();
  alignview->setDragCursor(getApp()->getDefaultCursor(DEF_MOVE_CURSOR));
  alignview->grab();
  // While the right-button mouse is down, do zooming
  panning=true;
  //startZX=alignview->getZoomX();
  //startZY=alignview->getZoomY();
  //zoomPt.x=ev->win_x;
  //zoomPt.y=ev->win_y;
  return 1;
  }

long MDIChild::onViewUpdate(FXObject*,FXSelector,void*){
 //on hscroll or zooming, update the range display
	//TODO: this should be triggered by FXClView somehow, when its content
	//is scrolled horizontally or zoomed in/out
	FXString region;
	FXint xstart=0, xend=0;
	alignview->getVisibleRegion(xstart, xend);
	region.format("%d - %d", xstart, xend);
  selregion->setText(region);
  return 1;
}

long MDIChild::onMouseMove(FXObject*, FXSelector, void* ptr){
  FXEvent *ev=(FXEvent*)ptr;
  if (zooming) {
    //compute distance:
    //int dx=(ev->win_x-zoomPt.x)*2;
    //int dy=(zoomPt.y-ev->win_y)/3;
    double newzx = startZX+((double)(ev->win_x-zoomPt.x)*2)/1000;
    double newzy = startZY+(((double)(zoomPt.y-ev->win_y))/3)/100;
    int minzx, maxzx;
    int minzy, maxzy;
    sliderZX->getRange(minzx, maxzx);
    sliderZY->getRange(minzy, maxzy);
    if (newzx*1000<minzx) newzx=(double)minzx/1000;
       else
         if (newzx*1000>maxzx) newzx=(double)maxzx/1000;
    if (newzy*100<minzy) newzy=(double)minzy/100;
       else
        if (newzy*100>maxzy) newzy=(double)maxzy/100;
    if (newzx!=alignview->getZoomX() ||
            newzy!=alignview->getZoomY()) {
      alignview->Zoom(newzx,newzy,
                        zoomPt.x, zoomPt.y);
      FXString s;
      sliderZX->setValue((FXint)(newzx*1000));
      sliderZX->forceRefresh();
      s.format("%2.3f",alignview->getZoomX());
      zxLabel->setText(s);
      sliderZY->setValue((FXint)(newzy*100));
      sliderZY->forceRefresh();
      s.format("%2.2f",alignview->getZoomY());
      zyLabel->setText(s);
      }
    }
  if (panning) {
    alignview->scrollBy((ev->win_x-ev->last_x)*2, (ev->win_y-ev->last_y)*2);
    }
  return 1;
  }


long MDIChild::onNoZoom(FXObject*,FXSelector sel,void*) {
  sliderZX->setValue(1000);
  zxLabel->setText("1.000");
  sliderZY->setValue(100);
  zyLabel->setText("1.00");
  alignview->Zoom(1,1);
  return 1;
}

long MDIChild::onMaxZoom(FXObject*,FXSelector sel,void*) {
  FXString s;
  double newzx, newzy;
  alignview->set1pxZoom();
  newzx=alignview->getZoomX();
  newzy=alignview->getZoomY();
  sliderZX->setValue((FXint)(newzx*1000));
  sliderZX->forceRefresh();
  s.format("%2.3f",newzx);
  zxLabel->setText(s);
  sliderZY->setValue((FXint)(newzy*100));
  sliderZY->forceRefresh();
  s.format("%2.2f",newzy);
  zyLabel->setText(s);
  return 1;
}

long MDIChild::onCmdGrps(FXObject*,FXSelector sel,void*) {
FXString ofile=FXFileDialog::getOpenFilename(this,"Open file", //NULL,
    "Contig files (*.ace,*.lyt)\nAny file (*)");
if (!ofile.empty()) loadGroups(ofile);
return 1;
}



long MDIChild::onRMouseUp(FXObject*,FXSelector,void* ptr){
  //FXEvent *ev=(FXEvent*) ptr;
  alignview->ungrab();
  zooming=false;
  if (!panning)
   alignview->setDragCursor(getApp()->getDefaultCursor(DEF_ARROW_CURSOR));
  return 1;
  }

long MDIChild::onMMouseUp(FXObject*,FXSelector,void* ptr){
  alignview->ungrab();
  panning=false;
  if (!zooming) 
    alignview->setDragCursor(getApp()->getDefaultCursor(DEF_ARROW_CURSOR));
  return 1;
  }


bool MDIChild::openFile(const char* filename) {
 //try to recognize the file type and call the appropriate loader
 FILE* f;
 if ((f=fopen(filename, "r"))==NULL) {
   FXMessageBox::error(this,MBOX_OK,"File error",
    "Cannot open file '%s'", filename);
   }
 fclose(f);
 FXString fname(filename);
 AceParser ap(fname.text());
 if (ap.open()) {
     delete layoutparser;
     layoutparser=NULL;
     layoutparser=new AceParser(fname.text());
     return loadLayoutFile(fname);
     }
  else { //assumes it's our layout format file
    delete layoutparser;
    layoutparser=NULL;
    layoutparser=new LayoutParser(fname.text());
    }
 return loadLayoutFile(fname);
}


bool MDIChild::loadLayoutFile(FXString& lytfile) {
 if (!layoutparser->open()) {
    FXMessageBox::error(this,MBOX_OK,"File error",
      "Cannot open file '%s'. Not a valid file?", lytfile.text());
    delete layoutparser;
    layoutparser=NULL;
    return false;
    }
 alignview->filetype=layoutparser->getFileType();
 isAce=(alignview->filetype=='A');
 if (!layoutparser->parseContigs() || layoutparser->getNumContigs()==0) {
 //if (!layoutparser->parse() || layoutparser->getNumContigs()==0) {
	FXMessageBox::error(this,MBOX_OK,"File error",
	 "Error encountered parsing file '%s'. Not a valid file?", lytfile.text());
	delete layoutparser;
	layoutparser=NULL;
	return false;
	}
 file=lytfile;
 //file validated. Load align data:
 //show all contig labels:
 layoutparser->contigsByNumSeqs();
 FXString name;
 for (int i=0; i<layoutparser->getNumContigs(); i++) {
   LytCtgData* ctg=layoutparser->getContig(i);
   name.format("%s (%d seqs, %d nts)", ctg->name, ctg->numseqs,
          ctg->len);
   selcontig->appendItem(name);
   }
 //show first contig:

 selcontig->setCurrentItem(0);
 selContig(0, cbgaps->getCheck());

 if (alignview->HasSeqs()) {
      clropt3->enable();
      cbgaps->show();
      }
    else {
      clropt3->disable();
      cbgaps->hide();
      }
 alignview->setColorStyle(FXClView::csDensity);
 return true;
}

void MDIChild::loadGroups(FXString& lytfile) {
 //assumes sequences already loaded here
 if (alignview->seqlist->Count()==0) {
   FXMessageBox::error(this,MBOX_OK,"Action error",
      "You need a valid layout file loaded before setting groups.");
   return;
   }
 LayoutParser* lp = new AceParser(lytfile.text());
 if (!lp->open()) { //try ACE first
     delete lp;
     lp=new LayoutParser(lytfile.text());
     }
    else lp->close();
 if (!lp->open() || 
          !lp->parse()) {
     FXMessageBox::error(this,MBOX_OK,"File error",
      "Error encountered opening/parsing layout file '%s'. Not a layout file?", lytfile.text());
     delete lp;
     return;
     }
 delete groupHolder;
 //groupHolder=NULL;
 groupHolder=lp;
 assignGroups();
 }
 
void MDIChild::assignGroups() {
 if (groupHolder==NULL) return;
 int numGroups=0;
 //we must count only contigs having at least one
 //sequence in common with the currently displayed layout
 //for each contig, test all component sequences
 //groupHolder->parseContigs();
 int nc=groupHolder->getNumContigs();
 for (int i=0; i<nc; i++) {
   LytCtgData* ctg=groupHolder->getContig(i);
   bool assigned=false;
   //for each sequence, see if it's in the viewer's current layout
   for (int j=0;j<ctg->numseqs;j++) {
     FXString s(ctg->seqs[j]->name);
     if (alignview->setSeqGrp(s, numGroups)) assigned=true;
     }
   if (assigned) numGroups++; //advance to next group  
   }
 if (numGroups>0) {
     //alignview->initGrpColors(numGroups);
     alignview->update();
     }
   else
     FXMessageBox::information(this,MBOX_OK,"Information",
      "No groups were found for sequences in the selected contig");
}

void MDIChild::selContig(int ctgno, bool hideGaps) {
 if (ctgno<0 || ctgno>=layoutparser->getNumContigs()) return;
 alignview->Clear();
 FXApp* fapp=getApp();
 fapp->beginWaitCursor();
 fapp->forceRefresh();

 layoutparser->loadContig(ctgno);
 LytCtgData* ctg=layoutparser->getContig(ctgno);
 char *ctgseq=layoutparser->getContigSeq(ctg);
 //now ctgseq contains '*'
 //for each position, record how many '*' chars we have on the left
 int* ctgrm=NULL;
 bool removeGaps=hideGaps;
 int ctglen=ctg->len;
 if (removeGaps && ctgseq!=NULL) {
    int len=strlen(ctgseq);
    GMALLOC(ctgrm, (len+1)*sizeof(int));
    int numstars=0;
    ctgrm[0]=0;
    for (int i=1;i<len;i++) {
      numstars += (ctgseq[i-1]=='*');
      ctgrm[i]=numstars;
      }
    ctgrm[len]=ctgrm[len-1];
    //now remove all the stars from the contig sequence
	if (numstars>0) {
       char* cleanseq=NULL;
       GMALLOC(cleanseq, len-numstars+1);
       ctglen-=numstars;
       int j=0;
       for (int i=0;i<=len;i++)
         if (ctgseq[i]!='*') {
               cleanseq[j] = ctgseq[i];
               j++;
               }
       GFREE(ctgseq);
       ctgseq=cleanseq;
       removeGaps=true;
	   }
	  else removeGaps=false;
    }
 alignview->addContig(ctg->name, ctglen, ctgseq, ctg->offs);
 GFREE(ctgseq);
 int numseqs=ctg->numseqs;
 for (int i=0;i<numseqs;i++) { //reload each read sequence
   LytSeqInfo* seq=ctg->seqs[i];

   char* s=layoutparser->getSeq(seq);
   //   GMessage("seq %d : %s \n", i, seq->name);

   int slen=0;
   if (isAceFile() && (s==NULL || (slen=strlen(s))!=(unsigned int)seq->seglen())) {
          FXMessageBox::error(this,MBOX_OK,"File error",
           "Error reading sequence '%s' (declared %d, loaded %d). Invalid ACE file?",
                 seq->name, seq->length(), slen);
          delete layoutparser;
          layoutparser=NULL;
          GFREE(s);
          return;
          }
   if (removeGaps && s!=NULL) {
	 //shift all coordinates and lengths accordingly, remove gaps
	 int clpL=seq->left-1;
	 int clpR=seq->length()-seq->right;
	 int old_start=seq->offs+clpL;
	 int old_end=seq->offs+seq->length()-clpR-1;
	 int new_start=old_start-ctgrm[old_start-1];
	 int new_end=old_end-ctgrm[old_end-1];
	 //fix interseg coordinates too, when introns are given
	 // otherwise the "remove gaps" is going to mess up the MSA !!
         
	 for (int g = 0; g < seq->numisegs; g++) {
	   int prevcoord=seq->intersegs[g].segEnd;
	   seq->intersegs[g].segEnd-=ctgrm[prevcoord-1];
	   prevcoord=seq->intersegs[g].nextStart;
	   seq->intersegs[g].nextStart-=ctgrm[prevcoord-1];
	   } 
	 //now create the cleaned copy of the sequence
         char* clnseq=NULL;
         int newlen=new_end-new_start+1+clpL+clpR;
         GMALLOC(clnseq, sizeof(char)*(newlen+1));
         clnseq[newlen]='\0';
         if (clpL>0) strncpy(clnseq, s, clpL);
         int j=clpL;
         int numdels=0;
         int numins=0;
         for (int b=clpL; b<seq->length()-clpR; b++) {
             int ctgpos=seq->offs+b;
             //delete only if there was also a '*' in the initial contig seq
             if (ctgrm[ctgpos]-ctgrm[ctgpos-1]==0) {
                   clnseq[j] = s[b];
                   if (s[b]=='*' || s[b]=='-') {
                      clnseq[j]='-';
                      numins++;
                      }
                   j++;
                   }
                 else {
                  if (s[b]!='*' && s[b]!='-')
                         numdels++;
                  }
             }
         clnseq[j]='\0';
         if (clpR>0) strncat(clnseq, s+(seq->length()-clpR), clpR);
         GFREE(s);
         s=clnseq;
         alignview->addSeq(seq->name, newlen, new_start-clpL, seq->reversed,
             clpL,clpR,s, numins, numdels);
         //alignview->addSeq(seq,s);
     } //gaps removal
   else { // raw sequence display (with all gaps in place)  
     /*alignview->addSeq(seq->name, seq->length(), seq->offs, seq->reversed,
         seq->left-1,seq->length()-seq->right,s);*/
       alignview->addSeq(seq,s);  
     }
   GFREE(s);
   }
 alignview->buildLayout();
 assignGroups();
 selseq->clearItems();
 showSeqInfo(NULL);
 for (int i=0;i<alignview->seqlist->Count();i++) {
   ClSeq* seq=alignview->seqlist->Get(i);
   if (seq==alignview->contig) continue;
   FXString name;
   name.format ("%s (%d nt)", seq->name.text(), seq->len);
   selseq->appendItem(name);
   }
 fapp->endWaitCursor();
 if (alignview->HasSeqs()) {
      clropt3->enable();
      cbgaps->show();
      }
    else {
      clropt3->disable();
      cbgaps->hide();
      }
}

long MDIChild::onSelContig(FXObject*,FXSelector,void* ptr) {
  selContig(selcontig->getCurrentItem(), cbgaps->getCheck());
  return 1;
  }

long MDIChild::onSelRegion(FXObject*,FXSelector,void* ptr) {
  //selContig(selcontig->getCurrentItem(), cbgaps->getCheck());
  return 1;
  }


long MDIChild::onHideGaps(FXObject*,FXSelector,void* ptr) {
  selContig(selcontig->getCurrentItem(), cbgaps->getCheck());
  return 1;
  }


long MDIChild::onSelSeq(FXObject*,FXSelector,void* ptr) {
  int seqidx=selseq->getCurrentItem();
  if (seqidx<0) return 1;
  ClSeq* seq=alignview->seqlist->Get(seqidx);
  alignview->selectSeq(seq, true);
  showSeqInfo(seq);
  return 1;
  }
/* fox 1.7 - seems to no longer work with this:
long MDIChild::onMaxRestore(FXObject* sender,FXSelector sel,void* ptr) {  
 FXbool maximized=isMaximized();
 if (wasMaximized!=maximized) {
     FXWindow* mainwin=getShell();
     mainwin->handle(mainwin, FXSEL(SEL_COMMAND, ID_MAXRESTORE),ptr);
     wasMaximized=maximized;
     }
 return FXMDIChild::handle(sender,sel,ptr); 
}
*/

MDIChild::~MDIChild() {
 delete pop;
 delete layoutparser;
 delete groupHolder;
 //== static members:
 numInst--;
 if (numInst==0) {
   delete ifont;
   delete icnDefCS;
   delete icnBaseCS;
   delete icnDensCS;
   }
 }


