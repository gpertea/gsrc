#include "fxkeys.h"
#include "fx.h"
#include <math.h>
#include "FXClView.h"
#include "clrutils.h"
#include "stdlib.h"
/*******************************************************************************/

// Map
FXDEFMAP(FXClView) FXClViewMap[]={
  FXMAPFUNC(SEL_PAINT,0,FXClView::onPaint),
  //FXMAPFUNC(SEL_CONFIGURE,0,FXClView::onConfigure),
  FXMAPFUNC(SEL_LEFTBUTTONPRESS,0,FXClView::onLeftBtnPress),
  FXMAPFUNC(SEL_LEFTBUTTONRELEASE,0,FXClView::onLeftBtnRelease),
  FXMAPFUNC(SEL_CLICKED,0,FXClView::onClicked),
  FXMAPFUNC(SEL_MOTION,0,FXClView::onMotion),
  FXMAPFUNC(SEL_COMMAND,0,FXClView::onCommand),
  FXMAPFUNC(SEL_QUERY_TIP,0,FXClView::onQueryTip),
  FXMAPFUNC(SEL_QUERY_HELP,0,FXClView::onQueryHelp),
  FXMAPFUNC(SEL_TIMEOUT, FXList::ID_TIPTIMER, FXClView::onTipTimer)
  };

// Object implementation
FXIMPLEMENT(FXClView,FXScrollArea,FXClViewMap,ARRAYNUMBER(FXClViewMap))

/*******************************************************************************/
#define CL_BORDER 6  //border size (around the whole cluster layout)

ClSeq::ClSeq(FXClView* v, LytSeqInfo* lytseq, const char* seq) {
    name=lytseq->name;
    len=lytseq->length(); //except when segmented
    sequence=Gstrdup(seq);
    view=v;
    flag=0;
    group = -1;
    cdsofs = 0;
    clipL=lytseq->left-1;
    clipR=lytseq->length()-lytseq->right;
    cl_xpos=lytseq->offs;
    cl_ypos=-1;
    ins=0;
    dels=0;
    numsegs=0;
    segs=NULL;
    isET= (name.find("np|")>=0 || name.find("et|")>=0 ||
            name.find("egad|")>=0)?1:0;
    reversed=(lytseq->reversed==1);
    //--- now, add the segment definitions, if any:
    if (lytseq->numisegs>0) {
     numsegs=lytseq->numisegs+1;
     GMALLOC(segs,numsegs*sizeof(ClSeqSegment));
	 int segStart=lytseq->offs;
	 char LSpl=0;
	 int segSeqOfs=0;
	 int segLClip=clipL; //now these "segment" coords DO include clipped ends
			 //while those coords stored in LytSeqInfo::segments EXCLUDE clipping
	 int remainingLen=len; //sequence "length" left (excluding introns!)
	 int nonclipStart=lytseq->offs+clipL;
	 for (int i=0;i<lytseq->numisegs;i++) {
	   segs[i].xL = segStart;
	   segs[i].xR = lytseq->intersegs[i].segEnd+lytseq->intersegs[i].segRClip;
	   segs[i].clipL = segLClip;
	   segs[i].clipR = lytseq->intersegs[i].segRClip;
	   segs[i].splL=LSpl;
	   segs[i].splR=lytseq->intersegs[i].segRSplice;
	   segs[i].seqofs=segSeqOfs;
	   segSeqOfs=lytseq->intersegs[i].nextSegSeq;
	   segStart=lytseq->intersegs[i].nextStart-lytseq->intersegs[i].nextLClip;
				 //+1; --not needed because interseg coords are in fact
				 //      given as seg-end coordinates!
	   segLClip=lytseq->intersegs[i].nextLClip;
	   LSpl=lytseq->intersegs[i].nextLSplice;
	   remainingLen-=(segs[i].xR-segs[i].xL+1);
	   nonclipStart=lytseq->intersegs[i].nextStart;
	   }//for each intersegment interval
	 //--add the final segment:
	 segs[numsegs-1].xL = segStart;
	 segs[numsegs-1].clipL = segLClip;
	 segs[numsegs-1].splL=LSpl;
	 if (v->filetype=='L') {
		segs[numsegs-1].xR = segStart+remainingLen-2;
		len=segStart+remainingLen-lytseq->offs-1;
		}
	  else {
		segs[numsegs-1].xR = lytseq->offs+len-1;
		}
	 segs[numsegs-1].clipR = clipR;
	 segs[numsegs-1].splR=0;
	 segs[numsegs-1].seqofs=segSeqOfs;
	 }//segmented sequence
	else { //not segmented: segment = whole sequence
	 addSegment(lytseq->offs, lytseq->offs+lytseq->length()-1, clipL, clipR);
	}
	}//FXClSeq constructor from LytSeqInfo

// Serialization
FXClView::FXClView():rows(false,true,false), columns(false,true,false) {
  vscroll=false;
  hscroll=false;
  flags|=FLAG_ENABLED;
  font=(FXFont*)-1;
  seqfont=(FXFont*)-1;
  gridH=12;
  showContig=false;
  seqfntH=12;
  seqfntW=8;
  seqH=seqfntH;
  seqW=seqfntW;
  grpColors[0]=NULL;
  grpColors[1]=NULL;
  XRight=0;
  XLeft=0;
  selSeq=NULL;
  content_h=0;
  content_w=0;
  selCol=-1;
  scale.sx=1;
  scale.sy=1;
  vSpace=seqH/4;
  backbuf=NULL;
  seqlist=NULL;
  filetype=0;
  //seqaligns=NULL;
  colorStyle=csDefault;
  canShowSeq=false;
  hasSeqs=false;
  maxThickness=0;
  contig=NULL;
  }

// Public constructor
FXClView::FXClView(FXComposite *p,FXObject* tgt,FXSelector sel,FXuint opts,FXint x,FXint y,FXint w,FXint h):
  FXScrollArea(p,opts,x,y,w,h), rows(false,true,false), columns(false,true,false){

  seqlist=new GList<ClSeq>(false,true,false);
  vscroll=false;
  hscroll=false;
  //seqaligns=new GList<ClAlign>(false,true,false);
  flags|=FLAG_ENABLED;
  font=new FXFont(getApp(), "helvetica", 8, FXFont::Normal);
  seqfont=new FXFont(getApp(),"courier", 10, FXFont::Normal);
  backbuf=NULL;
  gridH=12;
  grpColors[0]=NULL;
  grpColors[1]=NULL;
  showContig=false;
  seqfntH=12; //seqfont->getFontHeight();
  seqfntW=8; //seqfont->getFontWidth();
  seqH=seqfntH;
  seqW=seqfntW;
  XRight=0;
  XLeft=0;
  selSeq=NULL;
  selCol=-1;
  vSpace=seqH/4;
  scale.sx=1;
  scale.sy=1;
  target=tgt;
  message=sel;
  filetype=0;
  colorStyle=csDefault;
  //colorStyle=csDensity;
  canShowSeq=false;
  hasSeqs=false;
  contig=NULL;
  maxThickness=0;
  //initSeqColors();
  }

// Enter window
long FXClView::onEnter(FXObject* sender,FXSelector sel,void* ptr){
  FXScrollArea::onEnter(sender,sel,ptr);
  getApp()->addTimeout(this,ID_TIPTIMER,getApp()->getMenuPause());
  lastMousePos.x=-1;
  return 1;
  }
// Leave window
long FXClView::onLeave(FXObject* sender,FXSelector sel,void* ptr){
  FXScrollArea::onLeave(sender,sel,ptr);
  getApp()->removeTimeout(this,ID_TIPTIMER);
  lastMousePos.x=-1;
  return 1;
  }

// We timed out, i.e. the user didn't move for a while
long FXClView::onTipTimer(FXObject*,FXSelector,void*){
  FXTRACE((250,"%s::onTipTimer %p\n",getClassName(),this));
  flags|=FLAG_TIP;
  return 1;
  }

// We were asked about tip text
long FXClView::onQueryTip(FXObject* sender,FXSelector,void*){
  if((flags&FLAG_TIP)){   // No tip when autoselect!
    int col;
    ClSeq* seq=getSeqAt(lastMousePos.x, lastMousePos.y, col); //get item under cursor - this might take a while!
    if(seq!=NULL){
      FXString string=seq->name;
      sender->handle(this,MKUINT(ID_SETSTRINGVALUE,SEL_COMMAND),(void*)&string);
      return 1;
      }
    }
  return 0;
  }

// We were asked about status text
long FXClView::onQueryHelp(FXObject* sender,FXSelector,void*){
  if((flags&FLAG_HELP)) {
    FXTRACE((250,"%s::onQueryHelp %p\n",getClassName(),this));
    int col;
    ClSeq* seq=getSeqAt(lastMousePos.x, lastMousePos.y, col);
    //also compute the approximate layout position
    int xchpos = col+XLeft;
    FXString help;
    help.format("[ offset: %d ]   ",xchpos);
    if (seq!=NULL) {
      help+=seq->name;
      if (!seq->comment.empty())
          help+=" | "+seq->comment;
      }
    sender->handle(this,MKUINT(ID_SETSTRINGVALUE,SEL_COMMAND),(void*)&help);
    return 1;
    }
  return 0;
  }


// Command message
long FXClView::onCommand(FXObject*,FXSelector,void* ptr){
  return target && target->handle(this,MKUINT(message,SEL_COMMAND),ptr);
  }

// Mouse moved
long FXClView::onMotion(FXObject*,FXSelector,void* ptr){
  FXEvent* event=(FXEvent*)ptr;
  FXuint flg=flags;
  // Kill the tip
  flags&=~FLAG_TIP;
  // Kill the tip timer
  getApp()->removeTimeout(this,ID_TIPTIMER);
  lastMousePos.x=event->win_x;
  lastMousePos.y=event->win_y;
  // Reset tip timer if nothing's going on
  getApp()->addTimeout(this,ID_TIPTIMER,getApp()->getMenuPause());
  if (target)
    target->handle(this,MKUINT(message,SEL_MOTION),ptr);
  // Get item we're over
  //cursor=getItemAt(event->win_x,event->win_y);
  // Force GUI update only when needed
  //return (cursor!=oldcursor)||(flg&FLAG_TIP);
  return (flg&FLAG_TIP);
}


void FXClView::drawColumn(int colno) {
 if (colno<0) return;
 update(iround(seqW*colno)+CL_BORDER+getContentX(), 0, iround(seqW)+1, getVisibleHeight());
 }

// Pressed a button
long FXClView::onLeftBtnPress(FXObject*,FXSelector,void* ptr){
  FXEvent* event=(FXEvent*)ptr;
  flags&=~FLAG_TIP;
  handle(this,MKUINT(0,SEL_FOCUS_SELF),ptr);
  if (isEnabled()) {
    grab();
    flags&=~FLAG_UPDATE;
    // Locate seq under cursor
    btnSeq=getSeqAt(event->win_x,event->win_y, btnCol);
    if (btnCol>=0) {
      if (selCol==btnCol) { //deselect selCol
            int c=selCol;
            selCol=-1;
            drawColumn(c);
            }
          else {
          //select btnCol, deselect selCol
            int c=selCol; //previous selection
            selCol=btnCol; //new selection
            if (c>=0) drawColumn(c);
            drawColumn(selCol);
            }
      }//column selection
    // First change callback
    if (target && target->handle(this,MKUINT(message,SEL_LEFTBUTTONPRESS),ptr))
        return 1;
    // No seq hit
    flags|=FLAG_PRESSED;
    // Change current selection - not yet
    //selSeq=btnSeq;
    return 1;
    }
  return 0;
  }

//select a sequence and optionally bring it into view
void FXClView::selectSeq(ClSeq* seq, bool showIt) {
if (seq==NULL) {
    //deselect sequence
   ClSeq* p=selSeq;
   selSeq=NULL;
   if (p!=NULL) drawSeq(p);
   }
 else {//real selection
   //select btnSeq, deselect selSeq
   ClSeq* p=selSeq; //previous selection
   selSeq=seq;
   if (p!=NULL) drawSeq(p);
   drawSeq(selSeq);
   if (showIt)
     makeVisible(selSeq);
   }
}

void FXClView::makeVisible(ClSeq* seq) {
 CRect rect;
 seq->calcViewRect(rect); //the position of the seq rectangle in the view space
 //we'll try to determine how to scroll to get this sequence shown in the middle of the screen
 int yscroll =  getContentY();
 int xscroll =  getContentX();
 int cx=width/2 - xscroll; //view's x position of the center of the screen
 int cy=height/2 - yscroll; //view's y position of the center
 scrollBy(rect.ix()-cx,rect.iy()-cy);
}

// Released button
long FXClView::onLeftBtnRelease(FXObject*,FXSelector,void* ptr){
  FXEvent* event=(FXEvent*)ptr;
  if (!isEnabled()) return 0;
  ungrab();
  flags|=FLAG_UPDATE;
  flags&=~FLAG_PRESSED;
  // First chance callback
  if(target && target->handle(this,MKUINT(message,SEL_LEFTBUTTONRELEASE),ptr)) return 1;
  // No activity
  //if(!(flg&FLAG_PRESSED) && !(options&LIST_AUTOSELECT)) return 1;

  // Scroll to make item visibke
  //makeItemVisible(current);

  // Generate clicked callbacks, but only if a Seq was clicked:
  int col;
  ClSeq* seq=getSeqAt(event->win_x, event->win_y, col);
  if (btnSeq==NULL)
     selectSeq(NULL); //deselect
   else
     if (btnSeq==seq)
        selectSeq(btnSeq);
    if(event->click_count==1){
        handle(this,MKUINT(0,SEL_CLICKED),(void*)btnSeq);
        }
      else if(event->click_count==2){
        handle(this,MKUINT(0,SEL_DOUBLECLICKED),(void*)btnSeq);
        }
      else if(event->click_count==3){
        handle(this,MKUINT(0,SEL_TRIPLECLICKED),(void*)btnSeq);
        }
     handle(this,MKUINT(0,SEL_COMMAND),(void*)btnSeq);

  return 1;
 }

// Clicked in list
long FXClView::onClicked(FXObject*,FXSelector,void* ptr){
  return target && target->handle(this,MKUINT(message,SEL_CLICKED),ptr);
  }


/*
 void FXClView::initGrpColors(int num) {
 GFREE(grpColors);
 GMALLOC(grpColors, sizeof(FXColor)*num);
 FXColor c=FXRGB(0x80, 0,0);

 int step=255/(num+1);
 if (step==0) step=10;
 for (int i=0;i<num;i++) {
  c= modhls(c, step, 0, 0);
  if (i % 2)
     grpColors[num-i-1]=c;
    else
     grpColors[i]=c;
  }
} */

bool FXClView::setSeqGrp(FXString& seq, int grp) {
 //seqs must be sorted by name
 ClSeq test(this, seq.text(), 200);
 int sid=-1;
 if ((sid=seqlist->IndexOf(&test))<0)
     return false;
 seqlist->Get(sid)->group=grp;
 return true;
 }

void FXClView::initColors(FXColor* seqcolors, FXColor* basecolors,
      FXColor* matchcolors, FXColor* ctgquals, FXColor* grpcolors, FXColor* grpgapcolors) {
  gridColor=getApp()->getBaseColor();
  gridColorH=makeHiliteColor(makeHiliteColor(makeHiliteColor(gridColor)));
  gridColorSh=makeShadowColor(makeShadowColor(makeShadowColor(gridColor)));
  seqColors=seqcolors;
  baseColors=basecolors;
  matchColors=matchcolors;
  ctgQuals=ctgquals;
  grpColors[0]=grpcolors;
  grpColors[1]=grpgapcolors;
  /*
  //matchColors[16]=makeHiliteColor(makeHiliteColor(fxcolorfromname("Red")));
  matchColors[16]=mismatch;
  ctgQuals[16]=getApp()->getBaseColor();
  ctgQuals[15]=0xFF00FFFF;
  matchColors[15]=fxcolorfromname("darkBlue");
  for (int i=0;i<15;i++) {
    matchColors[14-i]=modhls(matchColors[15], 0, 10*i, -5*i);
    ctgQuals[14-i]=modhls(ctgQuals[15], 0, -4*i, 0);
    }

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
  /*
  //-- base sequence color:
  seqColors[0]=fxcolorfromname("darkBlue");
  // Hex Colors are BGR style !
  //-- highlight border:
  seqColors[1]=makeHiliteColor(fxcolorfromname("Gray50"));
  //--shadow border:
  seqColors[2]=makeShadowColor(fxcolorfromname("Gray20"));
  //seqColors[2]=0x0020A040;
  seqColors[3]=getApp()->getSelforeColor();
  //-- clipping color:
  //seqColors[4]=fxcolorfromname("Gray50"); //trimmed ends
  seqColors[4]=fxcolorfromname("MistyRose3"); //trimmed ends
  seqColors[5]=fxcolorfromname("Green");
  seqColors[6]=fxcolorfromname("Yellow");
  //-- NP color:
  seqColors[7]=fxcolorfromname("darkRed");
  //seqColors[7]=0x0000A000; //greenish
  //-- inter-segment gap line:
  //seqColors[8]=fxcolorfromname("Red");
  //seqColors[8]=0x006060F0;
  seqColors[8]=fxcolorfromname("LightSkyBlue3");
  //-- internal clipping color:
  //seqColors[9]=fxcolorfromname("DarkSalmon");
  seqColors[9]=fxcolorfromname("MistyRose1");
  //-- splice site consensus
  seqColors[10]=fxcolorfromname("Red"); //hard
  seqColors[11]=fxcolorfromname("DarkPink"); //soft

  baseColors[0]=fxcolorfromname("White"); //font color
  baseColors[1]=fxcolorfromname("DarkGreen");  //A
  baseColors[2]=fxcolorfromname("Blue");   //C
  baseColors[3]=fxcolorfromname("Orange"); //G
  baseColors[4]=fxcolorfromname("Red");    //T
  baseColors[5]=fxcolorfromname("Gray20");  // N,*, other codes
  */
  }

// Cleanup
FXClView::~FXClView(){
  getApp()->removeTimeout(this,ID_TIPTIMER);
  delete contig;
  delete font;
  font=(FXFont*)-1;
  //GFREE(grpColors);
  delete seqfont;
  delete backbuf;
  delete seqlist;
  //delete seqaligns;
}

// Create window
void FXClView::create(){
  FXScrollArea::create();
  font->create();
  seqfont->create();
  seqfntH=seqfont->getFontHeight()+2;
  seqfntW=seqfont->getFontWidth();
  if (backbuf==NULL) {
      backbuf=new FXImage(getApp(),
           NULL, IMAGE_SHMI|IMAGE_SHMP, width,height);
      }
  backbuf->create(); //needed?
  RecalcContent();
  }


// Detach window
void FXClView::detach(){
  FXScrollArea::detach();
  font->detach();
  seqfont->detach();
  backbuf->detach();
  }


// Move content
void FXClView::scrollBy(int dx,int dy){
   vertical->setPosition(-pos_y+dy);
   pos_y=-vertical->getPosition();
   horizontal->setPosition(-pos_x+dx);
   pos_x=-horizontal->getPosition();
   update();
  }

void FXClView::position(FXint x, FXint y, FXint w,FXint h) {
  //fprintf(stderr, "****** ::position() called!\n");
  bool resized=false;
  if(w!=width || h!=height){
    if (backbuf!=NULL) {
        backbuf->resize(w,h);
        resized=true;
    	}
  }
  FXScrollArea::position(x,y,w,h);
}

inline FXbool IntersectRect(const CRect& r1, const CRect& r2) {
 return
    !((r1.x + r1.w <= r2.x) || (r2.x + r2.w <= r1.x)
	|| (r1.y + r1.h <= r2.y) || (r2.y + r2.h <= r1.y));
 }


inline FXbool FXIntersectRect(const FXRectangle& r1, const FXRectangle& r2) {
 return
    !((r1.x + r1.w <= r2.x) || (r2.x + r2.w <= r1.x)
	|| (r1.y + r1.h <= r2.y) || (r2.y + r2.h <= r1.y));
 }

//computes the intersection of two rectangles
CRect RIntersection(CRect& r1, CRect& r2) {
 CRect R; //this is the result;
 if (IntersectRect(r1,r2)) {
     R.x = FXMAX (r1.x, r2.x);
     R.y = FXMAX (r1.y, r2.y);
     R.w = FXMIN(r1.x+r1.w, r2.x+r2.w) - R.x;
     R.h = FXMIN(r1.y+r1.h, r2.y+r2.h) - R.y;
     }
 return R;
}

//computes the intersection of two rectangles
FXRectangle FXRIntersection(FXRectangle& r1, FXRectangle& r2) {
  FXRectangle R(0,0,0,0); //this is the result;
  if (FXIntersectRect(r1,r2)) {
     R.x = FXMAX (r1.x, r2.x);
     R.y = FXMAX (r1.y, r2.y);
     R.w = FXMIN(r1.x+r1.w, r2.x+r2.w) - R.x;
     R.h = FXMIN(r1.y+r1.h, r2.y+r2.h) - R.y;
     }
  return R;
  }

inline void copyFXRect(CRect& rect, const FXRectangle& src) {
 rect.x=src.x;rect.y=src.y;
 rect.w=src.w;rect.h=src.h;
}

inline void offsetRect(CRect& rect, int offsx, int offsy) {
 rect.x+=offsx;
 rect.y+=offsy;
}

inline void offsetFXRect(FXRectangle& rect, int offsx, int offsy) {
 rect.x+=offsx;
 rect.y+=offsy;
}

//int ClSeq::xTran(int x) {
 //takes x as a sequence coordinate and computes the current view
//}


void ClSeq::calcViewRect(CRect& seqrect) {
 seqrect.x= view->seqW*(double)(cl_xpos - view->XLeft)+(double)CL_BORDER;
 seqrect.y=(view->seqH + view->vSpace)*cl_ypos + view->gridH+CL_BORDER;
 seqrect.w=view->seqW*(double)len+2; //?? is this 2 pixel correction really needed?
 seqrect.h=view->seqH;
}

void ClSeq::calcSegView(int segidx, CRect& segrect) {
//compute a seq segment's position
 //within the view
 segrect.x= view->seqW*(double)(segs[segidx].xL - view->XLeft)+(double)CL_BORDER;
 segrect.y=(view->seqH + view->vSpace)*cl_ypos + view->gridH+CL_BORDER;
 segrect.w=view->seqW*(double)(segs[segidx].xR-segs[segidx].xL+1)+2; //?? is this 2 pixel correction really needed?
 segrect.h=view->seqH;
}

/*FXbool ClSeq::InRect(const CRect& rect, CRect& seqrect) {
 //for testing if repainting is needed for this object
 //rect is in pixels, but already offset to the whole cluster area
 //also computes the real-world coordinates into seqrect
 //seqH and seqW must already be computed (acording to scale!)
 calcViewRect(seqrect);
 //to account for rounding errors, enlarge the rectangle a bit:
 //seqrect.x--;seqrect.y--;
 //seqrect.w++;seqrect.h++;
 //return IntersectRect(seqrect,rect);

 return (view->canShowSeq)?(seqrect.x+seqrect.w+view->seqfntW >rect.x
                   && seqrect.x-view->seqfntW<rect.x+rect.w):
     (seqrect.x+seqrect.w>rect.x && seqrect.x<rect.x+rect.w);
} */

static int nt2color(char c) {
  switch(c) {
    case 'a':
    case 'A':return 1;
    case 'c':
    case 'C':return 2;
    case 'g':
    case 'G':return 3;
    case 't':
    case 'T':return 4;
    }
  return 5;
  }


ClSeq* FXClView::getSeqAt(int x, int y, int& col) {
//locate the sequence that's on screen at coords x,y
//locate the row:
x-=getContentX();
y-=getContentY();
col=-1;
int r0=(y-gridH-CL_BORDER)/(seqH+vSpace); //layout row
int colno=(int)floor(((double)(x-CL_BORDER))/seqW);
if (colno>=0 && colno<columns.Count()) col=colno;
if (r0<0 || r0>=rows.Count())
    return NULL;

if ((y-gridH-CL_BORDER) % (seqH+vSpace) > seqH)
    return NULL; //hit between rows

int xchpos =  colno + XLeft; //layout column
//for each of the sequences on this row (they are sorted by cl_xpos)
for (int i=0;i<rows[r0]->seqs.Count();i++) {
 ClSeq* s=rows[r0]->seqs[i];
 if (xchpos>=s->cl_xpos && xchpos<=s->cl_xpos+s->len-1)
    return s;
 }
return NULL;
}


void FXClView::drawTriangle(FXDCWindow& dc, FXColor color, FXint l,FXint t,FXint r,FXint b){
  FXPoint points[3];
  int m=(t+b)/2;
  dc.setForeground(color);
  points[0].x=l;
  points[0].y=t;
  points[1].x=l;
  points[1].y=b;
  points[2].x=r;
  points[2].y=m;
  dc.fillPolygon(points,3);
 }


void FXClView::drawSeq(ClSeq* seq) {
 //redraw a particular sequence
 CRect rect;
 seq->calcViewRect(rect);
 offsetRect(rect, getContentX(), getContentY());
 update(rect.ix()-1,rect.iy()-1,rect.iw()+1,rect.ih()+1);
 /* directly on the surface dc ?
    seq->Paint(dc, rect);
 */
 }

void FXClView::paintSeq(FXDCWindow* dc, FXRectangle& ClpRct,
						ClSeq* seq, ColorStyle colorstyle) {
   //ClpRct is the rectangle that needs to be repainted
   //
   int xp = -getContentX();
   int yp = -getContentY();
   int margin=CL_BORDER<<1;
   int endx=ClpRct.x + ClpRct.w+1; //end X coordinate for the "dirty" rectangle
   //--first (leftmost) column of the layout which needs painting:
   int pLcol=(int)floor(((double)(ClpRct.x+xp-CL_BORDER))/seqW);
   //--last (rightmost) column of the layout which needs painting:
   int pRcol=(int)floor(((double)(endx+xp-CL_BORDER))/seqW);
   int dXL=(int)ClpRct.x-margin; //left px coord of dirty rectangle
   int dXR=(int)ClpRct.x+(int)ClpRct.w+margin; // right px coord of dirty
   int lastXend=INT_MIN; //for inter-segment connector drawing only
   int gapLx=0;
   int gapRx=0;
   char gapspL=0;//is there a splice site on the left side of the gap?
   char gapspR=0;//is there a splice site on the right side of the gap?
   FXColor gapcolor=seq->group>=0 ? grpColors[1][seq->group] : seqColors[8];
   //bool showSeq=(canShowSeq && !seq->sequence.empty());
   bool showSeq=canShowSeq && seq->sequence!=NULL;
  for (int iseg=0;iseg<seq->numsegs;iseg++) {
   CRect segrect;
   seq->calcSegView(iseg,segrect);
   offsetRect(segrect, -xp, -yp);
   FXint ix, iw; //last (rounded) pixel
         //coordinates used for rectangle painting
   // ---
   ix=0;iw=0;
   //-- first (leftmost) column of the segment (could be clipped)
   //      (it may not need painting)
   int lcol=seq->segs[iseg].xL-XLeft;
   //-- last (rightmost) column (base position) of the segment
   //      which is also in the dirty rect (so needs painting)
   int endcol=FXMIN(pRcol+1, seq->segs[iseg].xR+1-XLeft);

   //-- maxcol= maximum non-clipped column (letter)
   int maxcol=seq->segs[iseg].xR+1-XLeft-seq->segs[iseg].clipR;

   //compute seg_x and seg_w to be right outside the clipping range
   //when segrect.x and/or segrect.w are way out of drawable range
   int seg_w = iround(segrect.w);
   int seg_x = (int)floor(segrect.x);
   int seg_xr = seg_x+seg_w; //right end of segment (pixels)
   int seg_y = (int)floor(segrect.y);
   gapRx=seg_x;
   gapspR=seq->segs[iseg].splL;
   bool Lend_on_screen = true;
   bool Rend_on_screen = true;
   if (seg_xr > dXR) {
      //right end is outside the dirty area
      if (seg_x > dXR) {
         if (lastXend>INT_MIN) {
           paintGap(dc, gapcolor, ClpRct, seg_y, lastXend, dXR, gapLx, gapRx, gapspL, gapspR);
           lastXend=INT_MIN;
           }
         gapLx=(int)(segrect.x+segrect.w);
         gapspL=seq->segs[iseg].splR;
         continue; //skip this segment
         }
      Rend_on_screen=false;
      seg_xr = dXR;
      }
   if (seg_x < dXL) {
      //--left end is far away, outside the dirty area
      if (seg_xr < dXL) {
        lastXend=dXL;
        gapLx=(int)(segrect.x+segrect.w);
        gapspL=seq->segs[iseg].splR;
        continue; //skip this segment, it's not visible
        }
      Lend_on_screen=false;
      seg_x=dXL;
      }
   seg_w=seg_xr-seg_x;
   int xcol=pLcol-1;
   if (xcol<lcol) xcol=lcol;
   //xcol = first column of the segment
   //                in the dirty rectangle
   int xtocol=pRcol+1;
   if (xtocol>maxcol) xtocol=maxcol;
   ///xtocol = last non-clipped column of the segment
   //                in the directy rectangle
   //-----------------------------------------------------
   dc->setFont(seqfont);
   int txtY=0;
   if (showSeq) {
	 txtY=seg_y+iround((seqH-seqfntH)/2)+seqfont->getFontAscent()+1;
	 }
   //-- x is set to the px coordinate of the leftmost column to be painted
   double x=seqW*(double)xcol+CL_BORDER-xp;
   if (lastXend>INT_MIN) {
     if (xcol==lcol)
         paintGap(dc, gapcolor, ClpRct, seg_y, lastXend, iround(x), gapLx, gapRx,
                                                   gapspL, gapspR);
     lastXend=INT_MIN;
     }
   double newx=0; //-- keep adjusting newx to the next x
   //--- draw clipL:
   if (seq->segs[iseg].clipL>0 && lcol+seq->segs[iseg].clipL>xcol) {
     newx=seqW*((double)(lcol+seq->segs[iseg].clipL))+CL_BORDER-xp;
     //dc->setForeground(seqColors[4]);
     if (iseg == 0) dc->setForeground(seqColors[4]);
               else dc->setForeground(seqColors[9]);
     dc->fillRectangle(iround(x), iround(segrect.y), iround(newx-x)+1,seqH);
     //draw the text too:
     int tlen=lcol+seq->segs[iseg].clipL-xcol;
     if (showSeq) {
             dc->setForeground(baseColors[0]);
             dc->drawText((int)floor(x), txtY,
					// seq->sequence.text()+xcol-lcol, tlen);
					seq->sequence+seq->segs[iseg].seqofs+xcol-lcol, tlen);
             } //letter draw
     x=newx; // x is prepared for the next segment to be painted
     xcol+=tlen;
     }
  //=== here: x and xcol  are set correctly for the first used lettter
  if (xcol<XRight-XLeft) {
   int cth;
   int cmsm;
   int w;
   int cidx;
   int xlen=1; //run length for cth/cmsm
   int txtlen=1;
   int txtpos;
   int txtcol=xcol; //same as x but for columns
   bool mismatch=false;
   switch (colorstyle) {
    case FXClView::csDensity:
       cth=columns[xcol]->thickness;
       cmsm=columns[xcol]->mismatches;
       cidx=iround((15.0*(cth-cmsm))/maxThickness);
       while (xcol<xtocol) {
            xcol++;
            cth=columns[xcol]->thickness;
            cmsm=columns[xcol]->mismatches;
            newx=seqW*((double)xcol)+CL_BORDER-xp;
            mismatch= (seq->sequence!=NULL && columns[xcol]->letter!=0 &&
                  columns[xcol]->letter!=toupper(seq->sequence[seq->segs[iseg].seqofs+xcol-lcol]));
            int newcidx= mismatch? MISMATCH_CLRIDX : iround((15.0*(cth-cmsm))/maxThickness);
            if (cidx!=newcidx || xcol==xtocol) {
                //end previous run:
                FXColor color=matchColors[cidx];
                dc->setForeground(color);
                ix=iround(x);
                if (xcol==xtocol && Rend_on_screen) {
                  w=seg_xr-ix;
                  }
                 else {
                  w = iround(seqW*xlen);
                  if (w<iround(newx-x+1)) w++;
                  }
                dc->fillRectangle(ix,segrect.iy(),w,seqH);
                //now paint the text too in that area:
                if (showSeq) {
                    dc->setForeground(baseColors[0]);
                    dc->drawText((int)floor(x), txtY,
						 //seq->sequence.text()+txtcol-lcol, txtlen);
						 seq->sequence+seq->segs[iseg].seqofs+txtcol-lcol, txtlen);
                  } //letter draw
                xlen=1;
                txtlen=1;
                x=newx;
                txtcol=xcol;
                cidx=newcidx;
                }
              else {
                xlen++;
                txtlen++;
                }
             } //while column
    break;
   case FXClView::csBaseColor:
     txtlen=xtocol - xcol;
     txtcol=xcol;
     txtpos=(int)floor(x);
     while (xcol<xtocol) {
       newx=seqW*((double)(xcol+1))+CL_BORDER-xp;
       dc->setForeground(baseColors[nt2color(seq->sequence[xcol-lcol+seq->segs[iseg].seqofs])]);
       w = iround(newx-x);
       if (w<iround(newx-x+1)) w++;
       dc->fillRectangle(iround(x),segrect.iy(), w, seqH);
       x=newx;
       xcol++;
       }
    goto DRAW_TEXT;
    break;

   case FXClView::csDefault:
    //this also accounts for the subcoloring or ET coloring
    txtlen=xtocol - xcol;
    txtcol=xcol;
    txtpos=(int)floor(x);
    newx=seqW*((double)xtocol)+CL_BORDER-xp;
    if (seq->segs[iseg].clipR==0) {
       w = seg_xr-iround(x);
       }
      else {
       w=iround(seqW*(xtocol-xcol));
       if (w<iround(newx-x+1)) w++;
       }
    //if (w<=0) w=1;
    if (w>0) {
    if (seq->group>=0) {
       dc->setForeground(grpColors[0][seq->group]);
       }
      else {
       dc->setForeground(seq->isET?seqColors[7]:seqColors[0]);
       }
     dc->fillRectangle(iround(x),seg_y, w, seqH);
    //seq drawing:
    DRAW_TEXT:
     if (showSeq) {
       dc->setForeground(seqColors[3]);
	   dc->drawText(txtpos, txtY,
		   //seq->sequence.text()+txtcol-lcol, txtlen);
		   seq->sequence+seq->segs[iseg].seqofs+txtcol-lcol, txtlen);
       xcol+=txtlen;
       } //letter draw
    }
    x=newx;
    break;
   }  //switch colorstyle
  } //xcol within active range

  if (Rend_on_screen) lastXend=seg_xr;
   /* draw right end clipping ---- */
  if (seq->segs[iseg].clipR>0 && maxcol<endcol) {
         //dc->setForeground(seqColors[4]);
         if (iseg == (seq->numsegs - 1)) dc->setForeground(seqColors[4]);
                                    else dc->setForeground(seqColors[9]);
         //x=seqW*((double)(maxcol))+CL_BORDER-xp);
         //newx=seqW*((double)(endcol))+CL_BORDER-xp);
         ix=iround(x);
         iw=seg_xr-ix;
         if (iw>0) {
             dc->fillRectangle(ix, seg_y, iw ,seqH);
             //draw the text too:
             int tlen=endcol-maxcol;
             if (showSeq) {//letter drawing
                     dc->setForeground(baseColors[0]);
                     //dc->drawText((int)floor(x), txtY,
                       dc->drawText((int)floor(seqW*((double)(maxcol))+CL_BORDER-xp), txtY,
							//seq->sequence.text()+maxcol-lcol, tlen);
						seq->sequence+seq->segs[iseg].seqofs+maxcol-lcol, tlen);

                     }
             }
         }
   gapLx=(int)(segrect.x+segrect.w);
   gapspL=seq->segs[iseg].splR;
   } // segment loop

//-------
  paintSeqBorder(dc,ClpRct,seq,colorstyle);

}

void FXClView::paintGap(FXDCWindow* dc, FXColor gapcolor, FXRectangle& ClpRct, int y, int x1, int x2,
                               int gapl, int gapr, char splL, char splR) {
  if (x2<x1) return;
  if (x1<ClpRct.x) x1=ClpRct.x;
  if (x2>ClpRct.x+ClpRct.w) x2=ClpRct.x+ClpRct.w;
  //dc->setClipRectangle(ClpRct);
  int splh;
  if (seqH>=2) {
        dc->setLineWidth(2);
        y=iround(y+seqH/2);
        splh=seqH;
        }
    else {
      dc->setLineWidth(1);
      splh=1;
      }
  //dc->setForeground(seqColors[8]);
  dc->setForeground(gapcolor);
  dc->drawLine(x1,y, x2, y);
  dc->setLineWidth(1);
  int splw=iround(seqW*2.00);
  if (splw<1) splw=1;
  if (splh>5) splh=5;
  int sh = (splh>>1);
  if (splL>0) {
    FXRectangle rL(gapl,y-sh,splw,splh);
    FXRectangle ir=FXRIntersection(ClpRct,rL);
    if (ir.w>0) {
      dc->setForeground(splL=='S'?seqColors[10]:seqColors[11]);
      dc->fillRectangle(ir.x,ir.y,ir.w,ir.h);
      }
    }
  if (splR>0) {
    FXRectangle rR(gapr-splw,y-sh,splw,splh);
    FXRectangle ir=FXRIntersection(ClpRct,rR);
    if (ir.w>0) {
      dc->setForeground(splR=='S'?seqColors[10]:seqColors[11]);
      dc->fillRectangle(ir.x,ir.y,ir.w,ir.h);
      }
    }
 }

void FXClView::paintSeqBorder(FXDCWindow* dc, FXRectangle& ClpRct, ClSeq* seq,
                              ColorStyle colorstyle) {
    int xp = -getContentX();
    int yp = -getContentY();
    //int endx=ClpRct.x + ClpRct.w+1; //end X coordinate
    CRect seqrect;
    seq->calcViewRect(seqrect);
    offsetRect(seqrect, -xp, -yp);

    int margin=CL_BORDER<<1;
    //compute seq_x and seq_w to be right outside the clipping range
    //when seqrect.x and/or seqrect.w are way out of drawable range
    bool Lend_on_screen = true;
    bool Rend_on_screen = true;
    int seq_w = iround(seqrect.w);
    int seq_x = (int)floor(seqrect.x);
    int seq_xr = seq_x+seq_w;
    int seq_y = (int)floor(seqrect.y);
    if (seq_xr > (int)ClpRct.x+(int)ClpRct.w+margin) {
       //right end is far away, outside the clipping area
       Rend_on_screen=false;
       seq_xr = (int)ClpRct.x+(int)ClpRct.w+margin;
       }
    if (seq_x < (int)ClpRct.x-margin) {
       //left end is far away, outside the clipping area
       Lend_on_screen=false;
       seq_x=(int)ClpRct.x-margin;
       }
    seq_w=seq_xr-seq_x;
 //
  //FXColor color=seqColors[4];
  FXColor color=0xFF00F000;
  if (seqW>=CL_BORDER) { //draw triangle outside, we have room for that
          if (seq->reversed==0) {
               //draw on the right
               if (Rend_on_screen)
                 drawTriangle(*dc, color,
                      (int)(seqrect.x+seqrect.w), seq_y-1,
                      (int)(seqrect.x+seqrect.w)+CL_BORDER, seq_y+seqH);
               }
             else {
              //draw on the left
              if (Lend_on_screen)
               drawTriangle(*dc, color,
                   seqrect.ix(), seq_y-1,
                   seqrect.ix()-CL_BORDER, seq_y+seqH);
               }
            }
           else { //scaled down: draw the triangles inside
            int tw=iround(seqW);
            if (tw<3)
                  tw=3; //triangle width
            if (seqH>3 && seqrect.w>3) { //otherwise don't draw anything, it's too small scale
              if (seq->reversed==0) {
               //draw on the right
               if (Rend_on_screen)
                  drawTriangle(*dc, color,
                     (int)(seqrect.x+seqrect.w-tw)-2, seq_y+2,
                     (int)(seqrect.x+seqrect.w)-2, seq_y+seqH-2);
               }
             else  { //draw on the left
               if (Lend_on_screen)
                 drawTriangle(*dc, color,
                   seqrect.ix()+tw+2, seq_y+2,
                   seqrect.ix()+2, seq_y+seqH-2);
               }
             }
            }

  //border bevel/selection drawing
  //int endX = Rend_on_screen ? seq_x+seq_w:xlast;
  FXRectangle seqclp(seq_x, seq_y, seq_w, seqH);
  seqclp=FXRIntersection(seqclp, ClpRct);
  bool seqclpSet=false;
  /*
  if (seq->group>=0 && colorstyle==csDefault) { //group border coloring
       dc->setClipRectangle(seqclp);
       seqclpSet=true;
       dc->setLineWidth(4);
       dc->setForeground(makeHiliteColor(grpColors[seq->group]));
        dc->drawLine(seq_x,seq_y, seq_x+seq_w, seq_y);
        //left end
        if (Lend_on_screen)
          dc->drawLine(seq_x, seq_y, seq_x, seq_y+seqH);
       dc->setForeground(makeShadowColor(grpColors[seq->group]));
       dc->drawLine(seq_x,seq_y+seqH, seq_x+seq_w, seq_y+seqH);
        if (Rend_on_screen)
         dc->drawLine(seq_x+seq_w, seq_y, seq_x+seq_w, seq_y+seqH);
       dc->setLineWidth(1);
       } */
    // else if (colorstyle==csDefault) {
       /*dc->setClipRectangle(seqclp);
       seqclpSet=true;
       dc->setLineWidth(2);
       dc->setForeground(seqColors[1]);
        dc->drawLine(seq_x,seq_y, seq_x+seq_w, seq_y);
        //left end
        if (Lend_on_screen)
          dc->drawLine(seq_x, seq_y, seq_x, seq_y+seqH);
       dc->setForeground(seqColors[2]);
       dc->drawLine(seq_x,seq_y+seqH, seq_x+seq_w, seq_y+seqH);
        if (Rend_on_screen)
         dc->drawLine(seq_x+seq_w, seq_y, seq_x+seq_w, seq_y+seqH);
       dc->setLineWidth(1);
       */
     //  }
  if (selSeq==seq) {//draw selection rectangle
     if (!seqclpSet) {
        dc->setClipRectangle(seqclp);
        seqclpSet=true;
        }
     dc->setForeground(seqColors[6]);
     dc->setLineWidth(4);
     dc->drawRectangle(seq_x, seq_y, seq_w, seqH);
     dc->setLineWidth(1);
     }

  if (seqclpSet)
     dc->setClipRectangle(ClpRct);
}

void FXClView::getVisibleRegion(FXint& cstart, FXint& cend) {
	//this must be called AFTER seqW was updated
	cstart = XLeft+(FXint)floor(((double)(-getContentX()-CL_BORDER))/seqW);
	if (cstart<XLeft) cstart=XLeft;
	cend = XLeft+(FXint)iround(((double)(getVisibleWidth()-getContentX()-CL_BORDER))/seqW)-1;
	if (cend>XRight) cend=XRight;
	//fprintf(stderr, "cstart=%d, cend=%d\n", cstart, cend);
	fflush(stderr);
}

void FXClView::paintGrid(FXDCWindow* dc, FXRectangle& gridR) {
    dc->setForeground(gridColor);
    dc->fillRectangle(gridR.x,gridR.y,
                    gridR.w,gridR.h);
    int endx=gridR.x + gridR.w+1;
    int xp = -getContentX();
    //gridStep (in columns) must be computed after each zooming
    int xcol=(int)floor(((double)(gridR.x+xp-CL_BORDER))/seqW);
    int startcol=xcol-gridStep;
    int xtocol=(int)floor(((double)(endx+xp-CL_BORDER))/seqW);
    int endcol=xtocol+gridStep;
    for (; startcol<endcol && startcol%gridStep !=0 ; startcol++);

    int hf=showContig ? font->getFontHeight()+4: gridH;

    if (startcol<endcol) { //there is something to redraw here:
       dc->setFont(font);
       dc->setLineWidth(1);
       for (double x=seqW*(double)startcol+CL_BORDER-xp; x<=endx; x+=seqW*gridStep) {
         int col = iround(((double)(x+xp-CL_BORDER))/seqW)+XLeft;
         dc->setForeground(gridColorH);
         int xx=iround(x);
         dc->drawLine(xx+1,hf-4,xx+1,hf);
         dc->setForeground(gridColorSh);
         dc->drawLine(xx,hf-4,xx,hf);
         FXString s;
         s.format("%d",col);
         dc->drawText((int)x, hf-6, s.text(), s.length());
         }
      }
    //draw contig if there, with quality values:
    if (showContig && xcol<XRight-XLeft) { //xcol is the starting column
     xcol--;
     if (xcol<0) xcol=0;
     xtocol++;
     if (xtocol>=columns.Count()) xtocol=columns.Count()-1;
     double x=seqW*(double)xcol+CL_BORDER-xp;
     int cth=columns[xcol]->thickness;
     int cmsm=columns[xcol]->mismatches;
     int cidx=(xcol<contig->cl_xpos-XLeft || xcol>contig->cl_xpos-XLeft+contig->len-1) ?
                  NOCTG_CLRIDX : iround((15.0*(cth-cmsm))/maxThickness);
     int xlen=1; //run length for cth/cmsm
     dc->setFont(seqfont);
     int txtcol=xcol; //same as x but for columns
     int txtlen=1;
     while (xcol<xtocol) {
        xcol++;
        cth=columns[xcol]->thickness;
        cmsm=columns[xcol]->mismatches;
        double newx=seqW*((double)xcol)+CL_BORDER-xp;
        int newcidx=(xcol<contig->cl_xpos-XLeft || xcol>contig->cl_xpos-XLeft+contig->len-1) ?
                  NOCTG_CLRIDX : iround((15.0*(cth-cmsm))/maxThickness);
        if (cidx!=newcidx || xcol==xtocol) {
           // change of color; draw the previous rectangle
            FXColor color=ctgQuals[cidx];
            dc->setForeground(color);
            int w = iround(seqW*xlen);
            if (w<iround(newx-x+1)) w++;
            dc->fillRectangle(iround(x),gridH-seqH,w,seqH-2);
            //now paint the text too in that area:
            if (cidx!=NOCTG_CLRIDX && canShowSeq && contig->sequence!=NULL) {
                  dc->setForeground(FXRGB(0x00, 0x00, 0x00));
                  dc->drawText((int)floor(x), gridH-2,
                     contig->sequence+txtcol-(contig->cl_xpos-XLeft), txtlen);
              } //letter draw
            x=newx;
            txtcol=xcol;
            xlen=1;
            txtlen=1;
            cidx=newcidx;
            }
          else {
            xlen++;
            txtlen++;
            }
         } //while column
     //now for the last segment:
    } //if showcontig
    dc->setLineWidth(1);
    dc->setForeground(makeShadowColor(gridColor));
    dc->drawLine(gridR.x, gridH-1,
                 gridR.x+gridR.w, gridH-1);
    dc->setForeground(0xFF000000);
    dc->drawLine(gridR.x, gridH,
                 gridR.x+gridR.w, gridH);
}

// Widget has been resized
long FXClView::onConfigure(FXObject*,FXSelector,void*){
  resize(getWidth(), getHeight());
  return 1;
  }

// Draw item list
long FXClView::onPaint(FXObject*,FXSelector,void* ptr){
  FXEvent* event=(FXEvent*)ptr;
  //event holds the coordinates for the invalidated rectangle only!
  //avoiding unneccessary repainting of the whole canvas
  //cool!
  FXDCWindow dc(this,event);
  FXRectangle clip(event->rect);
  FXRectangle gridR(0,0,getVisibleWidth()+1, gridH+1);
  gridR=FXRIntersection(clip,gridR);
  //rectangle to be converted to global, "content" coordinates
  //back buffer stuff
  FXDCWindow *bufdc=new FXDCWindow(backbuf);
  if (gridR.h>0) {
    //paint the top grid:
    //gridR is 0,0 based
    paintGrid(bufdc, gridR);
    clip.h=clip.y+clip.h-gridH;
    clip.y=gridH+1;
    }

  if (rows.Count()>0) {
    /*
    if (clip.y<gridH) {//?!?!? strange bug sometimes this is just out of range (negative coords!)
       //GMessage("gridR non empty (gridR.h=%d, cliprect set to (%d, %d, %d, %d)\n",
       //          gridR.h, clip.x, clip.y, clip.w, clip.h);
       //clip.y=gridH+1;
       }
    */
    bufdc->setClipRectangle(clip);

    bufdc->setForeground(backColor);
    bufdc->fillRectangle(clip.x,clip.y,
                    clip.w,clip.h);
    //paint this by row - only visible rows
    //locate first and last visible rows:

    CRect rect(clip);
    //rect is offset to be positioned in the whole content space:
    offsetRect(rect,-getContentX(),-getContentY());
     //start-end columns to check
    int c1=(int)floor(((double)(rect.x-CL_BORDER))/seqW)-1;
    int c2=(int)floor(((double)(rect.x+rect.w-CL_BORDER))/seqW)+1;

    int r1=iround((rect.y-gridH-CL_BORDER)/(seqH+vSpace))-1;
    if (r1<0) r1=0;
    int r2=iround((rect.y+rect.h-gridH-CL_BORDER)/(seqH+vSpace))+1;
    if (r2>=rows.Count()-1) r2=rows.Count()-1;
    ColorStyle cs=colorStyle;
    if (seqW<1)
         cs=FXClView::csDefault;
    for (int i=r1; i<=r2; i++)
      for (int j=0;j<rows[i]->seqs.Count(); j++) {
         ClSeq* sq=rows[i]->seqs[j];
         if (sq->cl_xpos-XLeft>c2) continue;
         if (sq->cl_xpos+sq->len-XLeft<c1) continue;
         paintSeq(bufdc,clip, sq, cs);
         }
    if (selCol>=0) { //draw column selection
      bufdc->setClipRectangle(event->rect);
      int cx= iround(selCol*seqW)+CL_BORDER+getContentX();
      bufdc->setLineWidth(1);
      bufdc->setFillStyle(FILL_OPAQUESTIPPLED);
      bufdc->setStipple(STIPPLE_GRAY,cx,getContentY());
      //bufdc->setFunction((FXFunction)BLT_SRC_XOR_DST);
      bufdc->setForeground(getTextColor());
      bufdc->setBackground(seqColors[6]);
      bufdc->drawLine(cx,0, cx, getVisibleHeight() );
      bufdc->drawLine(cx+(int)seqW, 0, cx+(int)seqW, getVisibleHeight());
      //bufdc->drawRectangle(cx,0,seqW,getVisibleHeight());
      bufdc->setFillStyle(FILL_SOLID);
      bufdc->setStipple(STIPPLE_NONE);
      }
    delete bufdc; //? This is needed to actually flush the changes!?!?
    }
  else {
    dc.setForeground(backColor);
    dc.fillRectangle(clip.x,clip.y,
                    clip.w,clip.h);
    }
  //now draw the buffer:
  // Set font
  //dc.setForeground(fxcolorfromname("LightGray"));
  //dc.fillRectangle(0,0,
  //                 getWidth(),getHeight());
  //dc.drawImage(backbuf,event->rect.x,event->rect.y);
  dc.drawArea(backbuf, event->rect.x,event->rect.y, event->rect.w,event->rect.h,
              event->rect.x,event->rect.y);
  return 1;
  }

FXint FXClView::getContentWidth() {
  //if(flags&FLAG_RECALC) recompute();
  return content_w;
  }

// Determine content height of list
FXint FXClView::getContentHeight(){
  //if(flags&FLAG_RECALC) recompute();
  return content_h;
  }

// Change the font
void FXClView::setFont(FXFont* fnt){
  if(!fnt){ fxerror("%s::setFont: NULL font specified.\n",getClassName()); }
  if(font!=fnt){
    font=fnt;
    recalc();
    update();
    }
  }

// Set text color
void FXClView::setTextColor(FXColor clr){
  if(clr!=textColor){
    textColor=clr;
    update();
    }
  }

// Save data
void FXClView::save(FXStream& store) const {
  FXScrollArea::save(store);
  store << font;
  store << textColor;
  }
// Load data
void FXClView::load(FXStream& store){
  FXScrollArea::load(store);
  store >> font;
  store >> textColor;
  }

//-- compare functions for sequences
int compareXPos(void* p1, void* p2) {
 int c1=((ClSeq*)p1)->cl_xpos;
 int c2=((ClSeq*)p2)->cl_xpos;
 return (c1>c2)?1:((c1<c2)?-1:0);
}

int compareXRoom(void* p1, void* p2) {
 int c1=((ClSeq*)p1)->view->XRight-((ClSeq*)p1)->cl_xpos+((ClSeq*)p1)->len;
 int c2=((ClSeq*)p2)->view->XRight-((ClSeq*)p2)->cl_xpos+((ClSeq*)p2)->len;
 return (c1>c2)?-1:((c1<c2)?1:0);
}

int compareYPos(void* p1, void* p2) {
 int c1=((ClSeq*)p1)->cl_ypos;
 int c2=((ClSeq*)p2)->cl_ypos;
 return (c1>c2)?1:((c1<c2)?-1:0);
}

void FXClView::addContig(const char* name, FXint len,
     const char* sequence, FXint offs) {
  contig=new ClSeq(this, name,len, offs, false, 0, 0, sequence);
  //=addSeq(name, len, offs, false, 0,0,sequence);
  if (sequence!=NULL) hasSeqs=true;
  contig->cl_ypos=0;
  if (seqlist->Count()==0) {
    XLeft=offs;
    XRight=contig->cl_xpos+contig->len;
    }
   else {
     if (offs<XLeft) XLeft=offs;
     //also re-check current dimension
     if (contig->cl_xpos + contig->len > XRight)
       XRight=contig->cl_xpos+contig->len; //extreme right position
     }
  }

ClSeq* FXClView::addSeq(const char* name, FXint len, FXint offs,
 bool reversed, int trim_l, int trim_r, const char* sequence,
 unsigned int ins, unsigned int dels) {
 ClSeq* seq=new ClSeq(this, name,len, offs, reversed, trim_l, trim_r, sequence,
                              ins, dels);
 if (sequence!=NULL) hasSeqs=true;
 seqlist->Add(seq);
 //for now:
 seq->cl_ypos=seqlist->Count();
 if (seqlist->Count()==1) {
      XLeft=offs;
      XRight=seq->cl_xpos+seq->len;
      }
   else {
     if (offs<XLeft) XLeft=offs;
     //also re-check current dimension
     if (seq->cl_xpos + seq->len > XRight)
       XRight=seq->cl_xpos+seq->len; //extreme right position
     }
 return seq;
}

ClSeq* FXClView::addSeq(LytSeqInfo* seqinfo, const char* sequence) {
 ClSeq* seq=new ClSeq(this, seqinfo, sequence);
 if (sequence!=NULL) hasSeqs=true;
 seqlist->Add(seq);
 //--------
 seq->cl_ypos=seqlist->Count();
 if (seqlist->Count()==1) {
      XLeft=seqinfo->offs;
      XRight=seq->cl_xpos+seq->len;
      }
   else {
     if (seqinfo->offs<XLeft) XLeft=seqinfo->offs;
     //also re-check current dimension
     if (seq->cl_xpos + seq->len > XRight)
       XRight=seq->cl_xpos+seq->len; //extreme right position
     }
 return seq;
}

void FXClView::RecalcContent(bool re_paint) {
 seqH = iround((scale.sy * (double)seqfntH));
 if (seqH==0) seqH=1;
 vSpace=seqH/4;
 seqW = (scale.sx * (double)seqfntW);
 if (seqW==0) seqW=0.01;
 canShowSeq = (seqfntW<=seqW && seqfntH<=seqH);
 gridH=font->getFontHeight()+4;
 if (contig!=NULL) {
   gridH+=seqH+2;
   showContig=true;
   }
  else showContig=false;
 //getFontWidth() doesn't really work for ps fonts,
 //so I use 40 as the max width of a grid number in pixels
 for (gridStep=10; seqW*gridStep < 40; gridStep+=10) ;
 //gridStep = (seqW<4.0) ? 10 : 10*floor((double)4.0/seqW);
 //recompute the sequence height and vspace, based on scale.sy:
 if (re_paint) {
	 placeScrollBars(width,height);
   //layout();
   update();
   }
}

int qsortcmp(const void* p1, const void* p2) {
 NTColumn* c1 =(NTColumn*)p1;
 NTColumn* c2 =(NTColumn*)p2;
 return (c1->count>c2->count)? -1 : ((c1->count<c2->count)?1:0);
 }

void FXClView::buildLayout() {
 //try to place more than one sequence on a row
 seqlist->setSorted(compareXPos);
 rows.Clear();
 columns.Clear();
 for (int i=0;i<seqlist->Count();i++) {
   ClSeq* seq = seqlist->Get(i);
  //search back for possible room
   int r=-1;
   for (int j=0;j<rows.Count();j++)
	 if (seq->cl_xpos-2 > rows[j]->rpos) {
		  r=j;
		  break;
		  }
  if (r>=0)
	  rows[r]->addSeq(seq, r);
	else
	  rows.Add(new CLayoutRow(seq,rows.Count()));
  }

 //-- analyze columns -----:
  for (int c=XLeft; c<=XRight; c++) {
	ColumnData* column=new ColumnData();
	column->ntdata[0].letter='A';column->ntdata[0].count=0;
	column->ntdata[1].letter='C';column->ntdata[1].count=0;
	column->ntdata[2].letter='G';column->ntdata[2].count=0;
	column->ntdata[3].letter='T';column->ntdata[3].count=0;
	column->ntdata[4].letter='N';column->ntdata[4].count=0;
	columns.Add(column);
	}
 //populate column data based on row content;
 for (int r=0;r<rows.Count();r++) {
	for (int i = 0; i < rows[r]->seqs.Count(); i++) { //for each read
	  ClSeq* cs=rows[r]->seqs[i];
	  if (hasSeqs) {
	   // for each segment
	   for (int j = 0; j < cs->numsegs; j++) {
		 for (int k=cs->segs[j].xL;k<=cs->segs[j].xR;k++) {
			ColumnData* column = columns[k-XLeft];
			column->thickness++;
			char c=toupper(cs->sequence[k-cs->segs[j].xL+cs->segs[j].seqofs]);
			switch (c) {
				 case   0:break;
				 case ' ':column->thickness++;
						  break; //no sequence given
				 case 'A':
				 case 'a':column->ntdata[0].count++;column->thickness++;break;
				 case 'C':
				 case 'c':column->ntdata[1].count++;column->thickness++;break;
				 case 'G':
				 case 'g':column->ntdata[2].count++;column->thickness++;break;
				 case 'T':
				 case 't':column->ntdata[3].count++;column->thickness++;break;
				 case '*':
				 case '-': column->gap++;
						  column->thickness++;
						  break;
				 default: { column->ntdata[4].count++; column->nN++; column->thickness++; }
				 } //switch
			}
		 }
	   } // hasSeqs
	  else {
       for (int j = cs->cl_xpos; j < cs->cl_xpos+cs->len; j++) {
          columns[j-XLeft]->thickness++;
		  }
	   } //no seq
	  }//for each read
	}

 maxThickness=0;
 for (int c=XLeft; c<=XRight; c++) {
	ColumnData* column=columns[c-XLeft];
	if (column->thickness > maxThickness)
		maxThickness=column->thickness;
   //decide what's the letter to use:
	if (hasSeqs) {
	 char letter='\0';
	 column->letter='\0';
	 int max=0;
	 if (column->thickness > 0) {
	   qsort(column->ntdata,4, sizeof(NTColumn), qsortcmp);
	   max=FXMAX(column->ntdata[0].count, column->nN);
	   if (contig!=NULL && contig->sequence!=NULL && c>=contig->cl_xpos) {
		  column->letter=toupper(contig->sequence[c-contig->cl_xpos]);
		  }
	   }
	 if (column->letter) {
		column->mismatches=0;
		bool wasletter=false;
		for (int nt=0;nt<5;nt++) {
		  if (column->ntdata[nt].letter!=column->letter && column->ntdata[nt].letter!='N')
			   column->mismatches+=column->ntdata[nt].count;
		  if (column->ntdata[nt].letter==column->letter && column->ntdata[nt].letter!='N')
				wasletter=true;
		  }
		if (column->letter=='*' || column->letter=='-')
		   column->mismatches=column->thickness - column->gap;
		 else if (!wasletter)
		   column->mismatches=column->thickness - column->nN;

		if (column->mismatches<0) {
		   fprintf(stderr, "Warning: Column %d: %d mism, thickness=%d, gap=%d, nN=%d (%c)\n",
			   c, column->mismatches, column->thickness, column->gap, column->nN, wasletter?'1':'0');

		   for (int nt=0;nt<5;nt++) {
			 fprintf(stderr, "letter %c = %d\n", column->ntdata[nt].letter, column->ntdata[nt].count);
			 }
		   }
		}
	   else {
		if (column->gap>max) {
		   letter = '*';
		   column->mismatches = column->thickness - column->gap;
		   }
		 else {
		   letter = column->ntdata[0].letter;
		   column->mismatches=column->thickness - column->ntdata[0].count;
		   }
		}
	  }//hasSeqs
	} //for each column


/*
  maxThickness=0;
  for (int c=XLeft; c<=XRight; c++) {
	int gap=0, nN=0;
	ColumnData* column=new ColumnData();
   //TODO: use IUPAC ambiguity codes all over
   column->ntdata[0].letter='A';column->ntdata[0].count=0;
   column->ntdata[1].letter='C';column->ntdata[1].count=0;
   column->ntdata[2].letter='G';column->ntdata[2].count=0;
   column->ntdata[3].letter='T';column->ntdata[3].count=0;
   column->ntdata[4].letter='N';column->ntdata[4].count=0;
   for (int r=0;r<rows.Count();r++) {
	   switch (rows[r]->getNuc(c)) {
		 case   0:break;
		 case ' ':column->thickness++;
				  break; //no sequence given
		 case 'A':
		 case 'a':column->ntdata[0].count++;column->thickness++;break;
		 case 'C':
		 case 'c':column->ntdata[1].count++;column->thickness++;break;
		 case 'G':
		 case 'g':column->ntdata[2].count++;column->thickness++;break;
		 case 'T':
		 case 't':column->ntdata[3].count++;column->thickness++;break;
		 case '*':
		 case '-': gap++;
				  column->thickness++;
				  break;
		 default: { column->ntdata[4].count++; nN++; column->thickness++; }
		 } //switch
	   } //for each row
   if (column->thickness > maxThickness)
		maxThickness=column->thickness;
   //decide what's the letter to use:
   if (hasSeqs) {
	 char letter='\0';
	 column->letter='\0';
	 int max=0;
	 if (column->thickness > 0) {
	   qsort(column->ntdata,4, sizeof(NTColumn), qsortcmp);
	   max=FXMAX(column->ntdata[0].count, nN);
	   if (contig!=NULL && !contig->sequence.empty() && c>=contig->cl_xpos) {
		  column->letter=toupper(contig->sequence[c-contig->cl_xpos]);
		  //fprintf(stderr, "%c",column->letter);
		  }
	   }
	 if (column->letter) {
		column->mismatches=0;
		bool wasletter=false;
		for (int nt=0;nt<5;nt++) {
		  if (column->ntdata[nt].letter!=column->letter && column->ntdata[nt].letter!='N')
			   column->mismatches+=column->ntdata[nt].count;
		  if (column->ntdata[nt].letter==column->letter && column->ntdata[nt].letter!='N')
				wasletter=true;
		  }
		if (column->letter=='*' || column->letter=='-')
		   column->mismatches=column->thickness - gap;
		 else if (!wasletter)
		   column->mismatches=column->thickness - nN;

		if (column->mismatches<0) {
		   fprintf(stderr, "Warning: Column %d: %d mism, thickness=%d, gap=%d, nN=%d (%c)\n",
			   c, column->mismatches, column->thickness, gap, nN, wasletter?'1':'0');

		   for (int nt=0;nt<5;nt++) {
			 fprintf(stderr, "letter %c = %d\n", column->ntdata[nt].letter, column->ntdata[nt].count);
			 }
		   }
		}
	   else {
		if (gap>max) {
		   letter = '*';
		   column->mismatches = column->thickness - gap;
		   }
		 else {
		   letter = column->ntdata[0].letter;
		   column->mismatches=column->thickness - column->ntdata[0].count;
		   }
		}
	  }
   columns.Add(column);
   }//for each column
 */

 //sort seqlist by name!!! very important for group loading
 seqlist->setSorted(true);
 RecalcContent();
}


void FXClView::ZoomY(FXdouble zy, int fromY) {
 //save the row num of the center of the screen
 if (zy<=0) zy=0.01;
 int yscroll =  getContentY();
 int cy;
 if (fromY<0) {
	if (selSeq!=NULL) {
	  //zoom about the select sequence row
	  cy=gridH+CL_BORDER+(seqH+vSpace)*selSeq->cl_ypos;
	  }
	 else cy=height/2 - yscroll;
	}
   else
    cy=fromY-yscroll;

 double r0=((double)(cy - gridH-CL_BORDER))/(seqH+vSpace);
 scale.sy = zy;
 RecalcContent(true); //also recalc scrollbar ranges
 //new coordinates:
 //calcLayout();
 //try to offset the view with the number of rows
 reposVScroll(yscroll, r0, cy);
 update();
 }

 void FXClView::ZoomX(FXdouble zx, int fromX) {
  if (zx<=0) zx=0.001;
  int xscroll =  getContentX();
  int cx;
  if (fromX<0) { //zoom using a default center
    if (selCol>=0) {
      //zoom about the intersection nucleotide
      cx=CL_BORDER+(int)floor(seqW*(selCol+1));
      }
     else
       cx=width/2 - xscroll;
     }
   else
     cx=fromX-xscroll; //given center;
  scale.sx = zx;
  double c0=((double)(cx - CL_BORDER))/seqW;
  RecalcContent(true);
  //new coordinates:
  //try to offset the view back to it's original position
  //calcLayout();
  reposHScroll(xscroll, c0, cx);
  if(target) target->tryHandle(this,FXSEL(SEL_CHANGED,message),  NULL);
  update();
 }

/* Move content */
void FXClView::moveContents(FXint x,FXint y){
  FXint dx=x-pos_x;
  FXint dy=y-pos_y;
  pos_x=x;
  pos_y=y;
  FXint v_w=getVisibleWidth();
  scroll(0,gridH+1,v_w,getVisibleHeight()-gridH,dx,dy);
  if (dx!=0) { //scroll grid area too
    scroll(0,0,v_w,gridH,dx,0);
    if(target) target->tryHandle(this,FXSEL(SEL_CHANGED,message),  NULL); //(void*)(FXival)pos);
    }
  }

void FXClView::Zoom(FXdouble zx, FXdouble zy, int fromX, int fromY) {
 if (zy<=0) zy=0.01;
 if (zx<=0) zx=0.001;
 //if (zx<0.01 || zx>4 || zy<0.01 || zy>4) return;
 int yscroll =  getContentY();
 int xscroll =  getContentX();
 int cx;
 if (fromX<0) { //zoom using a default center
    if (selCol>=0) {
      //zoom about the intersection nucleotide
      cx=CL_BORDER+(int)floor(seqW*(selCol+1));
      }
    else
       cx=width/2 - xscroll;
    }
  else
    cx=fromX-xscroll; //given center;
 int cy;
 if (fromY<0) {
    if (selSeq!=NULL) {
      //zoom about the intersection nucleotide
      cy=gridH+CL_BORDER+(seqH+vSpace)*selSeq->cl_ypos;
      }
     else cy=height/2 - yscroll;
    }
   else
    cy=fromY-yscroll;
 scale.sy = zy;
 bool zxChanged = (scale.sx != zx);
 scale.sx = zx;
 double r0=((double)(cy - gridH-CL_BORDER))/(seqH+vSpace);
 double c0=((double)(cx - CL_BORDER))/seqW;
 RecalcContent(true); //recomputes seqW, seqH
 //new coordinates:
 reposVScroll(yscroll, r0, cy);
 reposHScroll(xscroll, c0, cx);
 if (zxChanged && target) target->tryHandle(this,FXSEL(SEL_CHANGED,message),  NULL);
 update();
}


void FXClView::placeScrollBars(FXint vw, FXint vh){
  //register FXint cw,ch;
  register FXint sh_h=0;
  register FXint sv_w=0;

  // Inviolate
  FXASSERT(pos_x<=0 && pos_y<=0);

  // Get content size
  //cw=getContentWidth();
  //ch=getContentHeight();
  // Get dimensions of the scroll bars
  if(!(options&HSCROLLER_NEVER)) sh_h=horizontal->getDefaultHeight();
  if(!(options&VSCROLLER_NEVER)) sv_w=vertical->getDefaultWidth();

  int minH=(seqH+vSpace)*rows.Count()+gridH+CL_BORDER*2;
  int minW=iround(seqW*(double)(XRight-XLeft))+2*CL_BORDER;
  //fprintf(stderr, "XRight=%d, minW=%d\n", XRight, minW);
  bool adjust_h=false, adjust_w=false;
  if (content_h!=minH+4) {
 	 content_h=minH+4;
  }
  if (content_w!=minW+4) {
 	 content_w=minW+4;
  }

  // Should we disable the scroll bars?  A bit tricky as the scrollbars
  // may influence each other's presence.  Also, we don't allow more than
  // 50% of the viewport to be taken up by scrollbars; when the scrollbars
  // take up more than 50% of the available space we simply turn them off.
  if(!(options&(HSCROLLER_ALWAYS|VSCROLLER_ALWAYS)) && (content_w<=vw) && (content_h<=vh)){sh_h=sv_w=0;}
  if(!(options&HSCROLLER_ALWAYS) && ((content_w<=vw-sv_w) || (0>=vh-sh_h-sh_h))) sh_h=0;
  if(!(options&VSCROLLER_ALWAYS) && ((content_h<=vh-sh_h) || (0>=vw-sv_w-sv_w))) sv_w=0;
  if(!(options&HSCROLLER_ALWAYS) && ((content_w<=vw-sv_w) || (0>=vh-sh_h-sh_h))) sh_h=0;

  // Viewport size with scroll bars taken into account
  vw-=sv_w;
  vh-=sh_h;

  // Adjust content size, now that we know about those scroll bars
  if((options&HSCROLLER_NEVER)&&(options&HSCROLLER_ALWAYS)) content_w=vw;
  if((options&VSCROLLER_NEVER)&&(options&VSCROLLER_ALWAYS)) content_h=vh;

  // Furthermore, content size won't be smaller than the viewport
  if(content_w<vw) content_w=vw;
  if(content_h<vh) content_h=vh;
  vscroll=(content_h>vh);
  hscroll=(content_w>vw);
  FXint new_x=pos_x;
  FXint new_y=pos_y;
  // Content size
  if (hscroll) {
	horizontal->show();
    horizontal->setRange(content_w);
    horizontal->setPage(vw);
    horizontal->setPosition(-pos_x);
    new_x=-horizontal->getPosition();
    }
  else horizontal->hide();
  if (vscroll) {
	vertical->show();
    vertical->setRange(content_h);
    // Page size may have changed
    vertical->setPage(vh);
    // Position may have changed
    vertical->setPosition(-pos_y);

    // Get back the adjusted position
    new_y=-vertical->getPosition();
  }
  else vertical->hide();
  // Scroll to force position back into range
  if(new_x!=pos_x || new_y!=pos_y){
    moveContents(new_x,new_y);
    }

  // Read back validated position
  if (hscroll)
	  pos_x=-horizontal->getPosition();
  if (vscroll)
	  pos_y=-vertical->getPosition();

  // Place horizontal scroll bar
  if (hscroll) {
    horizontal->position(0,height-sh_h,width-sv_w,sh_h);
    horizontal->raise();
  }

  // Place vertical scroll bar
  if (vscroll) {
    vertical->position(width-sv_w,0,sv_w,height-sh_h);
    vertical->raise();
  }
  // Place scroll corner
  if (vscroll && hscroll) {
	corner->show();
    corner->position(width-sv_w,height-sh_h,sv_w,sh_h);
    corner->raise();
  }
  else corner->hide();
}


// Recalculate layout

void FXClView::layout(){
  // Place scroll bars
  placeScrollBars(width,height);
  //fprintf(stderr, "layout triggered()!!!!\n");
  //fflush(stderr);
  if(target) target->tryHandle(this,FXSEL(SEL_CHANGED,message),  NULL);
  // Clean
  flags&=~FLAG_DIRTY;
}

// Return visible scroll-area width
FXint FXClView::getVisibleWidth() const {
  FXint w=width;
  if (vscroll)
	  w-=vertical->getWidth();
  return w;
}


// Return visible scroll-area height
FXint FXClView::getVisibleHeight() const {
  FXint h=height;
  if (hscroll) h-=horizontal->getHeight();
  return h;
}

void FXClView::reposVScroll(int yscroll, double r0, int cy) {
  if (vscroll) {
    double r1=(double)((cy - gridH - CL_BORDER))/(seqH+vSpace);
	//vertical->show();
    vertical->setPosition(-(yscroll+iround((r1-r0)*(seqH+vSpace))));
	vertical->raise();
    pos_y=-vertical->getPosition();
  }
  else {
	pos_y=0;
  }
}

void FXClView::reposHScroll(int xscroll, double c0, int cx) {
  if (hscroll) {
    double c1=(double)((cx - CL_BORDER))/seqW;
	//horizontal->show();
    horizontal->setPosition(-(xscroll+iround(((c1-c0)*seqW))));
    horizontal->raise();
    pos_x=-horizontal->getPosition();
  }
  else {
	//horizontal->hide();
	pos_x=0;
  }
}

