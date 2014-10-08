/********************************************************************************
*                                                                               *
*                         (light) Table       W i d g e t                       *
*                                                                               *
*********************************************************************************/
#ifndef FXCLVIEW_H
#define FXCLVIEW_H
#include "GList.hh"
#include <ctype.h>
#include "LayoutParser.h"
#define MISMATCH_CLRIDX 16
#define NOCTG_CLRIDX 16

struct FXTimer;

struct CScale {
 double sx;
 double sy;
};

struct CRect {
 double x,y;
 double w,h;
 CRect() {
  x=0;y=0;w=0;h=0;
  }
 CRect(FXint ax,FXint ay,FXint aw,FXint ah) {
   x=ax;y=ay;w=aw;h=ah;
   }
 CRect(const FXRectangle& r) {
   x=r.x;y=r.y;w=r.w;h=r.h;
   }
 int ix() { return (int)x; };
 int iy() { return (int)y; };
 int iw() { return (int)w; };
 int ih() { return (int)h; };
};

class FXClView;
class ClSeq;
/*
class ClAlign {
 public:
 ClSeq* seq1; //pointers to ClSeq objects sharing this overlap
 ClSeq* seq2;
 FXint l1,r1; // left/right overlap coordinates on seq1
 FXint l2,r2; // left/right overlap coordinates on seq2
 FXdouble pid;
 FXint score;
 FXdouble pvalue;
 ClAlign(ClSeq* nseq1, ClSeq* nseq2,
    FXint nl1, FXint nr1, FXint nl2, FXint nr2,
    FXdouble npid,FXint nscore,FXdouble npvalue) {
     seq1=nseq1;seq2=nseq2; l1=nl1; r1=nr1;
     l2=nl2; r2=nr2; pid=npid; score=nscore; pvalue=npvalue;
     }
 bool operator==(const ClAlign& s) const {
      return pvalue==s.pvalue;
     }
 bool operator<(const ClAlign& s) const {
      return pvalue<s.pvalue;
     }
};*/


class ClSeqSegment {
 public:
  char splL, splR;
  FXint xL,xR; // left/right column coordinates in the layout
  FXint clipL,clipR; // left/right end clipping length
  FXint len() { return xR-xL+1; }
  int seqofs; //offset to the first nucleotide of this segment in the sequence string
  ClSeqSegment(FXint xl, FXint xr, FXint cl=0, FXint cr=0,
					char sl=0, char sr=0, int segpos=0) {
	 xL=xl;xR=xr;
	 clipL=cl;clipR=cr;
	 splL=sl;splR=sr;
	 seqofs=segpos;
     }
 /*
 bool operator==(const ClSeqSegment& s) const {
      return xL==s.xL;
     }
 bool operator<(const ClSeqSegment& s) const {
      return xL<s.xL;
     }
  */
};



class ClSeq {
  public:
   FXClView* view;
   //GList<ClAlign> aligns; //list of pointers to ClAlign objects
   FXString name;
   //FXString sequence; //nt sequence itself!
   char* sequence; //nt sequence could be in fact pulled from the contig sequence
                  // e.g. in the case of GFF files
   FXString comment; //any annotation given
   FXint len; //total length of sequence [in chars = column span],
			  //includes clipping
			  //if segmented: total span of all segments+inter-segment gaps
   FXint slen; //like len, but without intron (sum of segment lengths)
   int numsegs; //numsegs MUST be at least 1
				// >1 means segmented!
   ClSeqSegment* segs; //array of seq segments -- with at least one element!
   FXint cl_xpos; //left position in cluster [in chars] (to be determined)
   FXint cl_ypos; //line order number in cluster (to be determined)
   unsigned int ins, dels; //gaps and deletions for this within the assembly
   FXint clipL;
   FXint clipR;
   short group; //or track_id for GFF files
   int cdsofs; //offset of cds start, relative to the beginning of this ClSeq
   unsigned char reversed;
   unsigned char isET;
   FXint flag;
   ClSeq(FXClView* v, const char* vname, FXint vlen, FXint offs=0,
		   bool minus=false, int l_clip=0, int r_clip=0, const char* seq=NULL,
		   unsigned int ains=0, unsigned int adels=0) {
	 name=vname;
	 len=vlen;
	 //sequence=seq;
	 sequence=Gstrdup(seq); //no need to make a copy?
	 view=v;
	 flag=0;
	 group = -1;
	 cdsofs=0;
	 clipL=l_clip;
	 clipR=r_clip;
	 cl_xpos=offs;
	 cl_ypos=-1;
	 ins=ains;
	 numsegs=0;
	 segs=NULL;
	 addSegment(offs,offs+len-1,clipL, clipR);
	 reversed=minus?1:0;

	 isET= (name.find("np|")>=0 || name.find("et|")>=0 ||
			 name.find("egad|")>=0)?1:0;
	 dels=adels;
	 }
   void addSegment(int xl, int xr, int clipl=0, int clipr=0,
					   char lspl=0, char rspl=0, int seqp=0) {
	 GREALLOC(segs,(numsegs+1)*sizeof(ClSeqSegment));
	 segs[numsegs].xL=xl; segs[numsegs].xR=xr;
	 segs[numsegs].clipL=clipl;
	 segs[numsegs].clipR=clipr;
	 segs[numsegs].splL=lspl;
	 segs[numsegs].splR=rspl;
	 segs[numsegs].seqofs=seqp;
     slen+=xr-xl+1;
	 numsegs++;
	 }
   ~ClSeq() {
	 GFREE(segs);
	 GFREE(sequence);
	 }
   ClSeq(FXClView* v, LytSeqInfo* seqinfo, const char* seq=NULL);
   //FXbool InRect(const CRect& rect, CRect& seqrect);
   //for testing if repainting is needed for this object
   //rect is in pixels, but already offset to the whole cluster area

   //sort by y coordinate
   void calcViewRect(CRect& seqrect); //compute the sequence position
									  //within the view
   void calcSegView(int segidx, CRect& segrect);//compute a seq segment's position
									  //within the view
   bool operator==(const ClSeq& s) const {
	  return (name==s.name);
	 }
   bool operator<(const ClSeq& s) const {
	  return name<s.name;
	 }
};

// utility class for layout compaction:
// a collection of these must always be stored and kept in cl_ypos order
class CLayoutRow {
 public:
  GList<ClSeq> seqs; //sequences shown on this row
  int rpos; //coordinate for the righmost end for any of seqs
            //(where the room starts)
  CLayoutRow(ClSeq* s, int idx): seqs(false,false,false) {
    s->cl_ypos=idx;
    seqs.Add(s);
    rpos=s->cl_xpos + s->len;
    }

  void addSeq(ClSeq* s, int idx) {
    s->cl_ypos=idx;
    seqs.Add(s);
    rpos=s->cl_xpos + s->len;
    }

  char getNuc(int c) {
    for (int i=0;i<seqs.Count(); i++) {
      ClSeq* s=seqs[i];
      if (c<s->cl_xpos) return 0;
      if (c>=s->cl_xpos+s->clipL && c<s->cl_xpos+s->len-s->clipR)
		//return s->sequence.empty() ? ' ' :
        return s->sequence==NULL ? ' ' :
		    toupper(s->sequence[c-s->cl_xpos]);
			// toupper(s->getNuc(c-s->cl_xpos));
      }
    return 0;
    }

 bool operator==(const CLayoutRow& s) const {
      return this==&s;
     }
 bool operator<(const CLayoutRow& s) const {
      return this<&s;
     }
 };

//layout column info:
//this is extended to include SNP candidate detection

struct NTColumn {
 char letter;
 int count;
 };

class ColumnData {
 public:
   char letter; //the consensus letter
   int  thickness; //how many active layers in this column
   int mismatches; //how many mismatches among these
   NTColumn ntdata[5]; //sorted - largest counts first!
   int nN;
   int gap;
  ColumnData() {
   letter=0;
   thickness=0;
   mismatches=0;
   memset(&ntdata,0,5*sizeof(NTColumn));
   nN=0;
   gap=0;
   }
  bool operator==(const ColumnData& s) const {
     return this==&s;
     }
  bool operator<(const ColumnData& s) const {
     return this<&s;
     }
 };

class FXAPI FXClView : public FXScrollArea {
  FXDECLARE(FXClView)
  friend class ClSeq;
protected:
  FXFont      *font;           // Font
  FXFont      *seqfont;           // Font
  FXColor     textColor;      // Text color
  FXImage     *backbuf;
  CScale      scale;
  FXColor*    grpColors[2];
  bool hasSeqs;
  bool canShowSeq;
  int maxThickness;
  // non-scrollable grid area at the top
  FXColor gridColor;
  FXColor gridColorSh; //shadow
  FXColor gridColorH; //highlight
  int content_w;
  int content_h;
  int gridStep;
  int gridH;
  /* -- this section was moved to mainwin.h :
  FXColor baseColors[6]; //base colors
  FXColor seqColors[12];
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

  FXColor matchColors[MISMATCH_CLRIDX+1]; //luminance ramp
  FXColor ctgQuals[NOCTG_CLRIDX+1]; //also a ramp, for contig
  */
  FXColor* baseColors;
  FXColor* seqColors;
  FXColor* matchColors; //luminance ramp for coverage
  FXColor* ctgQuals;  //luminance ramp for contig colors
protected:
  FXClView();
  //virtual void layout();
  void reposVScroll(int yscroll, double r0, int cy);
  void reposHScroll(int xscroll, double c0, int cx);
  bool vscroll;
  bool hscroll;
private:
  FXClView(const FXClView&);
  FXClView &operator=(const FXClView&);
  FXPoint lastMousePos;
public:
  enum {
    ID_TIPTIMER=FXScrollArea::ID_LAST,
    ID_LOOKUPTIMER,
//    ID_TREELIST,
    ID_REGIONCHANGE,
    ID_LAST
    };
  enum ColorStyle {
    csDefault,
    csBaseColor,
    csDensity,
    csMismatch //only if contig is available
    };
protected:
  ColorStyle colorStyle;
public:
  void initColors(FXColor* seqcolors, FXColor* basecolors,
      FXColor* matchcolors, FXColor* ctgquals, FXColor* grpcolors, FXColor* grpgapcolors);
  void scrollBy(int x,int y);
  //void initGrpColors(int num);
  bool setSeqGrp(FXString& seq, int grp);
  void paintGrid(FXDCWindow* dc, FXRectangle& gridR);
  void paintGap(FXDCWindow* dc, FXColor gapcolor, FXRectangle& ClpRct, int y, int x1, int x2,
                int gapl, int gapr, char splL, char splR);
  void paintSeq(FXDCWindow* dc, FXRectangle& clipR, ClSeq* seq,
                        ColorStyle colorstyle);
  void paintSeqBorder(FXDCWindow* dc, FXRectangle& clipR, ClSeq* seq,
                        ColorStyle colorstyle);
  long onConfigure(FXObject*,FXSelector,void*);

  long onPaint(FXObject*,FXSelector,void*);
  long onLeftBtnPress(FXObject*,FXSelector,void*);
  long onLeftBtnRelease(FXObject*,FXSelector,void*);
  long onClicked(FXObject*,FXSelector,void*);
  long onCommand(FXObject*,FXSelector,void*);
  long onMotion(FXObject*,FXSelector,void*);
  long onEnter(FXObject*,FXSelector,void*);
  long onLeave(FXObject*,FXSelector,void*);
  long onQueryTip(FXObject*,FXSelector,void*);
  long onQueryHelp(FXObject*,FXSelector,void*);
  long onTipTimer(FXObject*,FXSelector,void*);
  virtual bool canFocus() const {return true;}

  ClSeq* getSeqAt(int x, int y, int& colno);
  void drawColumn(int colno);
  void drawSeq(ClSeq* seq);
  void makeVisible(ClSeq* seq);
  void selectSeq(ClSeq* seq, bool showIt=false);
  void drawTriangle(FXDCWindow& dc, FXColor color, FXint l,FXint t,FXint r,FXint b);
public:
  /// Constructor
  GList<ClSeq>* seqlist;
  GList<CLayoutRow> rows;
  GList<ColumnData> columns;
  char filetype; //type of the file loaded in the view ('L'=*.lyt, 'A'=*.ace)
  bool showContig; //if true, only shown within the grid, if canShowSeq
  ClSeq* contig; // this will point to the sequence which is actually
                 // the contig (if loaded)
  ClSeq* selSeq; //selected sequence
  ClSeq* btnSeq; //seq which was pressed left btn on
  int selCol;
  int btnCol;
  //default values
  int seqfntW; //seqfont char width
  int seqfntH; //seqfont char height

  int seqH; //sequence bar height (pixels)
  double seqW; //sequence char width (pixels)
  int vSpace;
  FXint         XRight; //maximum right margin (in chars)
                        //given by max(cl_xpos+len)
  FXint         XLeft; //maximum left column (in chars)
                       //minimum offset
  ClSeq* addSeq(const char* name, FXint len, FXint offs, bool reversed=false,
                int trim_l=0, int trim_r=0, const char* sequence=NULL,
                unsigned int ains=0, unsigned int adels=0);
  //supporting segmented sequence addition:
  ClSeq* addSeq(LytSeqInfo* seqinfo, const char* sequence=NULL);
  void addContig(const char* name, FXint len, const char* sequence=NULL, FXint offs=0);
  void buildLayout();
  void getVisibleRegion(FXint& cstart, FXint& cend);
  FXClView(FXComposite *p,FXObject* tgt=NULL,FXSelector sel=0,FXuint opts=0,FXint x=0,FXint y=0,FXint w=0,FXint h=0);
  /// Destructor
  virtual ~FXClView();
   /// Save list to a stream
  virtual void save(FXStream& store) const;
  /// Load list from a stream
  virtual void load(FXStream& store);
  virtual void create();
  virtual void detach();
  //virtual void resize(FXint w,FXint h);
  virtual void position(FXint x, FXint y, FXint w,FXint h);
  virtual void moveContents(FXint x,FXint y);
  //---------------
  /// Change text font
  void setFont(FXFont* fnt);
  /// Return text font
  FXFont* getFont() const { return font; }
  /// Return normal text color
  FXColor getTextColor() const { return textColor; }
  /// Change normal text color
  void setTextColor(FXColor clr);
  //determine when to show scroll bars

  virtual FXint getContentWidth();
  virtual FXint getContentHeight();
  virtual void layout();
  /// Return visible scroll-area width
  virtual FXint getVisibleWidth() const;
  /// Return visible scroll-area height
  virtual FXint getVisibleHeight() const;

  void placeScrollBars(FXint vw, FXint vh);
  void RecalcContent(bool repaint=true);
  FXdouble getZoomX() { return scale.sx; }
  FXdouble getZoomY() { return scale.sy; }
  void ZoomX(FXdouble zx, int fromX=-1);
  void ZoomY(FXdouble zy, int fromY=-1);
  void Zoom(FXdouble zx, FXdouble zy, int fromX=-1, int fromY=-1);
  void setColorStyle(ColorStyle cs) { colorStyle=cs; update(); }
  void Clear() {
    rows.Clear();
    seqlist->Clear();
    //seqaligns->Clear();
    XRight=XLeft=0;
    columns.Clear();
    hasSeqs=false;
    selSeq=NULL;
    selCol=-1;
    RecalcContent();
    }
  //int getXLeft() { return XLeft; }
  bool HasSeqs() { return hasSeqs; }
  int getSelCol() { return selCol; }
  void set1pxZoom(int fromX=-1, int fromY=-1) {
     if (rows.Count()==0 || XRight-XLeft==0) return;
     Zoom(1.00/(double)seqfntW, 1.00/(double)seqfntH,
       fromX, fromY);
     }
  ColumnData* getColumnData(int cidx) { return columns[cidx]; }
  };

#endif
