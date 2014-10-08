#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "GArgs.h"
#include "GStr.h"
#include "GHash.hh"
#include "GList.hh"
#include "GCdbYank.h"
#include "GapAssem.h"
#define USAGE "Usage:\n\
 mblaor <nrcl_layouts_w_gapinfo> -d <fastadb.cidx> [-r <ref_db.cidx]\n\
    [-o <outfile_ace>] [-c <clipmax[%]>] [-p <ref_prefix>]\n\
    \n\
   <nrcl_layouts_w_gapinfo> is the .lyt file with alignment information\n\
               as produced by nrcl\n\
   -d cdb index (created with cdbfasta) for a multi-fasta file containing\n\
      the actual nucleotide sequences of the components (required)\n\
   -r cdb index for a multi-fasta file with the reference sequences\n\
      (only needed if they are not given in the input .lyt file nor\n\
      in <fastadb>)\n\
   -p only consider a nrcl cluster if the reference sequence name starts\n\
      with <ref_prefix>\n\
   -c maximum clipping allowed for component sequences when\n\
      mapped onto a reference; if it ends with '%' then maximum\n\
      clipping is estimated for each read as <clipmax> percent\n\
      of its length\n\
   -o the output acefile is written to <outfile.ace> instead of stdout\n\
   -G do not remove consensus gaps (default is to edit sequences\n\
      in order to remove gap-dominated columns in MSA)\n"

#define LOG_MSG_CLIPMAX "Overlap between %s and reference %s rejected due to clipmax=%4.2f constraint.\n"
#define LOG_MSG_OVLCLIP "Overlap between %s and %s invalidated by the \
%dnt clipping of %s at %d' end.\n"
//-------- global variables 
bool debugMode=false;
bool removeConsGaps=false;
char* ref_prefix=NULL;
bool verbose=false;
int rlineno=0;

static FILE* outf;
static GHash<GASeq> seqs(false);
static GList<GSeqAlign> alns(true, true,false);
                        // sorted, free element, not unique

// flag for reference status:
static const unsigned char flag_IS_REF=0;

//this is only for tree/transitive embedding:
static const unsigned char flag_HAS_PARENT=1;

float clipmax=0;
//--------------------------------
class RefAlign {
  char* linecpy;
  int lineno; //input file line#
  char* gpos;
  char* gpos_onref;
  char* gaps;//char string with gaps in this read
  char* gaps_onref;//char string with gaps in reference
 public:
  char* seqname; //aligned read
  int seqlen; //length of this read
  int offset; //offset of this read relative to reference
  int clip5;
  int clip3;
  char reverse; //is this reversed?
  void parseErr(int fldno);
  RefAlign(char* line, int len, int lno);
  ~RefAlign();
  int nextRefGap(int& pos);
  int nextSeqGap(int& pos);
 };

void loadAlnSeqs(GSeqAlign* aln, GCdbYank* cdbynk, GCdbYank* refcdb=NULL);
void printDebugAln(FILE* f, GSeqAlign* aln, int num, GCdbYank* cdbynk, GCdbYank* refcdb=NULL);

//returns the end of a space delimited token
void skipSp(char*& p) {
 while (*p==' ' || *p=='\t' || *p=='\n') p++;
}
char* endSpToken(char* str) {
 if (str==NULL || *str==0) return NULL;
 char* p=str;
 for (;*p!=0;p++)
  if (*p==' ' || *p=='\t' || *p=='\n') return p;   
 return p; // *p is '\0'
}

//-- prepareMerge checks clipping and even when no clipmax is given,
//   adjusts clipping as appropriately

//========================================================
//====================     main      =====================
//========================================================
int main(int argc, char * const argv[]) {
 //GArgs args(argc, argv, "DGvd:o:c:");
 GArgs args(argc, argv, "DGvd:p:r:o:c:");
 int e;
 if ((e=args.isError())>0)
    GError("%s\nInvalid argument: %s\n", USAGE, argv[e]);
 if (args.getOpt('h')!=NULL) GError("%s\n", USAGE);
 debugMode=(args.getOpt('D')!=NULL);
 removeConsGaps=(args.getOpt('G')==NULL);
 verbose=(args.getOpt('v')!=NULL);
 if (debugMode) verbose=true;
 MSAColumns::removeConsGaps=removeConsGaps;
 GStr infile;
 if (args.startNonOpt()) {
        infile=args.nextNonOpt();
        //GMessage("Given file: %s\n",infile.chars());
        }
 //==
 FILE* inf;
 if (!infile.is_empty()) {
    inf=fopen(infile, "r");
    if (inf==NULL)
       GError("Cannot open input file %s!\n",infile.chars());
    }
  else
   inf=stdin;
  GStr s=args.getOpt('c');
  if (!s.is_empty()) {
      bool ispercent=(s[-1]=='%');
      if (ispercent) s.trimR('%');
      int c=s.asInt();
      if (c<=0) GError("Error: invalid -c <clipmax> (%d) option provided (must be "
                             "a positive integer)!\n",c);
      if (ispercent && c>99)
                GError("Error: invalid percent value (%d) for -c option "
                       " (must be an integer between 1 and 99)!\n",c);
      if (ispercent) {
            clipmax = float(c)/100;
            int maxp=iround(100*clipmax);
            if (verbose) fprintf(stderr,
                  "Percentual max clipping set to %d%%\n",
                                               maxp);
            }
         else {
          clipmax=c;
          if (verbose) fprintf(stderr, "Max clipping set to %d bases\n",
                                               c);
          }

      } //clipmax option
  ref_prefix=args.getOpt('p');
  GStr outfile=args.getOpt('o');
  if (!outfile.is_empty()) {
     outf=fopen(outfile, "w");
     if (outf==NULL)
        GError("Cannot open file %s for writing!\n",outfile.chars());
     }
   else outf=stdout;
  //************************************************


  GStr dbidx=args.getOpt('d');
  if (dbidx.is_empty())
    GError("%sError: a cdb index of a fasta file must be provided!\n",USAGE);
  GCdbYank* cdbyank=new GCdbYank(dbidx.chars());
  GCdbYank* refcdb = NULL;
  dbidx=args.getOpt('r');
  if (!dbidx.is_empty())
    refcdb=new GCdbYank(dbidx.chars());

  GLineReader* linebuf=new GLineReader(inf);
  char* line;
  alns.setSorted(compareOrdnum);
  GASeq* refseq=NULL; // current reference sequence
  int ref_numseqs=0;
  int ref_len=0;
  int ref_lend=0, ref_rend=0;
  int gaplen,gappos;
  bool skipRefContig=false;
  while ((line=linebuf->getLine())!=NULL) {
   RefAlign* aln=NULL; 
   if (line[0]=='>') {
     //establish current reference
     char* ref_name=&line[1];
     char* p=endSpToken(ref_name);
     if (p==NULL) GError("Error parsing reference sequence name at:\n%s\n",line);
     if (ref_prefix!=NULL) {
          if (startsWith(ref_name,ref_prefix)) {
              skipRefContig=false;
              }
            else {
              skipRefContig=true;
              goto NEXT_LINE_NOCHANGE;
              }         
          }
     *p=0;
     if (seqs.Find(ref_name))
        GError("Error (line %d): reference sequence %s is not unique in the input data\n",
                             rlineno+1,ref_name);
     p++;
     if (!parseInt(p,ref_numseqs))
           GError("Error parsing the number of components for reference sequence %s\n",ref_name);
     if (!parseInt(p,ref_lend))
           GError("Error parsing the left end of reference sequence %s\n",ref_name);
     if (!parseInt(p,ref_rend))
           GError("Error parsing the right end of reference sequence %s\n",ref_name);
     ref_len=ref_rend;
     refseq=new GASeq(ref_name, 0, ref_len,
                       ref_lend-1, ref_len-ref_rend, 0);
     refseq->setFlag(flag_IS_REF);
     seqs.Add(ref_name,refseq);
     //may be followed by actual nucleotide sequence of this reference
     while (*p!=0 && *p!='\n') {
       if (*p!='\t' && *p!=' ') refseq->extendSeq(*p);
       p++;
       }
     if (refseq->len>0) {
        refseq->endSeq();
        if (refseq->len!=refseq->seqlen)
          GError("Error: reference sequence %s length mismatch "
                 "(declared %d, found %d)\n",
                  refseq->id, refseq->seqlen, refseq->len);
        }
     goto NEXT_LINE_NOCHANGE;
     } //reference line
   //-- component/read line --
   if (skipRefContig) goto NEXT_LINE_NOCHANGE;
     else {
   //else:
   //-------------------------------------------
   aln=new RefAlign(line, linebuf->length(), rlineno+1);
   if (strcmp(aln->seqname, refseq->id)==0)
     goto NEXT_LINE_NOCHANGE; //skip redundant inclusion of reference as its own child
   if (seqs.Find(aln->seqname))
      GError("Error (line %d): component sequence %s must be unique in the input data\n",
                    rlineno+1,aln->seqname);
   if (clipmax>0) {
     //check if any of the clipping involved is larger than clipmax
     int maxovh=(clipmax<1.00)?
            iround(clipmax * (float)aln->seqlen) :
            (int)clipmax;
     if (aln->clip3>maxovh || aln->clip5>maxovh) {
       if (verbose)
         fprintf(stderr, LOG_MSG_CLIPMAX,
                 aln->seqname, refseq->id, clipmax);
       goto NEXT_LINE_NOCHANGE;
       }// bad mismatching overhangs
     }
   GASeq* aseq=new GASeq(aln->seqname,aln->offset,aln->seqlen,
                           aln->clip5,aln->clip3,aln->reverse);
   GASeq* rseq;
   if (refseq->msa==NULL) {
     rseq=refseq;
     }
   else rseq=new GASeq(refseq->id,0,refseq->seqlen,
                            refseq->clp5,refseq->clp3,0);
   gappos=0;
   while ((gaplen=aln->nextRefGap(gappos))>0)
        rseq->setGap(gappos-1,gaplen);
   gappos=0;
   while ((gaplen=aln->nextSeqGap(gappos))>0)
        aseq->setGap(gappos-1,gaplen);
   //for mgblast alignment, only the query can be reversed
   if (aseq->revcompl==1)
         aseq->reverseGaps(); //don't update offset & reverse flags

   GSeqAlign *newaln=new GSeqAlign(rseq, aseq);
   if (rseq==refseq) {//first alignment with refseq
     newaln->incOrd();
     alns.Add(newaln);
     goto NEXT_LINE;
     }
   refseq->msa->addAlign(refseq,newaln,rseq);
   delete newaln;
   seqs.Add(aseq->name(),aseq);
   }// parsing new alignment
   //---------------------
 NEXT_LINE:
   /* debug print the progressive alignment */
   if (debugMode) {
    for (int a=0;a<alns.Count();a++) {
      printDebugAln(outf,alns.Get(a),a+1,cdbyank, refcdb);
      }
    }
 NEXT_LINE_NOCHANGE:
    //------------
   delete aln;
   if (linebuf->isEof()) break;
   rlineno++;
   //-------------
   if (verbose) {
     if (rlineno%1000==0) {
        fprintf(stderr, "..%d", rlineno);
        if (rlineno%8000==0) fprintf(stderr, "\n");
        fflush(stderr);
        }
     }
  }  //----> the big loop of parsing input lines from the .lyt files
  //*************** now build MSAs and write them
  delete linebuf;
  if (verbose) {
     fprintf(stderr, "\n%d input lines processed.\n",rlineno);
     fprintf(stderr, "Refining and printing %d MSA(s)..\n", alns.Count());
     fflush(stderr);
     }

  //print all the alignments
  for (int i=0;i<alns.Count();i++) {
   GSeqAlign* a=alns.Get(i);
   loadAlnSeqs(a,cdbyank, refcdb); //loading actual sequences for this cluster
   if (debugMode) { //write plain text alignment file
     fprintf(outf,">Alignment%d (%d)\n",i+1, a->Count());
     a->print(outf, 'v');
     }
   else {//write a real ACE file
     //a->buildMSA();
     GStr ctgname;
     ctgname.format("AorContig%d",i+1);
     a->writeACE(outf, ctgname.chars());
     a->freeMSA(); //free MSA and seq memory
     }
   } // for each PMSA cluster
  // oooooooooo D O N E oooooooooooo
  alns.Clear();
  seqs.Clear();
  fflush(outf);
  if (outf!=stdout) fclose(outf);
  if (inf!=stdin) fclose(inf);
  delete cdbyank;
  delete refcdb;
  //GFREE(ref_prefix);
  //GMessage("*** all done ***\n");
  #ifdef __WIN32__
  //getc(stdin);
  #endif
}

//---------------------- RefAlign class
void RefAlign::parseErr(int fldno) {
 fprintf(stderr, "Error parsing input line #%d (field %d):\n%s\n",
         lineno, fldno, linecpy);
 exit(3);
}

RefAlign::RefAlign(char* line, int len, int lno) {
  lineno=lno;
  linecpy=Gstrdup(line);
  seqname=line;
  char* p=endSpToken(seqname);
  if (p==NULL) parseErr(0);
  *p=0;
  p++;
  skipSp(p);
  if (*p!='-' && *p!='+') parseErr(1);
  reverse=0;
  if (*p=='-') reverse=1;
  p++;
  clip5=0;
  clip3=0;
  if (!parseInt(p, seqlen)) parseErr(2);  
  if (!parseInt(p, offset)) parseErr(3);
  offset--;
  skipSp(p);
  if (isdigit(*p)) {//clipping follows
      if (!parseInt(p, clip5)) parseErr(4);
      if (!parseInt(p, clip3)) parseErr(5);
      }
  if (reverse!=0) Gswap(clip5, clip3);
  skipSp(p);
  gaps=NULL;
  gaps_onref=NULL;
  gpos=NULL;
  gpos_onref=NULL;
  if (*p==0) return;
  while (*p!=0 && *p!='R') {
     while (*p!=' ' && *p!='\t' && *p!=0) p++;
     skipSp(p);
     }
  if (*p==0) return;
  // it is R
  p++;
  if (*p!=':') parseErr(6); //invalid field
  p++;
  if (isdigit(*p)) {//should be a number here
     gaps=p;gpos=p;
     }
    else if (*p!='/') parseErr(6);
  while (*p!=0 && *p!=' ' && *p!='/' && *p!='\t') p++;
  if (*p==0) return;
  if (*p=='/') {
      *p=0;//split here
      p++;
      gaps_onref=p;
      gpos_onref=p;
     }
  }
RefAlign::~RefAlign() {
  GFREE(linecpy); //that was just for error messages
  }

int RefAlign::nextSeqGap(int& pos) {
 int r=1;
 if (gpos==NULL || *gpos==0) return 0;
 if (!parseInt(gpos,pos)) parseErr(6);
 if (pos<=0) parseErr(6);
 if (*gpos=='+') { //gap length following
     gpos++;
     if (!parseInt(gpos,r)) parseErr(6);
     if (r<=0) parseErr(6);
     }
 if (*gpos==',') gpos++;
 return r;    
}

int RefAlign::nextRefGap(int& pos) {
 int r=1;
 if (gpos_onref==NULL || *gpos_onref==0) return 0;
 if (!parseInt(gpos_onref,pos)) parseErr(6);
 if (pos<=0) parseErr(6);
 if (*gpos_onref=='+') { //gap length following
     gpos_onref++;
     if (!parseInt(gpos_onref,r)) parseErr(6);
     if (r<=0) parseErr(6);
     }
 if (*gpos_onref==',') gpos_onref++;
 return r;    
}

void printDebugAln(FILE* f, GSeqAlign* aln, int num, GCdbYank* cdbynk, GCdbYank* refcdb) {
  fprintf(f,">[%d]DebugAlign%d (%d)\n",rlineno+1, num, aln->Count());
  loadAlnSeqs(aln,cdbynk, refcdb);
  aln->print(f,'=');
}

void loadAlnSeqs(GSeqAlign* aln, GCdbYank* cdbynk, GCdbYank* refcdb) {
  //FastaSeq seq;
  for (int i=0;i<aln->Count();i++) {
    GASeq* s=aln->Get(i);
    if (s->len==0) {//not loaded yet
      GCdbYank* yankdb=(s->hasFlag(flag_IS_REF) && refcdb!=NULL) ? refcdb : cdbynk;
      if (yankdb->getRecord(s->id, *s) <= 0)
        GError("Error retrieving sequence %s from database %s!\n", 
                                   s->id, yankdb->getDbName());
      if (s->seqlen!=s->len) 
        GError("Error: sequence %s length mismatch! Declared %d, retrieved %d\n",
                           s->id, s->seqlen, s->len);
      s->allupper();
      s->loadProcessing();
      }
    }
 }

