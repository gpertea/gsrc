#include <stdio.h>
#include <ctype.h>
#include <math.h>
//#define ALIGN_COVERAGE_DATA
#include "GArgs.h"
#include "GStr.h"
#include "GHash.hh"
#include "GList.hh"
#include "GCdbYank.h"
#include "GapAssem.h"
#define USAGE "Usage:\n\
 mblasm <mgblast_sortedhits> -d <fastadb.cidx>\n\
   [-c <clipmax>] [-o <outfile.ace>] [-l <dbload.ald>] [-G][-N][-v]\n\
   [-r <restrict_list>] [-x <exclude_list>]\n\n\
   <mgblast_sortedhits> is a mgblast tabulated data w/ gap info\n\
                        sorted by alignment score\n\
   -d multi-fasta database index (cdbfasta) for all the sequences\n\
      (required parameter)\n\n\
 Options:\n\
   -c maximum clipping allowed when forming clusters \n\
      (default: unlimited clipping, merge more sequences);\n\
      if followed by % then maximum clipping is estimated for\n\
      each read as <clipmax> percent of its length\n\
   -o the output acefile is written to <outfile.ace> instead of stdout\n\
   -l write a data file <dbload.ald> with all indels and other MSA\n\
      information for each read, suitable for database loading\n\
   -x ignore any hits involving reads listed in file <exclude_list>\n\
   -r only consider hits between reads listed in file <restrict_list>\n\
   -G do not remove consensus gaps in the ACE output\n\
   -N do not refine clipping at the end of each read in the MSA\n\
   -v verbose mode (report some progress)\n"

// -p a posteriori detection of chimeric reads and reporting
// -D debug mode: print only incremental alignments and exit
// (instead of ACE file)\n"

#define LOG_MSG_CLIPMAX "Overlap between %s and %s rejected due to clipmax=%d constraint.\n"
#define LOG_MSG_OVLCLIP "Overlap between %s and %s invalidated by the \
%dnt clipping of %s at %d' end.\n"
//-------- global variables 
bool debugMode=false;
bool removeConsGaps=false;
bool verbose=false;
int chimeraThreshold=0;
int rlineno=0;
GHash<int> xcludeList;
GHash<int> seqonlyList;

FILE* outf;
bool fltXclude=false;
bool fltRestrict=false;
GHash<GASeq> seqs(false);
GList<GSeqAlign> alns(true, true, false);
                // sorted, free element, not unique

float clipmax=0;
//--------------------------------
class MGPairwise {
  const char* linecpy;
  char* line; //local duplicate
  int lineno;
  char* gpos[2];
 public:
  bool rejected;
  char* seqname[2];
  int ovl[2]; //length of overlap
  int offset[2];
  int  seqlen[2];
  int ovlStart[2];
  int a5[2];//coords of overlap, as given, only always a5<a3
  int a3[2];
  int ovlEnd[2];
  int ovh5[2];
  int ovh3[2];
  int clip5[2];
  int clip3[2];
  int mism5[2]; //mismatching region at 5' end for each seq
  int mism3[2]; //mismatching region at 3' end for each seq
  int ext5[2]; //virtual "extension" at each end of each seq
  int ext3[2]; // within the mini-layout of THIS alignment
  double pid;
  int score1; 
  double score2;
  char reverse[2];
  char* gaps[2];
  void parseErr(int fldno);
  MGPairwise(const char* line_in, int len, int lno);
  ~MGPairwise();
  int nextGap(int seqidx, int& pos);
 };

int readNames(FILE* f, GHash<int>& xhash);
void loadAlnSeqs(GSeqAlign* aln, GCdbYank* cdbynk);
void printDebugAln(FILE* f, GSeqAlign* aln, int num, GCdbYank* cdbynk);

//-- prepareMerge checks clipping and even when no clipmax is given,
//   adjusts clipping as appropriately
bool prepareMerge(GASeq& lnkseq, MGPairwise& mgaln, int i, int j,
                                  GASeq& newseq, GASeq& alnseq);
bool checkAlnMerge(GASeq& h, GASeq& q, MGPairwise& mgaln,
                  int hi,     int qi, GASeq& hpw, AlnClipOps& qclipops);
void prepareAlnMerge(GASeq& h, GASeq& q, AlnClipOps& qclipops);

//========================================================
//====================     main      =====================
//========================================================
int main(int argc, char * const argv[]) {
 GArgs args(argc, argv, "ADGNvd:r:f:x:s:o:c:");
 int e;
 if ((e=args.isError())>0)
    GError("%s\nInvalid argument: %s\n", USAGE, argv[e]);
 if (args.getOpt('h')!=NULL) GError("%s\n", USAGE);
 debugMode=(args.getOpt('D')!=NULL);
 bool rawAlign=(args.getOpt('A')!=NULL);
 removeConsGaps=(args.getOpt('G')==NULL);
 verbose=(args.getOpt('v')!=NULL);
 if (debugMode) verbose=true;
 MSAColumns::removeConsGaps=removeConsGaps;
 MSAColumns::refineClipping=(args.getOpt('N')==NULL);
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

  // exclude hits to those involving reads in a list:
  s=args.getOpt('x');
  FILE* f=NULL;
  if (!s.is_empty()) {
   if ((f=fopen(s, "r"))==NULL)
      GError("Cannot open read exclusion file '%s'!\n", s.chars());
   int c=readNames(f, xcludeList);
   if (verbose) GMessage("Loaded %d sequences to be excluded.\n", c);
   fclose(f);
   fltXclude=(c>0);
   }

  // restrict hits to those between reads in a list:
  s=args.getOpt('r');
  f=NULL;
  if (!s.is_empty()) {
   if ((f=fopen(s, "r"))==NULL)
      GError("Cannot open read restriction file '%s'!\n", s.chars());
   int c=readNames(f, seqonlyList);
   if (verbose) GMessage("Loaded %d sequences to consider exclusively.\n", c);
   fclose(f);
   fltRestrict=(c>0);
   }
   
  GStr outfile=args.getOpt('o');
  if (!outfile.is_empty()) {
     outf=fopen(outfile, "w");
     if (outf==NULL)
        GError("Cannot open file %s for writing!\n",outfile.chars());
     }
   else outf=stdout;
  //************************************************

  s=args.getOpt('f');
  FILE* fltout=NULL;
  if (!s.is_empty()) {
    if ((fltout=fopen(s,"w"))==NULL) 
        GError("Cannot create file %s for writing!\n",s.chars());
    }

  GStr dbidx=args.getOpt('d');
  if (dbidx.is_empty())
    GError("%sError: a cdb index of a fasta file must be provided!\n",USAGE);

  GCdbYank* cdbyank=new GCdbYank(dbidx.chars());

  //TESTING -- start reading and print every alignment found
  GLineReader* linebuf=new GLineReader(inf);
  char* line;
  alns.setSorted(compareOrdnum);
  
  if (verbose) fprintf(stderr, "Processing hits..\n");
  //-==================== HITS PARSING =========================-
  while ((line=linebuf->getLine())!=NULL) {
   //parse fields as needed
   MGPairwise mgpw(line, linebuf->length(), rlineno+1);
   if (mgpw.rejected) continue;
   if (fltout!=NULL) 
     fprintf(fltout, "%s\n",line); 
     
   GASeq* s[2];
   GSeqAlign *pwaln;
   GASeq* lnkseq1;
   GASeq* lnkseq2;
   if (clipmax>0) {
     //check if any of the clipping involved is larger than clipmax
     // like mgblast, be flexible and allow overhangs up to 10% of the overlap size
     //  or 8% of the shorter sequence
     //
     int maxovh;
     if (clipmax<1.00) {
        int minl=GMIN(mgpw.seqlen[0],mgpw.seqlen[1]);
        maxovh = iround(clipmax * (float)minl);
        }
       else maxovh=(int)clipmax;
        
     if (mgpw.mism5[0]>maxovh || mgpw.mism3[0]>maxovh) {
       /*if (verbose) {
         fprintf(stderr, "[pairwise (%d)] ",maxovh);
         fprintf(stderr, LOG_MSG_CLIPMAX,
                 mgpw.seqname[0],mgpw.seqname[1], clipmax);
         }*/
       goto NEXT_LINE_NOCHANGE;
       }// bad mismatching overhangs
     }
   for (int seqidx=0;seqidx<2;seqidx++) {
      s[seqidx]=new GASeq(mgpw.seqname[seqidx],mgpw.offset[seqidx],
                       mgpw.seqlen[seqidx], mgpw.clip5[seqidx],
                       mgpw.clip3[seqidx],mgpw.reverse[seqidx]);
      s[seqidx]->ext5=mgpw.ext5[seqidx];
      s[seqidx]->ext3=mgpw.ext3[seqidx];
      int gaplen;
      int gappos;
      while ((gaplen=mgpw.nextGap(seqidx,gappos))>0) {
        s[seqidx]->setGap(gappos-1,gaplen);
        }
      }
   //for mgblast alignment, only the query can be reversed
   if (s[0]->revcompl==1)
             s[0]->reverseGaps(); //don't update offset & reverse flags

#ifdef ALIGN_COVERAGE_DATA
   pwaln=new GSeqAlign(s[0], mgpw.ovlStart[0]-1, mgpw.ovlEnd[0]-1,
                                  s[1], mgpw.ovlStart[1]-1, mgpw.ovlEnd[1]-1);
#else
   pwaln=new GSeqAlign(s[0], s[1]);
#endif
   //---------------------
   lnkseq1=seqs.Find(s[0]->name());
   lnkseq2=seqs.Find(s[1]->name());

   if (lnkseq1==NULL && lnkseq2==NULL) {
     //brand new sequences, not seen before
     seqs.Add(s[0]->name(),s[0]);
     seqs.Add(s[1]->name(),s[1]);
     pwaln->incOrd();
     alns.Add(pwaln);
     goto NEXT_LINE;
     }

   if (lnkseq1!=NULL && lnkseq2!=NULL) { //both sequences already in MSAs
      if (lnkseq1->msa==lnkseq2->msa) {//already in the same MSA
                                       
        //only transfer the coverage and minimum clipping information:
      #ifdef ALIGN_COVERAGE_DATA
        lnkseq1->addCoverage(s[0]);
        lnkseq2->addCoverage(s[1]);
      #endif
        delete pwaln;
        goto NEXT_LINE_NOCHANGE;
        }
      //both already in different MSA, merge is possible
      // --
      // checking if overlap isn't in a CLIPPED region already
      // a remaining overlap of at least 40 is required!            
      if (mgpw.a5[0]<=lnkseq1->clp5 && mgpw.a3[0]-lnkseq1->clp5<40) {
        /*if (verbose)
          fprintf(stderr, LOG_MSG_OVLCLIP,
            mgpw.seqname[0],mgpw.seqname[1], lnkseq1->clp5, lnkseq1->id, 5);
            */
         delete pwaln;
         goto NEXT_LINE_NOCHANGE;
         }
      int posclp3=lnkseq1->seqlen-lnkseq1->clp3-1;
      if (mgpw.a3[0]>=posclp3 && posclp3-mgpw.a5[0]<40) {
        /*if (verbose)
          fprintf(stderr, LOG_MSG_OVLCLIP,
            mgpw.seqname[0],mgpw.seqname[1], lnkseq1->clp3, lnkseq1->id, 3);*/
         delete pwaln;
         goto NEXT_LINE_NOCHANGE;
         }            
      //---- same verification for lnkseq2's cluster
      if (mgpw.a5[1]<=lnkseq2->clp5 && mgpw.a3[1]-lnkseq2->clp5<40) {
        /*if (verbose)
          fprintf(stderr, LOG_MSG_OVLCLIP,
            mgpw.seqname[0],mgpw.seqname[1], lnkseq2->clp5, lnkseq2->id, 5);*/
         delete pwaln;
         goto NEXT_LINE_NOCHANGE;
         }
      posclp3=lnkseq2->seqlen-lnkseq2->clp3-1;
      if (mgpw.a3[1]>=posclp3 && posclp3-mgpw.a5[1]<40) {
        /*if (verbose)
          fprintf(stderr, LOG_MSG_OVLCLIP,
            mgpw.seqname[0],mgpw.seqname[1], lnkseq2->clp3, lnkseq2->id, 3);*/
         delete pwaln;
         goto NEXT_LINE_NOCHANGE;
         }
       

      AlnClipOps clipops1; //possible clipping of lnkseq1->msa
      AlnClipOps clipops2; //possible clipping of linkseq2->msa
      bool msaSwap=false;//default: lnkseq2->msa swallows lnkseq1->msa;
      AlnClipOps* clipops=&clipops1; //default: 1 is clipped
      int seqidx1=1;
      int seqidx2=0;

      bool can1=checkAlnMerge(*lnkseq1,*lnkseq2, mgpw, 0, 1, *s[0], clipops2);
      bool can2=checkAlnMerge(*lnkseq2,*lnkseq1, mgpw, 1, 0, *s[1], clipops1);
      if (can1) {
        if (can2) { //both merging possible
           if (clipops1.total>clipops2.total) {
             //more clipping would be there if 2 swallows 1
             //merging into lnkseq1->msa is better
             msaSwap=true;
             }
            /*else {
             //more clipping would be there if 1 swallows 2
             //merging into lnkseq2->msa is better
             }*/
           }//both possible
          else {//cannot2, only merging into lnkseq1->msa is possible
           msaSwap=true;
           }
        }
       else //no can1, perhaps can2?
         if (!can2) {
            //no merging possible!
             /*if (verbose) {
               fprintf(stderr, "[%d aln-merge] ",rlineno);
               fprintf(stderr, "Couldn't merge MSAs of %s and %s\n",
                       mgpw.seqname[0],mgpw.seqname[1]);
               }*/
             clipops1.Clear();
             clipops2.Clear();
             delete pwaln;
             goto NEXT_LINE_NOCHANGE;
             }
      if (msaSwap) { //merging into lnkseq1->msa is better
           GASeq* ls=lnkseq1; lnkseq1=lnkseq2;
           lnkseq2=ls;
           seqidx1=0;
           seqidx2=1;
           clipops=&clipops2;
           }
      //lnkseq2's MSA is higher-scoring, merge all the others into it
      prepareAlnMerge(*lnkseq2, *lnkseq1, *clipops);
      clipops1.Clear();
      clipops2.Clear();
      //prepareAlnMerge also took care of clipping adjustments!
      lnkseq2->msa->addAlign(lnkseq2,pwaln,s[seqidx1]);
      delete pwaln;
      GSeqAlign* oaln=lnkseq1->msa;
      s[seqidx2]->clp3=lnkseq1->clp3;
      s[seqidx2]->clp5=lnkseq1->clp5;
      s[seqidx2]->ext3=lnkseq1->ext3;
      s[seqidx2]->ext5=lnkseq1->ext5;
      lnkseq2->msa->addAlign(s[seqidx2],oaln,lnkseq1);
      //replace sequence entry
      seqs.Add(s[seqidx2]->name(),s[seqidx2]);
      //alns.setFreeItem(false);
      alns.Remove(oaln);

      //alns.setFreeItem(true);
      goto NEXT_LINE;
      }// both sequences were already in other alignments
   //----------- only one seq was in a previous MSA
   int seqidx, seqidnew;
   GASeq* lnkseq;
   if (lnkseq1!=NULL) { //seq[0] was seen before
       seqidnew=1;
       seqidx=0;
       lnkseq=lnkseq1;
       }
   else { //lnkseq2!=NULL, seq[1] was seen before
       seqidnew=0;
       seqidx=1;
       lnkseq=lnkseq2;
       }
   //check layout consistency first
   if (!prepareMerge(*lnkseq,mgpw,seqidx, seqidnew, *s[seqidnew], *s[seqidx])) {
       delete pwaln;
       goto NEXT_LINE_NOCHANGE;
       }
   //prepareMerge also took care of clipping adjustments!
   lnkseq->msa->addAlign(lnkseq,pwaln,s[seqidx]);
   delete pwaln;       
   seqs.Add(s[seqidnew]->name(),s[seqidnew]);
   //---------------------
 NEXT_LINE:
   /* debug print the progressive alignment */
   if (debugMode) {
    for (int a=0;a<alns.Count();a++) {
      printDebugAln(outf,alns.Get(a),a+1,cdbyank);
      }
    }
 NEXT_LINE_NOCHANGE:
    //------------
   if (linebuf->isEof()) break;
   rlineno++;
   //-------------
   /*if (verbose) {
     if (rlineno%1000==0) {
        fprintf(stderr, "..%d", rlineno);
        if (rlineno%8000==0) fprintf(stderr, "\n");
        fflush(stderr);
        }
     }*/
   }  //-------- line parsing loop
  delete linebuf;
  if (verbose) {
     fprintf(stderr, "%d input lines processed.\n",rlineno);
     fprintf(stderr, "Refining and printing %d MSA(s)..\n", alns.Count());
     fflush(stderr);
     }

  //print all the alignments
  for (int i=0;i<alns.Count();i++) {
   GSeqAlign* a=alns.Get(i);
   loadAlnSeqs(a,cdbyank);
   if (debugMode || rawAlign) {
     fprintf(outf,">Alignment%d (%d)\n",i+1, a->Count());
     a->print(outf, 'v');
     }
   else {//write actual ACE file
     GStr ctgname;
     ctgname.format("MblContig%d",i+1);
     a->writeACE(outf, ctgname.chars()); 
     //writeACE() also calls buildMSA() and so
     a->freeMSA(); //free MSA and seq memory
     }
   }
  // oooooooooo D O N E oooooooooooo
  alns.Clear();
  seqs.Clear();
  if (fltout!=NULL) fclose(fltout);
  if (fltXclude) xcludeList.Clear();
  if (fltRestrict) seqonlyList.Clear();
  fflush(outf);
  if (outf!=stdout) fclose(outf);
  if (inf!=stdin) fclose(inf);
  delete cdbyank;

  //GMessage("*** all done ***\n");
  #ifdef __WIN32__
  //getc(stdin);
  #endif
}


//---------------------- MGPairwise class
void MGPairwise::parseErr(int fldno) {
 fprintf(stderr, "Error parsing input line #%d (field %d):\n%s\n",
         lineno, fldno, linecpy);
 exit(3);
}
MGPairwise::MGPairwise(const char* line_in, int len, int lno) {
  char* fields[16];
  int fldno=0;
  int fldstart=0;
  lineno=lno;
  fields[12]=NULL;
  fields[13]=NULL;
  line=Gstrdup(line_in);
  linecpy=line_in;
  
  for (int i=0;i<len;i++) {
     if (line[i]=='\t') {
       line[i]=0;
       fields[fldno]=line+fldstart;
       fldno++;fldstart=i+1;
       if (fldno>15) {
            fprintf(stderr,
              "Warning: too many tab delimited fields found in line %d (%s)\n",
              lineno, linecpy);
            break;
            }
       }//if tab
     }
   fields[fldno]=line+fldstart;
   rejected=true;
   //check input filters:
   if (fltXclude && (xcludeList.hasKey(fields[0]) ||
                     xcludeList.hasKey(fields[4]))) return;
                     
   if (fltRestrict && !(seqonlyList.hasKey(fields[0]) &&
                        seqonlyList.hasKey(fields[4]))) return;
   rejected=false;
   seqname[0]=fields[0];
   if (!parseInt(fields[1], seqlen[0]))   parseErr(1);
   if (!parseInt(fields[2], ovlStart[0])) parseErr(2);
   if (!parseInt(fields[3], ovlEnd[0]))   parseErr(3);
   seqname[1]=fields[4];
   if (!parseInt(fields[5], seqlen[1]))   parseErr(5);
   if (!parseInt(fields[6], ovlStart[1])) parseErr(6);
   if (!parseInt(fields[7], ovlEnd[1]))   parseErr(7);

   if (!parseNumber(fields[8], pid))   parseErr(8);
   if (!parseInt(fields[9], score1))   parseErr(9);
   if (!parseNumber(fields[10], score2))   parseErr(10);
   char c=(fields[11][0]);
   if (c!='-' && c!='+') parseErr(11);
   ovl[1]=ovlEnd[1]-ovlStart[1]+1;
   a5[1]=ovlStart[1];
   a3[1]=ovlEnd[1];
   if (ovlStart[0]>ovlEnd[0]) { // reversed alignment
     a5[0]=ovlEnd[0];
     a3[0]=ovlStart[0];
     ovl[0]=ovlStart[0]-ovlEnd[0]+1;
     int ostart=seqlen[0]-ovlStart[0]+1;
     int oend=ostart+ovlStart[0]-ovlEnd[0];
     ovlStart[0]=ostart;
     ovlEnd[0]=oend;
     }
    else {
      a5[0]=ovlStart[0];
      a3[0]=ovlEnd[0];
      ovl[0]=ovlEnd[0]-ovlStart[0]+1;
      }
   if (ovlStart[0]>=ovlStart[1]) {
      offset[0]=0;
      offset[1]=ovlStart[0]-ovlStart[1];
      }
     else {
      offset[0]=ovlStart[1]-ovlStart[0];
      offset[1]=0;
      }
   reverse[1]=0;
   ovh5[1]= ovlStart[1]-1;
   clip5[1]=ovh5[1];
   ovh3[1]=seqlen[1]-ovlEnd[1];
   clip3[1]=ovh3[1];
   //--only query can be reversed!
   ext3[0]=0;
   ext5[0]=0;
   ext3[1]=0;
   ext5[1]=0;
   //always clip the shorter mismatched end
   if (c=='-') { //reverse
     reverse[0]=1;
     ovh3[0]=ovlStart[0]-1;
     clip3[0]=ovh3[0];
     ovh5[0]=seqlen[0]-ovlEnd[0];
     clip5[0]=ovh5[0];
     /**/
     int d=clip3[0]-clip5[1];
     if (d>0) { //clip5[1] is smaller!
       mism5[1]=clip5[1];
       mism3[0]=clip5[1];
       ext5[1]=d;
       }
      else { //clip3[0] is smaller
       mism3[0]=clip3[0];
       mism5[1]=clip3[0];
       ext3[0]=-d;
       }
     d=clip5[0]-clip3[1];
     if (d>0) {
       mism3[1]=clip3[1];
       mism5[0]=clip3[1];
       ext3[1]=d;
       }
      else {
       mism5[0]=clip5[0];
       mism3[1]=clip5[0];
       ext5[0]=-d;
       }

     //no clipping if no "mismatched" end on the other sequence
     if (clip5[1]==0) clip3[0]=0;
     if (clip3[1]==0) clip5[0]=0;
     /**/
     // no clipping unless mismatch!
     /*if (clip5[0]==0) clip3[1]=0;
     if (clip3[0]==0) clip5[1]=0;*/

     //-- Don't clip the longer end for now
     if (clip3[0]>clip5[1]) clip3[0]=0;
                     else clip5[1]=0;
     if (clip5[0]>clip3[1]) clip5[0]=0;
                     else clip3[1]=0;
     }
    else { //forward
     reverse[0]=0;
     ovh5[0]=ovlStart[0]-1;
     clip5[0]=ovh5[0];
     ovh3[0]=seqlen[0]-ovlEnd[0];
     clip3[0]=ovh3[0];
     /**/
     int d=clip5[0]-clip5[1];
     if (d>0) { //clip5[1] is smaller
       mism5[1]=clip5[1];
       mism5[0]=clip5[1];
       ext5[1]=d;
       }
      else {
       mism5[0]=clip5[0];
       mism5[1]=clip5[0];
       ext5[0]=-d;
       }
     d=clip3[0]-clip3[1];
     if (d>0) { //clip3[1] is smaller
       mism3[1]=clip3[1];
       mism3[0]=clip3[1];
       ext3[1]=d;
       }
      else {
       mism3[0]=clip3[0];
       mism3[1]=clip3[0];
       ext3[0]=-d;
       }
     /*mismL=GMIN(clip5[1],clip5[0]);
     mismR=GMIN(clip3[1],clip3[0]);*/

     //don't clip unless mismatch
     if (clip5[1]==0) clip5[0]=0;
     if (clip3[1]==0) clip3[0]=0;
     // no clipping unless mismatch
     /*if (clip5[0]==0) clip5[1]=0;
     if (clip3[0]==0) clip3[1]=0;*/

     //--don't clip the longer end!
     if (clip5[0]>clip5[1]) clip5[0]=0;
                     else clip5[1]=0;
     if (clip3[0]>clip3[1]) clip3[0]=0;
                     else clip3[1]=0;
     }
   gaps[0]=fields[12];
   gaps[1]=fields[13];
   gpos[0]=gaps[0];
   gpos[1]=gaps[1];
   }

MGPairwise::~MGPairwise() {
  GFREE(line); //that was just for error messages
}

int MGPairwise::nextGap(int seqidx, int& pos) {
 int r=1;
 if (gpos[seqidx]==NULL || *gpos[seqidx]==0) return 0;
 if (!parseInt(gpos[seqidx],pos)) parseErr(12+seqidx);
 if (pos<=0) parseErr(12+seqidx);
 if (*gpos[seqidx]=='+') { //gap length following
     gpos[seqidx]++;
     if (!parseInt(gpos[seqidx],r)) parseErr(12+seqidx);
     if (r<=0) parseErr(12+seqidx);
     }
 if (*gpos[seqidx]==',') gpos[seqidx]++;
 return r;
}

void printDebugAln(FILE* f, GSeqAlign* aln, int num, GCdbYank* cdbynk) {
  fprintf(f,">[%d]DebugAlign%d (%d)\n",rlineno+1, num, aln->Count());
  loadAlnSeqs(aln,cdbynk);
  aln->print(f,'=');
}

void loadAlnSeqs(GSeqAlign* aln, GCdbYank* cdbynk) {
  //FastaSeq seq;
  for (int i=0;i<aln->Count();i++) {
    GASeq* s=aln->Get(i);
    if (s->len==0) {//not loaded yet
      if (cdbynk->getRecord(s->id, *s)<=0 || s->len!=s->seqlen)
              GError("Error retrieving sequence %s from database %s!\n", 
                                   s->id, cdbynk->getDbName());
      s->allupper();
      s->loadProcessing();
      }
    }
 }

// -- check for merging a new pairwise alignment into an existing msa
// clipmax checking but more importantly: clipping adjustment
bool prepareMerge(GASeq& lnkseq, MGPairwise& mgaln, int i, int j, GASeq& newseq, GASeq& alnseq) {
    //                   node h                       h        q     node q (j)          i
  GSeqAlign& msa=*lnkseq.msa;
  int d5, d3;
  int misovh5=mgaln.ovh5[j];
  int misovh3=mgaln.ovh3[j];
  //virtual mismatches for newseq (j)
  char q_rev=0;
  if (mgaln.reverse[0]!=0) { //reverse complement alignment
    d5= misovh5 - (lnkseq.seqlen+lnkseq.ext3-mgaln.a3[i]); //offset of q relative to h at 5' end of q
    d3=misovh3 - (lnkseq.ext5+mgaln.a5[i]-1); //offset of q relative to h at 3' end of q
    if (d5>0) misovh5 -=d5;//GMIN(q5-1, lnkseq.seqlen+lnkseq.ext3-mgaln.a3[i]);
    if (d3>0) misovh3 -=d3;//GMIN(newseq.seqlen-q3, lnkseq.ext5+mgaln.a5[i]-1);
    if (lnkseq.revcompl==0) q_rev=1;
    } // was reverse
 else {  //forward alignment
    d5= misovh5 - (mgaln.a5[i]+lnkseq.ext5-1); //offset of q relative to h at 5' end of q
    d3= misovh3 - (lnkseq.seqlen+lnkseq.ext3-mgaln.a3[i]); //offset of q relative to h at 3' end of q
    if (d5>0) misovh5 -=d5;
    if (d3>0) misovh3 -=d3;
    if (lnkseq.revcompl!=0) q_rev=1;
    }
 //if (hvcov<minscov || misovh5>maxovhang || misovh3>maxovhang)
 if (clipmax>1) {
       if (misovh5>clipmax || misovh3>clipmax) {
         /*if (verbose)
                 fprintf(stderr, LOG_MSG_CLIPMAX,
                       mgaln.seqname[0],mgaln.seqname[1], clipmax);*/
         return false;
         }
    }
 // checking if overlap isn't in a CLIPPED region already
 // a remaining overlap of at least 40 is required!
 if (mgaln.a5[i]<=lnkseq.clp5 && mgaln.a3[i]-lnkseq.clp5<40) {
   /*if (verbose)
     fprintf(stderr, LOG_MSG_OVLCLIP,
       mgaln.seqname[0],mgaln.seqname[1], lnkseq.clp5, lnkseq.id, 5);*/
   return false;
   }
 int posclp3=lnkseq.seqlen-lnkseq.clp3-1;
 if (mgaln.a3[i]>=posclp3 && posclp3-mgaln.a5[i]<40) {
   /*if (verbose)
     fprintf(stderr, LOG_MSG_OVLCLIP,
       mgaln.seqname[0],mgaln.seqname[1], lnkseq.clp3, lnkseq.id, 3);*/
   return false;
   }


 // it passed the virtual overhang tests, now adjusts ext5 and ext3
 // and clipping

 int alnclp5=-1;
 int alnclp3=-1; //for msa->evalClipping(&alnseq,alnclp5, alnclp3);
 AlnClipOps clipops;
 bool canClip5=true;
 bool canClip3=true;
 if (mgaln.reverse[0]!=0) {
   if (misovh5>0 || lnkseq.clp3>0) {
       alnclp3=GMAX(lnkseq.clp3, mgaln.mism3[i]);
       canClip3=alnseq.msa->evalClipping(&alnseq,alnclp5,alnclp3, clipmax,
                                  clipops);
       }

   if (misovh3>0 || lnkseq.clp5>0) {
       alnclp5=GMAX(lnkseq.clp5, mgaln.mism5[i]);
       canClip5=alnseq.msa->evalClipping(&alnseq,alnclp5,alnclp3, clipmax,clipops);
       }
   }
 else {
   if (misovh5>0 || lnkseq.clp5>0) {
      alnclp5=GMAX(lnkseq.clp5, mgaln.mism5[i]);
      canClip5=alnseq.msa->evalClipping(&alnseq,alnclp5,alnclp3,clipmax, clipops);
      }
   if (misovh3>0 || lnkseq.clp3>0) {
      alnclp3=GMAX(lnkseq.clp3, mgaln.mism3[i]);
      canClip3=alnseq.msa->evalClipping(&alnseq,alnclp5,alnclp3,clipmax, clipops);
      }
  }

 if (canClip5 && canClip3) alnseq.msa->applyClipping(clipops);
                      else return false;
 //-- adjust virtual extensions ext3 and ext5
 newseq.ext3=-d3;
 newseq.ext5=-d5;
 //--
 //adjust end-ext for all sequences in the old msa
 if (d5>0) {
     newseq.ext5+=d5;
     /*if (newseq.revcompl==q_rev) newseq.ext5+=d5;
                            else newseq.ext3+=d5;*/
     for (int i=0;i<msa.Count();i++) {
       GASeq* n=msa.Get(i);
       if (n->revcompl==q_rev) n->ext5+=d5;
                          else n->ext3+=d5;
       }
     }
 if (d3>0) {
     newseq.ext3+=d3;
     /*if (newseq.revcompl==q_rev) newseq.ext3+=d3;
                          else newseq.ext5+=d3;*/
     for (int i=0;i<msa.Count();i++) {
        GASeq* n=msa.Get(i);
        if (n->revcompl==q_rev) n->ext3+=d3;
                        else n->ext5+=d3;
       }
   }
 return true;
}

bool checkAlnMerge(GASeq& h, GASeq& q, MGPairwise& mgaln,
                  int hi,     int qi, GASeq& hpw, AlnClipOps& qclipops)       {
 /* checking if h (hmsa) can swallow q (qmsa)
    -- by trimming q and qmsa appropriately */
 //GSeqAlign& hmsa=*h.msa; //the "better" MSA ?
 GSeqAlign& qmsa=*q.msa;
 //int hvcov = (int) iround(((double)(h3-h5)*100)/(hlen+h->ovx3+h->ovx5));
 //int qvcov = (int) iround(((double)(q3-q5)*100)/(qlen+q->ovx3+q->ovx5));
 //int misovh5, misovh3, //mismatched overhangs at 5' and 3' of q
 int misovh5=q.ext5+mgaln.ovh5[qi];
 int misovh3=q.ext3+mgaln.ovh3[qi];
 qclipops.q_rev=0; //indicates if one of the clusters should be reversed
 if (mgaln.reverse[0]!=0) { //reverse complement alignment
   qclipops.d5 = misovh5 - (mgaln.seqlen[hi]+h.ext3-mgaln.a3[hi]); //offset of q's cluster relative to
                                  //h's cluster at 5' end of q
   qclipops.d3 = misovh3 - (h.ext5+mgaln.a5[hi]-1); //offset of q's cluster relative to
                                  //h's cluster at 3' end of q
   if (qclipops.d5>0) misovh5 -=qclipops.d5;//GMIN(q5-1, mgaln.seqlen[hi]+h.ext3-mgaln.a3[hi]);
   if (qclipops.d3>0) misovh3 -=qclipops.d3;//GMIN(qlen-q3, h.ext5+mgaln.a5[hi]-1);
   if (h.revcompl==q.revcompl) qclipops.q_rev=1;
   }
 else { //forward strand alignment
   qclipops.d5= misovh5 - (mgaln.a5[hi]+h.ext5-1); //offset of q->cluster relative to h->cluster at 5' end of q
   qclipops.d3= misovh3 - (mgaln.seqlen[hi]+h.ext3-mgaln.a3[hi]); //offset of q->cluster relative to h at 3' end of q
   if (qclipops.d5>0) misovh5 -=qclipops.d5;
   if (qclipops.d3>0) misovh3 -=qclipops.d3;
   if (h.revcompl!=q.revcompl) qclipops.q_rev=1;
   }
 //--it passed the virtual overhang tests
 //if (hvcov<minscov || qvcov<minscov || misovh5>maxovhang || misovh3>maxovhang)
 if (clipmax>1) {
       if (misovh5>clipmax || misovh3>clipmax) {
         /*if (verbose)
                 fprintf(stderr, LOG_MSG_CLIPMAX,
                       mgaln.seqname[0],mgaln.seqname[1], clipmax);*/
         return false;
         }
    }

 //--- adjust clipping -- keep in mind that h is to be taken as
 //    the "better" sequence, so we're only supposed to trim q and qmsa
 // If h is clipped in hmsa, there is possible EXTRA clipping of q/qmsa
 // inherited from h through the pairwise alignment hpw!
 int hclp5=-1;
 int hclp3=-1;
 if (h.clp3>0) hclp3=h.clp3;
 if (h.clp5>0) hclp5=h.clp5;
 //clip the pairwise alignment by h clipping in hmsa
 //
 // we should only do this if there is any intersection
 // between the clipped ends of h and qmsa through h->q pairwise alignment?
 if (!hpw.msa->evalClipping(&hpw,hclp5,hclp3,0,qclipops))
            return false;
 int qovh3=(misovh3>0) ? mgaln.ovh3[qi]: mgaln.mism3[qi];
 int qovh5=(misovh5>0) ? mgaln.ovh5[qi]: mgaln.mism5[qi];
 if (qclipops.Count()>0) {
    //hpw.msa->applyClipping(qclipops);
    for (int c=0;c<qclipops.Count();c++) {
      SeqClipOp& op=*qclipops.Get(c);
      if (op.seq==&hpw) continue; //shouldn't be the case
      if (qovh5<op.clp[0]) qovh5=op.clp[0];
      if (qovh3<op.clp[1]) qovh3=op.clp[1];
      }
    qclipops.Clear();//but it keeps the total!
    }
 /*int alnclp5=-1;
 int alnclp3=-1;*/
 bool canClip5=true;
 bool canClip3=true;
 int qmax3=GMAX(q.clp3, qovh3);
 int qmax5=GMAX(q.clp5, qovh5);

 //-- just to avoid some special cases when q is clipped
 //   but q's msa does not exhibit any "mismatch" with h's msa:
 // **TODO: this should only be the case if there is any overlap between
 // **the two MSAs which would INTERSECT any of these clipped regions of q!
 if (q.ext5+qmax5 > qclipops.d5) {
  if (misovh5==0) misovh5=q.clp5;
  if (misovh5>0 || qmax5>q.clp5) //{
     //alnclp5=qmax5;
     canClip5=qmsa.evalClipping(&q, qmax5,-1,clipmax,qclipops);
     //}
  }
 if (q.ext3+qmax3 > qclipops.d3) {
   if (misovh3==0) misovh3=q.clp3;
   if (misovh3>0 || qmax3>q.clp3) //{
     //alnclp3=qmax3;
     canClip3=qmsa.evalClipping(&q,-1,qmax3,clipmax,qclipops);
     //this will adjust newseq clipping also!
     //}
  }
return (canClip5 && canClip3);
}

void prepareAlnMerge(GASeq& h, GASeq& q, AlnClipOps& qclipops) {
//h (hmsa) is to swallow q (qmsa)
GSeqAlign& hmsa=*h.msa; //the "better" MSA
GSeqAlign& qmsa=*q.msa;
qmsa.applyClipping(qclipops);
//char prevqrev=q->nrev; //current state of q
for (int i=0;i<qmsa.Count();i++) {
 GASeq* qn=qmsa.Get(i);
 if (qn->revcompl==q.revcompl) { //same orientation with q
    qn->ext3+=-qclipops.d3;
    qn->ext5+=-qclipops.d5;
    }
   else { //rev.compl. relative to q
    qn->ext5+=-qclipops.d3;
    qn->ext3+=-qclipops.d5;
    }
 //if (q_rev) qn->nrev=!qn->nrev; -- IMPORTANT FOR FURTHER q considerations!
 // -- also adjust ext5/ext3 for this cluster
 }
//char newqrev = (q_rev==1) ? (q.revcompl^1) : q.revcompl;
char newqrev = qclipops.q_rev ^ q.revcompl; //XOR
//-- adjust ovx (ext)
 if (qclipops.d5>0) {
   //adjust query cluster ext:
   for (int i=0;i<qmsa.Count();i++) {
     GASeq* qn=qmsa.Get(i);
     if ((qclipops.q_rev ^ qn->revcompl)==newqrev) qn->ext5+=qclipops.d5;
                                              else qn->ext3+=qclipops.d5;
      }
   for (int i=0;i<hmsa.Count();i++) {
     GASeq* hn=hmsa.Get(i);
     if (hn->revcompl==newqrev) hn->ext5+=qclipops.d5;
                           else hn->ext3+=qclipops.d5;
     }
  }//d5>0
 if (qclipops.d3>0) {
    for (int i=0;i<qmsa.Count();i++) {
       GASeq* qn=qmsa.Get(i);
       if ((qclipops.q_rev ^ qn->revcompl)==newqrev) qn->ext3+=qclipops.d3;
                                                else qn->ext5+=qclipops.d3;
        }
    for (int i=0;i<hmsa.Count();i++) {
      GASeq* n=hmsa.Get(i);
      if (n->revcompl==newqrev) n->ext3+=qclipops.d3;
                           else n->ext5+=qclipops.d3;
      }
  }
}

int readNames(FILE* f, GHash<int>& xhash) {
  int c;
  int count=0;
  char name[256]; 
  int len=0;
  while ((c=getc(f))!=EOF) {
    if (isspace(c)) {
      if (len>0) {
        name[len]='\0';
        xhash.Add(name, NULL);
        count++;
        len=0;
        }
      continue;  
      }
    //a non-space
    name[len]=(char) c;
    len++;
    if (len>255) {
      name[len-1]='\0';
      GError("Error reading file: name too long ('%s') !\n",name);      
      }
    }
 return count;   
}

