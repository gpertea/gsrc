#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "GBase.h"
#include "GArgs.h"
#include "GStr.h"
#include "GHash.hh"
#include "GList.hh"

#define usage "Perform containment clustering by filtering tabulated hits.\n\
The 'representative' (longest) sequences are listed first within each cluster.\n\
\n\
Usage:\n\
 nrcl [<hits_file>] [-H] [-t] [-n] [-S] [-o <out_file>] [-y <layouts_file>]\n\
 [-L] [-l <parent_list>] [-p <parent_prefix] [SEQFLT={ET|EST|ET2EST}]\n\
 [-c <clsfile>] [-x <excludefile>] [-f <flthits_file>] [-d <debug_tree>]\n\
 [SCOV=xx] [OVHANG=xx] [SCORE=xx] [PID=xx]\n\
 Options:\n\
 -n  : non-unique sequence representation across clusters; allow clusters \n\
       to share sequences\n\
 -t  : assume containment transitivity (tree structure)\n\
 -o  : write output in <out_file> instead of stdout\n\
 -y  : write cluster layouts to file <layouts_file>\n\
 -f  : write to <flthits_file> all the lines that passed the filters\n\
 -S  : write the 'singletons' to the result file too\n\
 -d  : write the containment tree into file <debug_tree>\n\
 -L  : sequence length, instead of overlap score, is used for\n\
       parenting decisions\n\
 -w  : file with names of parent names and their relative weights\n\
       to be used for parenting decisions (supercedes -L or scoring)\n\
 -p  : sequences whose names start with <parent_prefix> (default: \"np|\")\n\
       are given parenting  priority over those without such prefix\n\
 -l  : sequences whose names are found in file <parent_list> are given\n\
       parenting priority over other sequences\n\
Filters:\n\
     -c   load cluster file <clsfile> and only consider containment within\n\
          the same cluster. This option causes the creation of two extra files:\n\
          a containment layout file <clsfile>.lyt and a container-formatted \n\
          cluster file (one container per line) <clsfile>.scls\n\
     -x   discard any pairs/hits related to sequences listed in <excludefile>\n\
SEQFLT=   filter sequence names for specific sequence type clustering: \n\
          EST  - EST-only clustering (use only ESTvsEST hits) \n\
          ET   - ETs/NPs only clustering\n\
          ET2EST - ignore EST-EST pairs, uses only EST-ET or ET-ET pairs\n\
SCOV=xx   is the minimum coverage of the shorter sequence (percentage,\n\
          default value 90)(overlap as percentage of the shorter sequence)\n\
OVHANG=xx is the maximum overhang (in nucleotides) to be allowed for the \n\
          shorter sequence only (default:unlimited, use SCOV value)\n\
PID=xx    is the minimum percent of identity\n\
SCORE=xx  is the minimum score for a match to be considered\n\
\n\
If <hits_file> is not given, hits are expected at stdin\n\
Each input line must have these tabulated fields:\n\n\
q_name q_len q_n5 q_n3 hit_name hit_len hit_n5 hit_n3 pid score xscore strand\n\
 [q_gapinfo hit_gapinfo]\n"

#define ERR_INVALID_PAIR "Invalid input line encountered:\n%s\n"
#define ERR_CL_PARSE "Error parsing cluster file at line:\n%s\n"


bool by_length=false;
bool wclusters=false;
bool npWeights=false;
bool namesLoaded=false;
int inbuf_len=1024; //starting inbuf capacity
char* inbuf=NULL; // incoming buffer for sequence lines.
const char* prefix="np|";

int base_ordnum=0; //node creation order

FILE* fscls=NULL;

int readClusters(FILE* f);

class GCluster;

class CNode {
 public:
   char* id; //actual node ID
   bool isFull;
   int pweight;
   int ordnum; //creation order#
   int len; //sequence length
   char* gaps;  //list of gaps, as given by mgblast -D5
   char* pgaps; //list of gaps on parent, ..
   unsigned char lytminus;
   int lytpos;
   //=== if it's a master sequence:
   int numseqs; //final number of sequences "swallowed" by this sequence
                //(across the whole subtree, if transitivity is allowed)
   //int l_start; // the leftmost coordinate of an overlap with a child
   //int r_end;   // the rightmost coordinate of an overlap with a child
   //=== if its a child sequence:  
   int poffs; //offset relative to left end of the parent sequence
   short pclipL, pclipR; // overhangs for the overlap of this with the parent
   unsigned char pminus; //1 if complement overlap to parent sequence
   //CNode* longest; //longest sequence found to include this one
   CNode* parent; //best parent sequence for this sequence (best scoring overlap)
   //====
   //the relation is validated at the end only if parent is the same with cluster->master
   int pscore; //current parent overlap score 
   //
   //unsigned char seqtype;
   GCluster* kids; //list of all child nodes (seqs) found to be contained into this
   GCluster* cluster; //for the wcluster option
   GCluster* contained; //just at the end this will collect all contained sequences
   CNode(char* name, int slen, int ordn=-1);

   /*-- node creation when loading from cluster file --*/
   CNode(char* name, GCluster*& cl, int ordn);

   //directly set the parent as specified, with no other checking
   void setParent(CNode* np, int offs, int clpL, int clpR,
                    int score, bool minus, char* g=NULL, char* pg=NULL);

   //check if a sequence np is to be assigned as "parent" to this
   //returns true if parenting works according to the sequence AND/OR link weights
   //or false if the parenting is reversed!
   bool assignParent(CNode* np, int score, int offs, int clpL, int clpR,
                                  bool minus, char* g=NULL, char* pg=NULL);
   ~CNode();

  //-- container comparison operators: just pointer based
    bool operator==(CNode& d){
     return (this==&d);
     }
    bool operator>(CNode& d){
     return (this>&d);
     }
    bool operator<(CNode& d){
     return (this<&d);
     }
 };

class GCluster :public GList<CNode> {
 public:
  int total;
  bool operator==(GCluster& d){
     return (this==&d);
     }
  bool operator>(GCluster& d){
     return (this>&d);
     }
  bool operator<(GCluster& d){
     return (this<&d);
     }
  GCluster():GList<CNode>(true,false,false) {
    //default is: sorted by Node pointer, don't free nodes, unique
    //size=0;
    total=0;
    }
  GCluster(bool sorted, bool free_elements=false, bool beUnique=false):
     GList<CNode>(sorted,free_elements,beUnique) {
    //size=0;
    total=0;
    }
  void addNode(CNode* c) {
    Add(c);
    //size++;
    total++;
    }
  };

// when use clusters, keep them here:
GList<GCluster> tclusters(false, false,false);
                     //not sorted, not free, non-unique


bool transitive=false;
bool asfasta=true;
bool wSinglets=false; //write singletons too
bool do_share=true; //allow clusters sharing sequences

GHash<int> excludeList;
GHash<int> fullGenes;

//-- global hash containing all the nodes with structured links info
// holds all the nodes
//if cluster-targeted containment, 
GHash<CNode> all;

GList<CNode> seqs(false,false,false);

#ifdef DBGTRACE
   GShMem mem("nrcl", 4096);
   GShMem memsum("nrclsum", 1024);
   GStr slog;
#endif

//keeps info about an overlap
class COverlap {
 public:
 char* id1;
 char* id2;
 CNode* n1;
 CNode* n2;
 int len1;
 int len2;
 int o1start; //overlap start coordinate (leftmost) on n1
 int o2start; //overlap start coordinate (leftmost) on n2
 int ovl1;
 int ovl2;
 bool minus;
 int score;
 char* gaps1;
 char* gaps2;
 int pid;
 void Swap() {
  Gswap(id1, id2);
  CNode* sw=n1;
  n1=n2;n2=sw;
  //swap(n1,  n2);
  Gswap(len1, len2);
  Gswap(o1start, o2start);
  Gswap(ovl1,ovl2);
  Gswap(gaps1, gaps2);
  }
  
 COverlap(char* seq2, char* seq1, int l2, int l1, int os2, int os1, 
          int o2, int o1, int sc, int p_id, bool m, char* g2, char* g1) {
   id1=seq1;id2=seq2;len1=l1;len2=l2;
   o1start=os1;
   o2start=os2;
   ovl1=o1;ovl2=o2;score=sc; pid=p_id;
   gaps1=g1;
   gaps2=g2;
   n1=NULL;n2=NULL;
   minus=m;
   }
 };
 
bool flt_ET_only=false;
bool flt_EST_only=false;
bool flt_EST2ET=false;
bool flt_Exclude=false; //if an exclusion list was given
bool ETcheck=false;

int minpid=0, minscov=90, maxovhang = INT_MAX, minscore=0;


FILE* outf=NULL;

FILE* fTrees=NULL;



//=========== functions used:

CNode* addNode(char* id, int len, int ordn=-1);
void addPair(COverlap& o);
CNode* mergeLinks(CNode* n1, CNode* n2);
int getNextValue(char*& p, char* line);
double getNextNumber(char*& p, char* line);
//bool seq_filter(char* s1, char* s2); //returns true if the pair can pass through the filter
int readNames(FILE* f, GHash<int>& xhash);
int readWeights(FILE* f, GHash<int>& xhash);
void collectNodes(GCluster* C, CNode* n, bool dbg=false);

//============ comparison functions:
int comparePtr(void* p1, void* p2);
//compare sequences by length
int compareLen(void* p1, void* p2);
//compare sequences by num of contained
int compareContained(void* p1, void* p2);
//=====
//compare cluster sizes:
int compareSize(void* p1, void* p2);

//=====
int compareID(void* p1, void* p2);


//========================================================
//====================     main      =====================
//========================================================
int main(int argc, char* const argv[]) {
 GMessage("Command line was:\n ");
 for (int i=0;i<argc-1;i++) 
   GMessage("%s ", argv[i]);
 GMessage("%s\n", argv[argc-1]);
 GArgs args(argc, argv, "hntSFLl:p:c:f:o:y:x:d:w:SEQFLT=PID=SCOV=OVHANG=SCORE=");
 int e;
 

 if ((e=args.isError())>0)
    GError("%s\nInvalid argument: %s\n", usage, argv[e]);

 if (args.getOpt('h')!=NULL) GError("%s\n", usage);

 GStr infile;
 if (args.startNonOpt()) {
        infile=args.nextNonOpt();
        GMessage("Given file: %s\n",infile.chars());
        }

 asfasta=(args.getOpt('H')==NULL);
 by_length=(args.getOpt('L')!=NULL);
 do_share=(args.getOpt('n')!=NULL);
 //fullNP=(args.getOpt('F')!=NULL);
 if (args.getOpt('p')!=NULL) prefix=(char*)args.getOpt('p');
                        //else prefix=NULL;
 transitive=(args.getOpt('t')!=NULL);
 wSinglets=(args.getOpt('S')!=NULL);
 GStr s=args.getOpt('x');
 FILE* fxclude=NULL;
 GMALLOC(inbuf, inbuf_len); 
 if (!s.is_empty()) {
   if ((fxclude=fopen(s, "r"))==NULL)
      GError("Cannot open exclusion file '%s'!\n", s.chars());
   int c=readNames(fxclude, excludeList);
   fclose(fxclude);
   GMessage("Loaded %d sequences to be excluded.\n", c);
   flt_Exclude=(c>0);
   }   
 GStr clsfile=args.getOpt('c');
 FILE* fclusters=NULL;
 if (!clsfile.is_empty()) {
   if ((fclusters=fopen(clsfile, "r"))==NULL)
      GError("Cannot open input cluster file '%s'!\n", clsfile.chars());
   int c=readClusters(fclusters);
   if (c<1) GError("Error: no clusters can be loaded from file '%s'!\n", clsfile.chars());
   wclusters=true;
   GMessage("Loaded %d sequence within %d clusters\n", all.Count(), c);
   s=clsfile+".scls";
   if ((fscls=fopen(s, "w"))==NULL)
      GError("Cannot open seeded cluster file '%s' for writing!\n", s.chars());
   }
 GStr lytfile=args.getOpt('y');
 if (lytfile.is_empty() && wclusters) {
   lytfile=clsfile+".lyt";
   }
 FILE* flyt=NULL;
 if (!lytfile.is_empty()) {
   if ((flyt=fopen(lytfile, "w"))==NULL)
      GError("Cannot open layout file '%s' for writing!\n", lytfile.chars());
   }
   
 FILE* fgenes=NULL;
 s=args.getOpt('l');
 if (!s.is_empty()) {
   if ((fgenes=fopen(s, "r"))==NULL)
      GError("Cannot open full length sequences list file '%s'!\n", s.chars());
   int c=readNames(fgenes, fullGenes);   
   fclose(fgenes);
   if (c>0) {
     namesLoaded=true;
     GMessage("Loaded a list of %d full-length sequences\n", c);
     }
   }
 
 s=args.getOpt('w');
 if (!s.is_empty()) {
   if ((fgenes=fopen(s, "r"))==NULL)
      GError("Cannot open full length sequences list file '%s'!\n", s.chars());
   int c=readWeights(fgenes, fullGenes);
   fclose(fgenes);
   GMessage("Loaded a list of %d gene seqs and their weights\n", c);
   if (c>1) npWeights=true;
   }

 
 /*s=args.getOpt("SEQFLT");
 if (!s.is_empty()) {
   if (s=="ET") {
     flt_ET_only=true;
     s=" ETs/NPs only (discard ESTs)";
     }
   else if (s=="EST") {
     flt_EST_only=true;
     s=" EST only (discard ETs/NPs)";
     }
   else if (s=="EST2ET" || s=="ET2EST")  {
     flt_EST2ET=true;
     s=" links to ETs/NPs only (discard EST-EST links)";
     }
   else GError("%s\nIncorrect SEQFLT option ('%s' not recognized)", s.chars());
   GMessage("Using sequence filter: %s\n", s.chars());
   ETcheck=true;
   }*/
 s=args.getOpt("PID");
 if (!s.is_empty()) {
    minpid = s.asInt();
   }
 if (minpid>0) GMessage("PID=%d\n",  minpid);
 s=args.getOpt("SCOV");
 if (!s.is_empty()) {
    minscov = s.asInt();
   }
 if (minscov>0) GMessage("SCOV=%d\n",  minscov);
 s=args.getOpt("OVHANG");
 if (!s.is_empty()) {
    maxovhang = s.asInt();
    if (maxovhang<INT_MAX) GMessage("OVHANG=%d\n",  maxovhang);
   }
 s=args.getOpt("SCORE");
 if (!s.is_empty()) {
    minscore = s.asInt();
   }
 if (minscore>0) GMessage("SCORE=%d\n",  minscore);
 GStr outfile=args.getOpt('o');
 if (!outfile.is_empty()) {
     outf=fopen(outfile, "w");
     if (outf==NULL)
        GError("Cannot open file %s for writing!\n",outfile.chars());
     }
   else outf = (wclusters) ? NULL : stdout;

 s="Type of clustering: ";
 s+=(do_share)?" shared inclusions, ": "unique container, ";
 s+=(transitive)?" transitive containment. ":" direct inclusions only.";
 s+="\n";
 GMessage(s.chars());
 //==
 FILE* inf;
 if (!infile.is_empty()) {
    inf=fopen(infile, "r");
    if (inf==NULL)
       GError("Cannot open input file %s!\n",infile.chars());
    }
  else
   inf=stdin;
   
 GStr pfile=args.getOpt('f'); //write filtered pairs here
 FILE* fpairs=NULL;
 if (!pfile.is_empty()) {
   if (pfile == "-") 
      fpairs=stdout;
    else 
     if ((fpairs=fopen(pfile, "w"))==NULL)
      GError("Cannot write filtered hits file '%s'!", pfile.chars());
   }
 s=args.getOpt('d'); //write filtered pairs here
 if (!s.is_empty()) {
     if ((fTrees=fopen(s, "w"))==NULL)
      GError("Cannot write debug tree file '%s'!", s.chars());
   }
  
 //======== main program loop
 char* line;
 long fpos;
 while ((line=fgetline(inbuf, inbuf_len, inf, &fpos))!=NULL) {
        //GMessage("line is '#%s#'\n", line);        
        GStr linecpy(line);
        if (strlen(line)<4 || line[0]=='#') continue;
        char* tabpos=line;
        //find the 1st tab
        while (*tabpos != '\t' && *tabpos!='\0') tabpos++;
        if (*tabpos=='\0' || tabpos==line)
            GError(ERR_INVALID_PAIR, line);
        *tabpos='\0'; //so line would be the first node name
        if (flt_Exclude && excludeList.hasKey(line)) continue;
        tabpos++; //tabpos is now on the first char of the second field (q_len)
        int score, scov, ovhang, pid;
        double score2;
        //skip 3 other tabs delimited
        //read the query length:
        int qlen=getNextValue(tabpos, line);
        tabpos++;
        int q5=getNextValue(tabpos,line);
        tabpos++;
        int q3=getNextValue(tabpos,line);
        tabpos++;
        if (q5==q3) GError(ERR_INVALID_PAIR, line);
        bool minus=false;
        if (q5>q3) {
          minus=true;          
          Gswap(q5,q3);
          }

        //now we should be on the first char of the hitname field
        while (isspace(*tabpos)) tabpos++; //skip any spaces in this second node name field
        if (*tabpos=='\0') GError(ERR_INVALID_PAIR, line);
        //add a string termination after this field
        char* p=tabpos; while (!isspace(*p) && *p!='\0') p++;
        *p='\0';
        //now tabpos contains the exact second sequence string
        if (strcmp(line, tabpos)==0) {
          //GMessage("Warning: self pairing found for node %s\n",line);
          continue;
          }
        if (flt_Exclude && excludeList.hasKey(tabpos)) continue;
        //if (!seq_filter(line, tabpos)) continue;
        p++; //move on the first char of the hitlen 
        int hitlen=getNextValue(p,line);
        p++;
        int h5=getNextValue(p,line);
        p++;
        int h3=getNextValue(p,line);
        p++;
        if (h5==h3) GError(ERR_INVALID_PAIR, line);
        if (h5>h3) {
          minus = !minus;
          Gswap(h5,h3);
          }
        char* qgaps=NULL;
        char* hgaps=NULL;  
        pid=getNextValue(p,line); p++;
        score=getNextValue(p,line);p++;
        //parse the 2nd score, the orientation and the gapping info if any!
        score2=getNextNumber(p,line);p++;
        if (flyt!=NULL && !transitive) {        
          char sense=*p;
          if (sense!='+' && sense!='-')
              GError("Error parsing hit orientation at (p='%s'): %s\n",sense, p, linecpy.chars());
          p++;
          if (*p=='\t') {//gapinfo present
            p++;qgaps=p;
            while (*p!='\t' && *p!='\n' && *p!=0) p++;
            if (*p=='\t') {
                if (qgaps==p) qgaps=NULL;
                         else *p=0;
                p++;
                }
            hgaps=p;
            while (*p!='\t' && *p!='\n' && *p!=0) p++;
            if (hgaps==p) hgaps=NULL;
                     else *p=0;
            }
          }//gapinfo can be stored
        
        //compute coverages:
        if (hitlen>qlen) { //query is shorter
          scov = iround(((double)(q3-q5)*100)/qlen);
          //lcov = (int) rint(((double)(h3-h5)*100)/hitlen);
          ovhang=GMAX(q5-1, qlen-q3);
          }
        else {//hit sequence is shorter
          scov = iround(((double)(h3-h5)*100)/hitlen);
          ovhang=GMAX(h5-1, hitlen-h3);
          }
        //GMessage("scov=%d, pid=%d, score=%d, ovhank=%d\n", scov, pid, score, ovhang);
        CNode* n1=NULL;
        CNode* n2=NULL;
        if (wclusters){
          if ((n1=all.Find(line))==NULL || (n2=all.Find(tabpos))==NULL) continue; //ignore pair!
          if (n1->cluster!=n2->cluster) continue; //ignore pair, sequences not in the same cluster!
          }
        if (scov>=minscov && pid>=minpid && score>=minscore && ovhang <= maxovhang) {
           //GMessage("%s(%d) %d-%d  | %s(%d) %d-%d, pid=%d%%, score=%d, scov=%d%%, ovhang=%d\n",
           //          line,qlen,q5,q3,tabpos,hitlen, h5,h3,pid, score, scov, ovhang);
           if (fpairs!=NULL) fprintf(fpairs, "%s\n", linecpy.chars());
           COverlap ovl(line, tabpos, qlen, hitlen, q5, h5, 
              q3-q5+1, h3-h5+1, score, pid, minus, qgaps, hgaps);
              
           addPair(ovl);
           }
         /*else { //yet, if any of these are full genes, keep them here:
           if (!all.hasKey(line) && (wSinglets || fullGenes.hasKey(line)))
              addNode(line, qlen);
           if (!all.hasKey(tabpos) && (wSinglets || fullGenes.hasKey(tabpos)))
              addNode(tabpos, hitlen);
           }*/
    } //while lines are coming
    
  if (inf!=stdin) fclose(inf);
  if (fpairs!=NULL && fpairs!=stdout) fclose(fpairs);
  //==========
  //   
 GMessage("Total sequences analyzed: %d\n", seqs.Count());
 GMessage("Sorting clusters by size..\n");
 //sort by cluster size (largest first)
 //if no transitivity is assumed, parent will be set to NULL
 //for any grand-child 
 for (int i=0; i<seqs.Count(); i++) {
    //if (seqs[i]->links.size<=1) break; //the rest are already included in other clusters 
    CNode* seq=seqs[i];
    if (!transitive && seq->parent!=NULL  && seq->parent->parent!=NULL)
        seq->parent=NULL; //set free if it's a second generation
    if (seq->parent!=NULL)
      continue; //only process sequences with no parent at this stage
    GCluster* C = new GCluster(true,false,false);
    collectNodes(C, seqs[i]);
    seqs[i]->numseqs=C->Count();
    delete C;
    }
 seqs.setSorted(&compareSize);
 if (wSinglets) 
   GMessage("Writing clusters (and singletons)..\n");
  else  
   GMessage("Writing clusters..\n");
  //defline: >NRCL# <numseqs> <refseq> <refseqlen> <flags>
  //<flags> can be: c (possibly chimeric parent), F (contains full NPs) | G (contains non-full NPs)
  int sngnum=0;
  int clnum=0;
  for (int i=0; i<seqs.Count(); i++) {
    CNode* seq=seqs[i];
    //GMessage("Collecting all children for '%s'\n", seq->id);
    if (seq->parent!=NULL) //sequence was in a cluster, skip it
                   continue; 
    if (seq->numseqs==1 && !wSinglets)
        break; //don't write "singletons"
    GCluster* C=new GCluster(true,false,false); //sorted by pointer
    //descend the full children tree, if containment is transitive
    collectNodes(C, seq, true);
    //----
    //GStr flags;
    if (C->Count()>=1) {
      GStr refSeq(seq->id);
      int refLen=seq->len;
      C->setSorted(&compareLen);//sort by sequence length (largest first)
      //decide what is the refseq to use
      /*bool isEST=(seq->seqtype==0);
      //scan all sequences for NPs/ETs
      bool hasFull=false;
      bool hasET=false;
      //scan to see if there are any NPs/ETs here
      for (int k=0; k<C->Count(); k++) {
        CNode* n=C->Get(k);
        if (n->seqtype>0) { //NP or ET!
          hasET=true;
          }
        if (!hasFull && fullGenes.hasKey(n->id)) hasFull=true;
        if (hasFull && hasET) break; //no need to check further
        }
      //set the flags accordingly
      if (C->Count()>1 && isEST && refLen - C->Get(1)->len > 50 &&
              (C->Count()>=20 || hasET))
                 flags="C"; //EST much longer than the rest of the herd
      if (hasFull) flags+="F";
              else if (hasET) flags+="G";*/
      if (asfasta) {
         if (seq->numseqs>1) {
             clnum++;
             if (outf!=NULL) fprintf(outf,">NRCL%d %d %s %d\n", clnum,
                                seq->numseqs, refSeq.chars(), refLen);
             }
           else {
             sngnum++;
             if (outf!=NULL) fprintf(outf,">SNGT%d %d %s %d\n", sngnum,
                           seq->numseqs, refSeq.chars(), refLen);
             }
         }
         /*if (outf!=NULL)  {
            if (flags.length()>0) fprintf(outf, " %s\n", flags.chars());
                             else fprintf(outf, "\n");
             }*/
      if (flyt!=NULL && C->Count()>1)
        fprintf(flyt, ">%s %d 1 %d\n", seq->id, C->Count(), seq->len);
      
      for (int j=0; j<C->Count();j++) {
          CNode* c=C->Get(j);
          if (wclusters && seq!=c)
               seq->cluster->Remove(c);
           if (outf!=NULL) {    
             if (j==0) 
               fprintf(outf,"%s",c->id);
             else 
               fprintf(outf," %s",c->id);
             }
          if (flyt!=NULL && C->Count()>1) {
            fprintf(flyt, "%s %c %d %d %d %d",
               c->id, c->lytminus?'-':'+',
               c->len, c->lytpos+1, c->pclipL, c->pclipR);
            if (c->gaps!=NULL || c->pgaps!=NULL) {
              fprintf(flyt," R:");
              if (c->gaps!=NULL) fprintf(flyt,"%s",c->gaps);
              if (c->pgaps!=NULL) fprintf(flyt,"/%s",c->pgaps);
              }
            fprintf(flyt,"\n"); 
            }
          } //for each component sequence
      if (outf!=NULL) fprintf(outf,"\n");
      }
    if (wclusters) {
             //C->setSorted(&comparePtr);
             //C->RemovePtr(seq);
             seq->contained=C;
             //seq->contained->setSorted(&compareSize);
             }
        else delete C;
    } //for each sequence
 if (wSinglets)
   GMessage("%d clusters written (and %d singletons).\n", clnum, sngnum);
  else  
   GMessage("%d clusters written.\n", clnum);
  //---
  if (wclusters) { //write the clusters into scls file
    for (int i=0;i<tclusters.Count();i++) {
      GCluster* c=tclusters[i];
      c->setSorted(&compareContained);
      fprintf(fscls,">SCL%d\t%d\t%d\n", i+1, c->total, c->Count());
      for (int j=0; j<c->Count(); j++) {
         CNode* n=c->Get(j);
         if (n->contained!=NULL) {
           //fprintf(fscls, "%s :%d:", n->id, n->contained->Count());
           fprintf(fscls, "%s", n->id);
           for (int k=0;k<n->contained->Count();k++) {
              CNode* cn=n->contained->Get(k);
              if (cn!=n)
                fprintf(fscls, " %s", cn->id);
              }
           fprintf(fscls, "\n");
           }
          else
           fprintf(fscls, "%s\n", n->id);
         }
      }
   fflush(fscls);
   fclose(fscls);
   }
  if (outf!=NULL) fflush(outf);
  if (flyt!=NULL) {
   fflush(flyt);
   fclose(flyt);
   }
  seqs.Clear();
  tclusters.Clear();
  all.Clear(); //this frees the nodes
  //GMessage("the nrcls list was cleared!\n");
  if (outf!=stdout && outf!=NULL) fclose(outf);
  GFREE(inbuf);
  //GMessage("*** all done ***\n");
  //getc(stdin);
}

// ====================================== CNODE methods

CNode::CNode(char* name, int slen, int ordn) {
    kids=new GCluster(true,false,false);
    id=Gstrdup(name);
    //seqtype=(isET(id) ? 1 : 0);; //1=ET 0=EST
    pweight=0;
    isFull=startsWith(name, prefix);
    if (namesLoaded && fullGenes.hasKey(name)) isFull=true;
    int *pw=NULL;
    if (npWeights && (pw=fullGenes.Find(name))!=NULL) {
      isFull=true;
      pweight=*pw;
      }
    numseqs=1; //self, by default being a singleton
    ordnum = (ordn>=0)? ordn : base_ordnum++;
    len=slen;
    parent=NULL;
    cluster=NULL;
    gaps=NULL;
    pgaps=NULL;
    poffs=0;
    lytminus=0;
    lytpos=0;
    pclipL=0;pclipR=0;
    pscore=0;
    contained=NULL;
 }
/*-- node creation when loading from cluster file --*/

CNode::CNode(char* name, GCluster*& cl, int ordn) {
    kids=new GCluster(true, false, false);
    id=Gstrdup(name);
    pweight=0;
    //seqtype=(isET(id) ? 1 : 0);; //1=ET 0=EST
    isFull=startsWith(name, prefix);
    if (namesLoaded && fullGenes.hasKey(name)) isFull=true;
    int *pw=NULL;
    if (npWeights && (pw=fullGenes.Find(name))!=NULL) {
      isFull=true;
      pweight=*pw;
      }
    numseqs=1; //self, by default being a singleton
    len=0; /* not initialized yet! */
    parent=NULL;
    gaps=NULL;
    pgaps=NULL;
    poffs=0;
    lytminus=0;
    lytpos=0;
    pclipL=0;pclipR=0;
    ordnum=ordn;
    pscore=0;
    if (cl==NULL) {
       cl=new GCluster(true,false,false);
       tclusters.Add(cl);
       cluster=cl;
       }
    cluster=cl;
    cluster->addNode(this);
    contained=NULL;
    }
     
   //directly set the parent as specified, with no other checking

void CNode::setParent(CNode* np, int offs, int clpL, int clpR,
                    int score, bool minus, char* g, char* pg) {
   if (parent!=NULL) {
         if (!do_share) //the previous parent should ditch this child
            parent->kids->Remove(this);
         /* if transitivity is not
         if (parent->parent!=NULL && !transitivity)
            return; //will not set the parent */
         } //switch parent
   parent=np;
   poffs=offs;
   pclipL=clpL;
   pclipR=clpR;
   pscore=score;
   replaceStr(gaps, g);
   replaceStr(pgaps, pg);
   pminus=minus?1:0;
   parent->kids->Add(this);
   /* when no transitivity is allowed,
      any immediate kids should be let go */
   /*
   if (!transitive)
     for (int i=0;i<kids.Count();i++)
       kids[i]->parent=NULL;
   */
   }

//check if a sequence np is to be assigned as "parent" to this
//returns true if parenting works according to the sequence AND/OR link weights
//or false if the parenting is reversed!
bool CNode::assignParent(CNode* np, int score, int offs, int clpL, int clpR,
                          bool minus, char* g, char* pg) {
   //returns true if the parent was set to the new node np
    if (np->parent==this || parent==np) return false; //avoid circular or double relations
    //GMessage(" assign '%s' as parent for '%s'\n", np->id, id);
    
    //do not allow an ET to be parented by an EST - only the other way around!
    //if (np->seqtype < seqtype) return false;
    if (parent==NULL || do_share) {
      setParent(np, offs, clpL, clpR, score, minus, g, pg);
      return true;
      }
    //switch to a better scoring parent link
    bool switchParent=false;
    //int *pw=NULL, *npw=NULL;
    
    //bool isFullparent=startsWith(parent->id, "np|");
    //bool isFullnp=startsWith(np->id, "np|");
    
    if (parent->pweight>0 || np->pweight>0) {
     switchParent= (np->pweight>parent->pweight);
     }
    else if (parent->isFull ^ np->isFull) {
      switchParent= (np->isFull && !parent->isFull);
      }
    else if (by_length && np->len!=parent->len) { // the length of the parent take precendence!
       switchParent = (np->len>parent->len);
       }
     else { // the higher similarity
      if (score == pscore) {
       switchParent = np->kids->Count()==parent->kids->Count() ?
           (np->len > parent->len) : (np->kids->Count()>parent->kids->Count());
       }
      else switchParent = score>pscore;
      }
       
    if (switchParent) { //better parent found:
      setParent(np, offs, clpL, clpR, score, minus, g, pg);
      return true;
      }
    return false;   //no reason to change the parent!
    }
    
CNode::~CNode() {
    if (gaps!=NULL) GFREE(gaps);
    if (pgaps!=NULL) GFREE(pgaps);
    if (contained!=NULL)
       delete contained;
    delete kids;   
    GFREE(id);
    }



//
int getNextValue(char*& p, char* line) {
 //expects p to be on the first char of a numeric field
 char* start, *end;
 while (isspace(*p) && *p!='\0') p++; //skip any spaces;
 if (*p=='\0') GError(ERR_INVALID_PAIR, line);
 start=p;
 while (*p!='\t' && *p!='\0') p++;
 //now p should be on the next tab or eos;
 end=p-1;
 while (isspace(*end)) end--;
 //char c=*(end+1);
 *(end+1)='\0';
 double d=atof(start);
 return (int)d;
 //backup
 }

double getNextNumber(char*& p, char* line) {
 //expects p to be on the first char of a numeric field
 while (isspace(*p) && *p!='\0') p++; //skip any spaces;
 if (*p=='\0') GError(ERR_INVALID_PAIR, line);
 char* endptr;
 double d=strtod(p, &endptr);
 if (endptr==p) GError(ERR_INVALID_PAIR, line);
 while (isspace(*endptr)) endptr--;
 //char c=*(end+1);
 p=endptr+1;
 *p='\0'; //destructive ending
 return d;
 //backup
 }


void showCl(GCluster* C) {
 GMessage("{%s",C->Get(0)->id);
  for (int i=1;i<C->Count();i++) 
     GMessage(" %s",C->Get(i)->id);
 GMessage("}\n");
}


void collectNodes(GCluster* C, CNode* n, bool dbg) {
 if (transitive) {
   // recursively descend to all links and add all the nodes to C
   if (C->IndexOf(n)<0) C->Add(n); 
                   else return; 
   if (dbg && fTrees!=NULL)
       if (n->kids->Count()>0) 
               fprintf(fTrees, "%s:%d{", n->id, n->kids->Count());
          else fprintf(fTrees, "%s ", n->id);
     //do not traversed a branch already traversed
   for (int i=0;i<n->kids->Count(); i++) { 
       CNode* k=n->kids->Get(i);
       if (k->parent!=n)
         GError("Error: kid '%s' should have '%s' as parent!\n",
             k->id, n->id);
       k->lytminus=(n->lytminus ^ k->pminus);
       k->lytpos=n->lytpos + 
             ((n->lytminus) ? n->len-k->len-k->poffs : k->poffs );
       if (n!=k) collectNodes(C, k,dbg);
       }
   if (dbg && fTrees!=NULL) {
     if (n->kids->Count()>0) fprintf(fTrees, "}");
     if (n->parent==NULL) fprintf(fTrees, "\n");
     }
   } //transitive case
  else {//not transitive, collect only n and its immediate children
   C->Add(n);
   for (int i=0;i<n->kids->Count(); i++) {
       CNode* k=n->kids->Get(i);
       if (k->parent!=n && k->parent!=NULL) {
           if (!do_share) GMessage(" Error: kid '%s' should have '%s' as parent!\n"
                "(but the parent is set to '%s')\n", k->id, n->id, k->parent->id);
         }
       if (k->parent==NULL) continue;
       k->lytpos=k->poffs;
       k->lytminus=k->pminus;
       if (n!=k) C->Add(k);
       }
   
   }
}


//add a new node with name id1, with length len1 
CNode* addNode(char* id1, int len1, int ordn) {
 //use the same cluster
 CNode* n1=new CNode(id1,len1, ordn); //this automatically adds n1 to the shared t-cluster 
 all.shkAdd(n1->id,n1);
 seqs.Add(n1);
 return n1;
}
/*
//merge the links for two existing sequences, a longer one against a shorter one
CNode* mergeLinks(CNode* nshort, CNode* nlong) {
//nlong has always more links than the nshort
GCluster* C=NULL;
//if ((C=nshort->links.isIncluded(&nlong->links))!=NULL)
  //C->size=0; //meaning: it was included in a larger cluster, 
                   //so it's not the non-redundification result
return ((C==&nshort->links) ? nshort : nlong);
}
*/

/*bool seq_filter(char* s1, char* s2) {
 bool isET1, isET2;
 if (ETcheck) {
   isET1=isET(s1);
   isET2=isET(s2);
   if (flt_ET_only) 
     return (isET1 && isET2); //both are et
   else if (flt_EST_only)
     return (!isET1 && !isET2);
   else if (flt_EST2ET)
     return (isET1 || isET2);
   }
 return true; //passed
 }*/
 
void addPair(COverlap& o) {
 o.n1=all.Find(o.id1);
 o.n2=all.Find(o.id2);
 if (wclusters) {
   if (o.n1->len==0 && o.n1!=NULL) o.n1->len=o.len1;
   if (o.n2->len==0 && o.n2!=NULL) o.n2->len=o.len2;
   }
 bool new_n1=false;  
 if (o.n1==NULL) {
        o.n1 = addNode(o.id1, o.len1);
        new_n1=true;
        }
    else 
        if (o.n1->len!=o.len1) GError("Length mismatch for '%s'\n", o.id1);

 if (o.n2==NULL) {
       o.n2 = new_n1 ? addNode(o.id2, o.len2, o.n1->ordnum): addNode(o.id2, o.len2);
       }
    else 
       if (o.n2->len!=o.len2) GError("Length mismatch for '%s'\n", o.id2);

 if (o.n1==o.n2) return;
 // Now the big parenting decision:
 //-- will arrange n1 and n2 as if n1 is the parent of n2:
 int offs;
 int clpL;
 int clpR;
 /*int o1start=o.o1start, o2start=o.o2start, 
     ovl1=o.ovl1, ovl2=o.ovl2;*/
 //decide who's the boss (parent), and clip as needed
 //int *w1=NULL, *w2=NULL; //weights
 bool nswap=false; //swap such that n1 is the "parent"
 //bool isFull1=startsWith(o.n1->id, "np|");
 //bool isFull2=startsWith(o.n2->id, "np|");
 int n1kids = o.n1->kids->Count();
 int n2kids = o.n2->kids->Count();
 if (o.n1->pweight>0 || o.n2->pweight>0) {
    nswap=(o.n1->pweight<o.n2->pweight);
    }
   else if (o.n2->isFull ^ o.n1->isFull) {
    nswap=(o.n2->isFull && !o.n1->isFull);
    }
   else if (by_length) { //pick by length
    nswap=(o.len1<o.len2);
    }
   else if (n1kids!=n2kids) {//by the number of children they already have
    nswap=(n2kids>n1kids);
    }
   else {
    nswap=(o.n1->ordnum==o.n2->ordnum) ? (o.len1>o.len2) :
                               (o.n1->ordnum>o.n2->ordnum);
     //the earliest encountered node should parent
    }
 if (nswap) o.Swap();
   /*{ n1=o.n2; n2=o.n1;
   o1start=o.o2start; o2start=o.o1start;
   ovl1=o.ovl2;ovl2=o.ovl1; }*/
   
  //clipping -- always clip the child, not the parent 
 //GMessage("Add %s[%d] to %s[%d] (gaps: %s/%s)\n", o.n2->id, o.n2->len, o.n1->id, o.n1->len, o.gaps2,o.gaps1);
 if (o.minus) {
   clpL=o.n2->len-(o.o2start-1+o.ovl2);
   offs=o.o1start-1-clpL;
   clpR=o.o2start-1;
   }
  else {
   clpL=o.o2start-1;
   clpR=o.n2->len-(o.o2start-1+o.ovl2);
   offs=o.o1start-o.o2start;
   }
 if (((o.len2-clpL-clpR)*100)/o.len2 >= minscov)
   o.n2->assignParent(o.n1, o.score, offs, clpL, clpR, o.minus, o.gaps2, o.gaps1);
 //TRY to set n1 as the parent of n2, but with the provision for the EXISTING parent of n2!
 //(if decided, the previous parent of n2, if any, should abandon n2
 }

#ifdef DBGTRACE
void memlog(int level, const char* msg1,
         GCluster* C=NULL,
         const char* msg2=NULL) {
 GStr s;
 if (level>0) {
   s.format("%d>", level);
   s.padR(level+2);
   }
 if (msg1!=NULL) s+=msg1;
 if (C!=NULL) {
     GStr s2;
     s2.format("(%d){%s",C->Count(), C->Get(0)->id);
     s+=s2;
     for (int i=1;i<C->Count();i++)  {
        s2.format(" %s",C->Get(i)->id);
        s+=s2;
        }
     s+="}";
     }
 if (msg2!=NULL) s+=msg2;
 mem.log(s.chars()); 
 }
#endif  

int readNames(FILE* f, GHash<int>& xhash) {
  int c;
  int count=0;
  char name[256]; 
  int len=0;
  while ((c=getc(f))!=EOF) {
    if (isspace(c)) {
      if (len>0) {
        name[len]='\0';      
        xhash.Add(name, new int(1));
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

//-- reads seq_name followed by an integer, per each line
int readWeights(FILE* f, GHash<int>& xhash) {
  int count=0;
  char* line;
  long fpos;
  while ((line=fgetline(inbuf, inbuf_len,f, &fpos))!=NULL) {
   if (strlen(line)<=1 || line[0]=='#') continue;
   char* spcpos=line;
   //find the 1st space to split there
   while (*spcpos!='\t' && *spcpos!=' ' && *spcpos!='\0') spcpos++;
    if (*spcpos=='\0' || spcpos==line) {
             GMessage(ERR_INVALID_PAIR, line); continue; 
             }
   *spcpos='\0';
   spcpos++;
   int weight=getNextValue(spcpos,line);
   xhash.Add(line, new int(weight));
   count++;
  }
 return count;   
}


//============ comparison functions:
int comparePtr(void* p1, void* p2) {
  return (p1>p2) ? -1 : ((p1<p2)?1:0);
  }
//compare sequences by length

int compareSize(void* p1, void* p2) { //by valid cluster size (all nodes)
 int c1=((CNode*)p1)->numseqs;
 int c2=((CNode*)p2)->numseqs;
 return (c1>c2)?-1:((c1<c2)?1:0);
}

int compareLen(void* p1, void* p2) { //by sequence length (links)
 int c1=((CNode*)p1)->len;
 int c2=((CNode*)p2)->len;
 return (c1 > c2) ? -1 : ((c1<c2)? 1:0);
}

//=====
int compareID(void* p1, void* p2) {
 return strcmp(((CNode*)p1)->id,((CNode*)p2)->id);
}

int compareContained(void* p1, void* p2) { //by number of contained seqs
 GCluster* c1=((CNode*)p1)->contained;
 GCluster* c2=((CNode*)p2)->contained;
 int v1=(c1==NULL) ? 0: c1->Count();
 int v2=(c2==NULL) ? 0: c2->Count();
 return (v1 > v2) ? -1 : ((v1<v2)? 1:0);
}

int readClusters(FILE* f) {
  char *p;
  char buf[256];
  int c;
  int len, count=0, scount=0;
  int numseqs;
  //assumes: first line is always a defline 
  //         deflines are always shorter than 255
  //         sequence names are always shorter than 255 chars
  CNode* seq=NULL;
  while (!feof(f)) {
   if (fgets(buf,255,f)==NULL) break;
   len=strlen(buf);
   if (buf[len-1]!='\n') GError(ERR_CL_PARSE, buf); //very long defline ?!?
   buf[len-1]='\0';
   if (buf[0]!='>') GError(ERR_CL_PARSE,buf);
   p=buf+1;
   p=strpbrk(p,"\t "); //skip NRCL#
   if (p==NULL) GError(ERR_CL_PARSE,buf);
   while (isspace(*p) && *p!='\0') p++;
   //s to the first nonspace
   if (*p=='\0') GError(ERR_CL_PARSE,buf);
   numseqs=0;
   numseqs=getNextValue(p, buf); //should return number of seqs within this cluster
   if (numseqs<=0) GError(ERR_CL_PARSE,buf);
   count++;
   //--now read all component sequences for this cluster
   len=0;
   scount=0;
   GCluster* curcl=NULL;
   while ((c=getc(f))!=EOF) {
    if (isspace(c)) { //word separator found
      if (len>0) { //token found: it must be a sequence name
        buf[len]='\0';
        scount++;
        seq=new CNode(buf, curcl, base_ordnum);
        all.shkAdd(seq->id, seq); //add it to the global hash too
        seqs.Add(seq);

        len=0;
        }
      if (c=='\n') {
        c=getc(f);
        if (c==EOF) break;
        if (c=='>') {
          ungetc(c, f);
          break;
          }
        }
      continue;
      } //isspace
    //a non-space: append to the current name
    buf[len]=(char) c;
    len++;
    if (len>255) {
      buf[len-1]='\0';
      GError("Error reading seqname: too long ('%s') !\n",buf);
      }
    } //while getc
   base_ordnum++; 
   //GMessage("numseqs=%d, scount=%d\n", numseqs, scount);
   if (scount!=numseqs) GError(ERR_CL_PARSE, seq->id);
   } //while
 fclose(f);
 return count;
}
