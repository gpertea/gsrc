#include <ctype.h>
#include "gcl/GBase.h"
#include "gcl/GArgs.h"
#include "gcl/GString.h"
#include "gcl/GHash.hh"
#include "gcl/GList.hh"

#define usage "Perform transitive-closure clustering with layout consistency \n\
 by filtering tabulated hits. Usage:\n\
 lclust [<hits_file>] [-H] [-f <flthits_file>] [-o <out_file>] [-y <lyt_file>]\n\
 [-x <xcludelist>] [-s <seqlist>] [-r <restrictlst>] [-c <clone_lines>]\n\
 [SEQFLT={ET|EST|ET2EST}] [-x <excludefile>] [SCORE=xx] [PID=xx] [OVL=xx]\n\
 [COV=xx] [OVHANG=xx] \n\
 Options:\n\
 -H  : disable the fasta-style header (with node count info) for output clusters\n\
 -o  : write output in <out_file> instead of stdout\n\
 -f  : write to <flthits_file> all the lines that passed the filters\n\
 -y  : write cluster layouts to <lyt_file>\n\
 \n\
 Filters (optional, can be combined):\n\
 -s use only pairs involving at least one sequence from <seqlist> \n\
 -r use only hits between sequences from list <restrictlst> \n\
 -c use clone sets from file <clone_lines>\n\
 -x load sequence names from <excludefile> and discard any pairs related to them \n\
 SEQFLT=   filter sequence names for specific sequence type clustering: \n\
           EST  - EST-only clustering (use only ESTvsEST hits) \n\
           ET   - ETs/NPs only clustering\n\
           ET2EST - ignore EST-EST pairs, uses only EST-ET or ET-ET pairs\n\
 OVL=xx    minimum overlap (bp) (default=20)\n\
 COV=xx    is the minimum coverage of a sequence  \n\
           (minimum overlap length as percentage of the shorter sequence)\n\
 OVHANG=xx is the maximum overhang (in nucleotides) to be allowed\n\
           for any overlap(default 1000)\n\
 PID=xx   is the minimum percent of identity\n\
 SCORE=xx is the minimum score\n\
 \n\
 If <hits_file> is not given, hits are expected at stdin\n\
 The program expects each input line to have these tabulated fields:\n\n\
 q_name q_len q_n5 q_n3 hit_name hit_len hit_n5 hit_n3 pid score1 score2 strand\n\
"
#define BUF_LEN 8192 //maximum input line length 
#define ERR_INVALID_PAIR "Invalid input line encountered:\n%s\n"

/*---
 OVHCOV=xx is the maximum percent length of an overhang \n\
           (maximum overhang length as percentage of participant sequence length)\n\

*/


//============ structures:
//define a cluster object as a list of node names (GStrings)

// dumb node sorting function -- by pointer value (!)
int comparePtr(void* p1, void* p2) {
  return (p1>p2) ? -1 : ((p1<p2)?1:0);
  }

class CNodeBase {
 public:
   char* id; //actual node ID
   int ovx5; // left overhang extension
   int ovx3; // right overhang extension
   int trim5; // left overhang extension
   int trim3; // right overhang extension
   // positioning:
   //CNode* nto; // node which the offset nofs is relative to
   bool nrev; // reverse complement positioning to nto?
   int seqlen; // read length
   //int nofs; // offset in the current cluster
   //
public:
    CNodeBase() {
      id=NULL;
      ovx5=0;
      ovx3=0;
      trim5=0;
      trim3=0;
      seqlen=0;
      //nofs=0;
      //nto=NULL;
      nrev=false;
      }
  //comparison: just pointer based:
    bool operator==(CNodeBase& d){
     return (this==&d);
     }
    bool operator>(CNodeBase& d){
     return (this>&d);
     }
    bool operator<(CNodeBase& d){
     return (this<&d);
     }

};

class GCluster :public GList<CNodeBase> {
 public:
  bool operator==(GCluster& d){
     return (this==&d);
     }
  bool operator>(GCluster& d){
     return (this>&d);
     }
  bool operator<(GCluster& d){
     return (this<&d);
     }
  GCluster():GList<CNodeBase>(true,false,true) {
      }
  GCluster(bool sorted, bool free_elements=false, bool beUnique=true):
     GList<CNodeBase>(sorted,free_elements,beUnique) {
     
      }
  };

class CNode:public CNodeBase {
 public:
   GCluster* cluster; /* pointer to a list of nodes (the t-cluster
                         containing this sequence) */
   //GCluster links; //list of all nodes linked to this node
   CNode(char* name) {
    //links are sorted, not free, unique
    id=Gstrdup(name);
    ovx5=0;
    ovx3=0;
    trim5=0;
    trim3=0;
    seqlen=0;
    //nofs=0;
    //nto=NULL;
    nrev=false;
    cluster=new GCluster;
    cluster->Add(this);
    }
   ~CNode() { GFREE(id); }
   // create & add to an existing cluster
   CNode(char* name, CNode* pairnode) {
    id=Gstrdup(name);
    ovx5=0;
    ovx3=0;
    trim5=0;
    trim3=0;
    seqlen=0;
    //nofs=0;
    //nto=NULL;
    nrev=false;
    cluster=pairnode->cluster; //same cluster
    cluster->Add(this); //add self to the cluster
    }
 };

//-- global hash containing all the nodes
GHash<CNode> all;

bool asfasta;
bool flt_ET_only=false;
bool flt_EST_only=false;
bool flt_EST2ET=false;
bool flt_Exclude=false; //if an exclusion list was given
bool flt_seqOnly=false;
bool flt_seqRestrict=false;
bool ETcheck=false;
bool do_debug=false;

int minpid=0, minscov=0, minovl=20, maxovhang=100000, minscore=0;

static char inbuf[BUF_LEN]; // incoming buffer for sequence lines.
static int lineno;
FILE* outf;

//all clusters are stored in a unique sorted list
GHash<int> excludeList;
GHash<int> seqonlyList;
GList<GCluster> tclusters(true,true,true);

//=========== functions used:
bool addPair(char* s1, char* s2, bool minus, int qlen, int q5, int q3, int hlen, int h5, int h3);
int compareCounts(void* p1, void* p2);
int compareID(void* p1, void* p2); //compare node IDs
int getNextValue(char*& p, char* line);
bool seq_filter(char* s1, char* s2); //returns true if the pair can pass through the filter
int readNames(FILE* f, GHash<int>& xhash);
int readClones(FILE* f);

void printClusters() {
 for (int i=0; i<tclusters.Count(); i++) {
    GCluster* c=tclusters[i];
    c->setSorted(compareID);
    fprintf(stderr, ">Cluster[%d]%d\t%d\n", lineno, i, c->Count());
    for (int j=0;j<c->Count();j++) {
      CNode* n = (CNode*)c->Get(j);
      char sc=n->nrev?'-':'+';
      char ch = j%4?' ':'\n';
      fprintf(stderr, "(%c%s %d,%d {%d,%d})%c",
         sc,n->id,n->ovx5,n->ovx3,n->trim5, n->trim3, ch);
      }
    c->setSorted(true);  
    fprintf(stderr,"\n");
    }
}

//#define _DEBUG
//========================================================
//====================     main      =====================
//========================================================
int main(int argc, char * const argv[]) {
 GArgs args(argc, argv, "hHDf:c:o:y:s:r:x:SEQFLT=PID=COV=OVHANG=OVL=SCORE=DEBUG=");
 int e;
 if ((e=args.isError())>0)
    GError("%s\nInvalid argument: %s\n", usage, argv[e]);
 if (args.getOpt('h')!=NULL) GError("%s\n", usage);

 GString infile;
 if (args.startNonOpt()) {
        infile=args.nextNonOpt();
        GMessage("Given file: %s\n",infile.chars());
        }
 asfasta=(args.getOpt('H')==NULL);
 do_debug=(args.getOpt("DEBUG")!=NULL) || (args.getOpt('D')!=NULL);
 GString s=args.getOpt('x');
 FILE* fxclude=NULL;

 if (!s.is_empty()) {
   if ((fxclude=fopen(s, "r"))==NULL)
      GError("Cannot open exclusion file '%s'!\n", s.chars());
   int c=readNames(fxclude, excludeList);
   GMessage("Loaded %d sequences to be excluded.\n", c);
   fclose(fxclude);
   flt_Exclude=(c>0);
   }
 FILE* flyt=NULL;
 s= args.getOpt('y');
 if (!s.is_empty()) {
   if (s == "-")
      flyt=stdout;
    else
     if ((flyt=fopen(s, "w"))==NULL)
      GError("Cannot created layout file '%s'!", s.chars());
   }


 FILE* fseqonly=NULL;
 s=args.getOpt('s');
 if (!s.is_empty()) {
   if ((fseqonly=fopen(s, "r"))==NULL)
      GError("Cannot open exclusion file '%s'!\n", s.chars());
   int c=readNames(fseqonly, seqonlyList);
   GMessage("Loaded %d sequence names - only consider pairs including them.\n", c);
   flt_seqOnly=(c>0);
   fclose(fseqonly);
   }
 s=args.getOpt('r');
 if (!s.is_empty()) {
   if ((fseqonly=fopen(s, "r"))==NULL)
      GError("Cannot open exclusion file '%s'!\n", s.chars());
   seqonlyList.Clear();
   int c=readNames(fseqonly, seqonlyList);
   GMessage("Loaded %d sequence names - only consider hits between them.\n", c);
   flt_seqRestrict=(c>0);
   if (flt_seqRestrict) flt_seqOnly=false;
   fclose(fseqonly);
   }

 s=args.getOpt('c');
 if (!s.is_empty()) {
   if ((fxclude=fopen(s, "r"))==NULL)
      GError("Cannot open clone list file '%s'!\n", s.chars());
   int c=readClones(fxclude);
   GMessage("Loaded %d clones.\n", c);
   fclose(fxclude);
   }
 
 s=args.getOpt("SEQFLT");
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
   }
 s=args.getOpt("PID");
 if (!s.is_empty()) {
    minpid = s.asInt();
    if (minpid>0) GMessage("PID=%d\n",  minpid);
   }
 s=args.getOpt("COV");
 if (!s.is_empty()) {
    minscov = s.asInt();
    if (minscov>0) GMessage("SCOV=%d\n",  minscov);
   }

 s=args.getOpt("OVHANG");
 if (!s.is_empty()) {
    maxovhang = s.asInt();
   }
 if (maxovhang<1000) GMessage("OVHANG=%d\n",  maxovhang);   

 s=args.getOpt("OVL");
 if (!s.is_empty()) {
    minovl = s.asInt();
   }
 if (minovl>0) GMessage("OVL=%d\n",  minovl);   

 /*s=args.getOpt("OVHCOV");
 if (!s.is_empty()) {
    maxovhcov = s.asInt();
    if (maxovhcov<100) GMessage("OVHCOV=%d\n",  maxovhcov);
   }*/
 s=args.getOpt("SCORE");   
 if (!s.is_empty()) {
    minscore = s.asInt();
    if (minscore>0) GMessage("SCORE=%d\n",  minscore);
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
   
 GString pfile=args.getOpt('f'); //write filtered pairs here
 FILE* fpairs=NULL;
 if (!pfile.is_empty()) {
   if (pfile == "-") 
      fpairs=stdout;
    else 
     if ((fpairs=fopen(pfile, "w"))==NULL)
      GError("Cannot write filtered hits file '%s'!", pfile.chars());
   }
  
 //======== main program loop
 char* line;
 lineno=0;
 while ((line=fgets(inbuf, BUF_LEN-1,inf))!=NULL) {
        int l=strlen(line);
        if (line[l-1]=='\n') line[l-1]='\0';
        GString linecpy(line);
        if (line[0]=='#' || strlen(line)<=1) continue;
        char* tabpos=line;
        //find the 1st tab
        while (*tabpos != '\t' && *tabpos!='\0') tabpos++;
        if (*tabpos=='\0' || tabpos==line)
            GError(ERR_INVALID_PAIR, line);
        *tabpos='\0'; //so line would be the first node name
        if (flt_Exclude && excludeList.hasKey(line)) continue;
        tabpos++; //tabpos is now on the first char of the second field (q_len)
        int score, scov, ovh_r, ovh_l, lcov, pid;
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
          Gswap(q5,q3);
          minus=true;
          }

        //now we should be on the first char of the hitname field
        while (isspace(*tabpos)) tabpos++; //skip any spaces in this second node name field
        if (*tabpos=='\0') GError(ERR_INVALID_PAIR, line);
        //add a string termination after this field
        char* p=tabpos; while (!isspace(*p) && *p!='\0') p++;
        *p='\0';
        //now tabpos contains the exact second sequence string
        if (flt_seqOnly && seqonlyList.Find(line)==NULL && seqonlyList.Find(tabpos)==NULL)
             continue;
        if (flt_seqRestrict && (seqonlyList.Find(line)==NULL || seqonlyList.Find(tabpos)==NULL))
             continue;

        if (strcmp(line, tabpos)==0) {
          GMessage("Warning: self pairing found for node %s\n",line);
          continue;
          }
        if (flt_Exclude && excludeList.hasKey(tabpos)) continue;
        if (!seq_filter(line, tabpos)) continue;
        p++; //move on the first char of the hitlen
        int hitlen=getNextValue(p,line);
        p++;
        int h5=getNextValue(p,line);
        p++;
        int h3=getNextValue(p,line);
        p++;
        if (h5==h3) GError(ERR_INVALID_PAIR, line);
        if (h5>h3) {//never the case with mgblast, actually..
          Gswap(h5,h3);
          minus=!minus;
          }
        pid=getNextValue(p,line); p++;
        score=getNextValue(p,line);
        //compute coverages:
        ovh_r=minus ?(GMIN(q5-1, hitlen-h3)) :(GMIN(hitlen-h3, qlen-q3));
        ovh_l=minus ?(GMIN(h5-1, qlen-q3)) :(GMIN(h5-1, q5-1));
        int overlap = GMAX(q3-q5+1, h3-h5+1);
        if (hitlen>qlen) { //query is shorter
          scov = (int) iround(((double)(q3-q5)*100)/qlen);
          lcov = (int) iround(((double)(h3-h5)*100)/hitlen);
          }
        else {
          lcov = (int) iround(((double)(q3-q5)*100)/qlen);
          scov = (int) iround(((double)(h3-h5)*100)/hitlen);
          }
        lineno++;
        if (scov>=minscov && pid>=minpid && overlap>=minovl
            && score>=minscore && ovh_r <= maxovhang && ovh_l <= maxovhang) {
           //GMessage("%s(%d) %d-%d  | %s(%d) %d-%d, pid=%d%%, score=%d, scov=%d%%, lcov=%d%%\n",
           //          line,qlen,q5,q3,tabpos,hitlen, h5,h3,pid, score, scov, lcov);

           //if (fpairs!=NULL)
           //   fprintf(fpairs, "%s\t%s\t%d\t%d\n",line, tabpos, ovh_l, ovh_r);
           if (addPair(line, tabpos, minus, qlen, q5, q3, hitlen, h5,h3)) {
                  if (fpairs!=NULL) fprintf(fpairs, "%s\n", linecpy.chars());
                  #ifdef _DEBUG
                  if (do_debug) printClusters();
                  #endif
                  }
           }
      } //while hits
  if (inf!=stdin) fclose(inf);
  if (fpairs!=NULL && fpairs!=stdout) fclose(fpairs);
  //==========
  //
  GString outfile=args.getOpt('o');
  if (!outfile.is_empty()) {
     outf=fopen(outfile, "w");
     if (outf==NULL)
        GError("Cannot open file %s for writing!\n",outfile.chars());
     }
   else outf=stdout;
  GMessage("Total t-clusters: %d \n", tclusters.Count());
  tclusters.setSorted(compareCounts);
  if (tclusters.Count()>0) {
     GMessage("Largest cluster has %d nodes\n", tclusters[0]->Count());
     for (int i=0; i<tclusters.Count(); i++) {
        GCluster* c=tclusters[i];
        c->setSorted(compareID);
        if (asfasta) fprintf(outf,">CL%d\t%d\n", i+1, c->Count());
        CNode* n=(CNode*)c->Get(0);fprintf(outf,"%s",n->id);
        for (int j=1; j<c->Count();j++)
          fprintf(outf,"\t%s",c->Get(j)->id);
        fprintf(outf,"\n");
        }
     fflush(outf);
     if (flyt!=NULL) {
      //------- print layouts
       for (int i=0; i<tclusters.Count(); i++) {
          GCluster* c=tclusters[i];
          int clen=0;
          CNode* a=(CNode*)c->Get(0);
          clen=a->ovx3+a->seqlen+a->ovx5;
          fprintf(flyt,">LCL%d %d 1 %d\n", i+1, c->Count(), clen);
          for (int j=0; j<c->Count();j++) {
             CNode* n=(CNode*)c->Get(j);
             char sc;
             int lstart,ltrim, rtrim;
             if (n->nrev) {
               sc='-';lstart=n->ovx3+1;
               ltrim=n->trim3;
               rtrim=n->trim5;
               }  else {
               sc='+';lstart=n->ovx5+1;
               ltrim=n->trim5;
               rtrim=n->trim3;
               }
             fprintf(flyt, "%s %c %d %d %d %d\n",
                      n->id,sc,n->seqlen, lstart, ltrim, rtrim);
             } //for each CNode
          } //for each cluster
      }// if layout file    
     }// clusters found
  if (outf!=stdout) fclose(outf);
  if (flyt!=NULL && flyt!=stdout) fclose(flyt);
  all.Clear(); //this frees the nodes
  tclusters.Clear(); /* this frees the clusters,
     but should not try to free the nodes inside */
  //GMessage("the tclusters list was cleared!\n");
  //GMessage("*** all done ***\n");
  #ifdef __WIN32__
   getc(stdin);
  #endif
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
 return (int)iround(d);
 //backup
 }
//=====
int compareCounts(void* p1, void* p2) {
 int c1=((GCluster*)p1)->Count();
 int c2=((GCluster*)p2)->Count();
 return (c1>c2)?-1:((c1<c2)?1:0);
}

//=====
int compareID(void* p1, void* p2) {
 return strcmp(((CNode*)p1)->id,((CNode*)p2)->id);
}


//-----
//both nodes are new - create global entries for them
// and a new shared t-cluster containing them
void newPair(char* id1, char* id2, CNode*& q, CNode*& h, bool minus,
             int qlen, int q5, int q3, int hlen, int h5, int h3) {
 q=new CNode(id1); /* this creates a new cluster with only self in it */
 h=new CNode(id2, q); //h also added to q->cluster

 all.shkAdd(q->id,q);
 q->seqlen=qlen;
 h->seqlen=hlen;
 all.shkAdd(h->id,h);
 q->trim5=q5-1; q->trim3=qlen-q3;
 h->trim5=h5-1; h->trim3=hlen-h3;

 // now we have to update ovhang boundaries:
 if (minus) { //reverse complement alignment
  // relative offset to be made positive:
  int qq5=qlen-q3+1;
  if (h5<qq5) { // q is leftmost (offset 0)
     //h->nofs=qq5-h5;
     //h->cluster->len=GMAX((hlen+qq5-h5),qlen);
     /*h->nto=q; */
     h->nrev=true; //h is considered "reversed"
     // left extension for h
     h->ovx5=qq5-h5;
     }
   else { // h is the leftmost (offset 0)
     //q->nofs=h5-qq5;
     //h->cluster->len=GMAX((qlen+h5-qq5),hlen);
     /*q->nto=h;   */
     q->nrev=true; //q is "reversed" in this alignment
     //extension:
     q->ovx3=h5-qq5;
     }
   /*if (h5>qq5) q->ovx3=h5-qq5;
          else   h->ovx5=qq5-h5; */
   int hh3=hlen-h3+1;
   if (hh3>q5)
               q->ovx5=hh3-q5;
          else h->ovx3=q5-hh3;
   }
 else {  //forward strand alignment
  // relative offset to be made positive:
  if (h5<q5) { // q=q is leftmost
     //h->nofs=q5-h5;
     //h->cluster->len=GMAX((hlen+q5-h5),qlen);
     //h->nto=q;
     h->ovx5 = q5-h5;
     }
   else { // h=h is the leftmost
     //q->nofs=h5-q5;
     //h->cluster->len=GMAX((qlen+h5-q5),hlen);
     //q->nto=h;
     q->ovx5 = h5-q5;
     }
 /* if (h5>q5) q->ovx5 = h5-q5;
        else h->ovx5 = q5-h5;*/
  int hh3=hlen-h3+1;
  int qq3=qlen-q3+1;
  if (hh3>qq3)
               q->ovx3 = hh3-qq3;
          else h->ovx3 = qq3-hh3;
  }
 //one cluster created: add it to the global list of t-clusters:
 #ifdef _DEBUG
 if (do_debug) {
   CNode* a=q;CNode* b=h;
   if (strcmp(a->id,b->id)>0) { a=h;b=q; }
   char sa=a->nrev?'-':'+';
   char sb=b->nrev?'-':'+';
   fprintf(stderr, ">[%d]NewPair: [%c%s %d,%d ~ %c%s %d,%d]\n",
              lineno, sa, a->id, a->ovx5, a->ovx3, sb, b->id, b->ovx5, b->ovx3);
   }
 #endif  
 tclusters.Add(q->cluster);
}

//add a new node q to be paired with an existing node h
// which is already in a cluster
// must check if ovx5/ovx3 allow this
CNode* addNode(char* qname, CNode* h, bool minus, int qlen, int q5, int q3, int hlen, int h5, int h3) {
 //virtual h coverage to limit the virtual ovhang
 int hvcov = (int) iround(((double)(h3-h5)*100)/(hlen+h->ovx3+h->ovx5));
 int misovh5, misovh3, //mismatched overhangs at 5' and 3' of q
     d5, d3;
 bool q_rev=false;
 if (minus) { //reverse complement alignment
  //misovh5=GMIN(q5-1,    hlen+h->ovx3-h3);
  //misovh3=GMIN(qlen-q3, h->ovx5+h5-1);
  misovh5=q5-1;
  misovh3=qlen-q3;
  d5= misovh5 - (hlen+h->ovx3-h3); //offset of q relative to h at 5' end of q
  d3=misovh3 - (h->ovx5+h5-1); //offset of q relative to h at 3' end of q
  if (d5>0) misovh5 -=d5;//GMIN(q5-1, hlen+h->ovx3-h3);
  if (d3>0) misovh3 -=d3;//GMIN(qlen-q3, h->ovx5+h5-1);
  if (hvcov<minscov || misovh5>maxovhang || misovh3>maxovhang)
     return NULL;
  if (!h->nrev) q_rev=true;
  } // was reverse
 else {  //forward alignment
  //misovh5=GMIN(q5-1,    h5+h->ovx5-1);
  //misovh3=GMIN(qlen-q3, hlen+h->ovx3-h3);
  misovh5=q5-1;
  misovh3=qlen-q3;
  d5= misovh5 - (h5+h->ovx5-1); //offset of q relative to h at 5' end of q
  d3= misovh3 - (hlen+h->ovx3-h3); //offset of q relative to h at 3' end of q
  if (d5>0) misovh5 -=d5;
  if (d3>0) misovh3 -=d3;
  if (hvcov<minscov || misovh5>maxovhang || misovh3>maxovhang)
     return NULL;
  if (h->nrev) q_rev=true;
  }
 // it passed the virtual overhang tests, store it
 CNode* q=new CNode(qname,h); //this automatically adds q to the cluster of h
 q->trim5=q5-1; q->trim3=qlen-q3; //virtual trim to the ends of the alignment
 q->nrev=q_rev;
 q->ovx3=-d3;
 q->ovx5=-d5;
 if (d5>0)
   for (int i=0;i<h->cluster->Count();i++) {
     CNode* n=(CNode*) h->cluster->Get(i);
     if (n->nrev==q->nrev) n->ovx5+=d5;
                      else n->ovx3+=d5;
     }
 if (d3>0)
   for (int i=0;i<h->cluster->Count();i++) {
     CNode* n=(CNode*)h->cluster->Get(i);
     if (n->nrev==q->nrev) n->ovx3+=d3;
                      else n->ovx5+=d3;
     }

 //positioning is relative to h as reference
 //q->nto=h;
 //adjust h trims only if shorter..
 if (h->trim5>h5-1) h->trim5=h5-1;
 if (h->trim3>hlen-h3) h->trim3=hlen-h3;
 all.shkAdd(q->id,q);
 q->seqlen=qlen;
 #ifdef _DEBUG
 if (do_debug) {
   CNode* a=q;CNode* b=h;
   if (strcmp(a->id,b->id)>0) { a=h;b=q; }
   char sa=a->nrev?'-':'+';
   char sb=b->nrev?'-':'+';
   fprintf(stderr, ">[%d]addNode: [%c%s %d,%d ~ %c%s %d,%d]\n",
              lineno, sa, a->id, a->ovx5, a->ovx3, sb, b->id, b->ovx5, b->ovx3);
   }
 #endif  
 return q;
}
//override for readClones compatibility!
CNode* addClone(char* qname, CNode* h) {
 //use the same cluster
 CNode* q=new CNode(qname,h); //this automatically adds q to the cluster of h
 all.shkAdd(q->id,q);
 q->seqlen=600; //just for consistency..
 return q;
}


//create a link between two existing nodes..
// if they don't violate the virtual overhangs!
bool linkNodes(CNode* q, CNode* h, bool minus,
             int qlen, int q5, int q3, int hlen, int h5, int h3) {
 if (q->cluster==h->cluster) return false;
 int hvcov = (int) iround(((double)(h3-h5)*100)/(hlen+h->ovx3+h->ovx5));
 int qvcov = (int) iround(((double)(q3-q5)*100)/(qlen+q->ovx3+q->ovx5));
 int misovh5, misovh3, //mismatched overhangs at 5' and 3' of q
     d5, d3;
 bool q_rev=false; //indicates if one of the clusters should be reversed
 if (minus) { //reverse complement alignmnet
  misovh5=q->ovx5+q5-1;
  misovh3=qlen+q->ovx3-q3;
  d5 = misovh5 - (hlen+h->ovx3-h3); //offset of q's cluster relative to
                                 //h's cluster at 5' end of q's cluster
  d3 = misovh3 - (h->ovx5+h5-1); //offset of q's cluster relative to
                                 //h's cluster at 3' end of q's cluster
  if (d5>0) misovh5 -=d5;//GMIN(q5-1, hlen+h->ovx3-h3);
  if (d3>0) misovh3 -=d3;//GMIN(qlen-q3, h->ovx5+h5-1);
  if (hvcov<minscov || qvcov<minscov || misovh5>maxovhang || misovh3>maxovhang)
     return false;
  if (h->nrev==q->nrev) q_rev=true;
   }
 else { //forward strand alignment
  misovh5=q5+q->ovx5-1;
  misovh3=qlen+q->ovx3-q3;
  d5= misovh5 - (h5+h->ovx5-1); //offset of q->cluster relative to h->cluster at 5' end of q
  d3= misovh3 - (hlen+h->ovx3-h3); //offset of q->cluster relative to h at 3' end of q
  if (d5>0) misovh5 -=d5;
  if (d3>0) misovh3 -=d3;
  if (hvcov<minscov || qvcov<minscov || misovh5>maxovhang || misovh3>maxovhang)
     return false;
  if (h->nrev!=q->nrev) q_rev=true;
  }
 // it passed the virtual overhang tests, merge clusters

if (q->trim5>q5-1) q->trim5=q5-1;
if (q->trim3>qlen-q3) q->trim3=qlen-q3;
if (h->trim5>h5-1) h->trim5=h5-1;
if (h->trim3>hlen-h3) h->trim3=hlen-h3;

//we have to merge the two clusters
// MERGE q->cluster into h->cluster ! q=>h
//-- adjust all ovx and reverse flags
//   for q->cluster
bool prevqrev=q->nrev; //current state of q
for (int i=0;i<q->cluster->Count();i++) {
 CNode* qn=(CNode*)q->cluster->Get(i);
 if (qn->nrev==prevqrev) { //same orientation with q
    qn->ovx3+=-d3;
    qn->ovx5+=-d5;
    }
   else { //rev.compl. relative to q
    qn->ovx5+=-d3;
    qn->ovx3+=-d5;
    }
 if (q_rev) qn->nrev=!qn->nrev;
}

//join the two clusters (add q->cluster to h->cluster)
 GCluster* qcl=q->cluster;
 for (int i=0; i<qcl->Count(); i++) {
     CNode* n= (CNode*)qcl->Get(i);
     n->cluster=h->cluster;
     h->cluster->Add(n);  /* check duplicates? they should never be there! */
     //qcl->Forget(i);
     }
//-- adjust ovx
 if (d5>0)
   for (int i=0;i<h->cluster->Count();i++) {
     CNode* n=(CNode*)h->cluster->Get(i);
     if (n->nrev==q->nrev) n->ovx5+=d5;
                      else n->ovx3+=d5;
     }
 if (d3>0)
   for (int i=0;i<h->cluster->Count();i++) {
     CNode* n=(CNode*)h->cluster->Get(i);
     if (n->nrev==q->nrev) n->ovx3+=d3;
                      else n->ovx5+=d3;

     }
 tclusters.Remove(qcl);
 #ifdef _DEBUG
 if (do_debug) {
   CNode* a=q;CNode* b=h;
   if (strcmp(a->id,b->id)>0) { a=h;b=q; }
   char sa=a->nrev?'-':'+';
   char sb=b->nrev?'-':'+';
   char sc=q_rev? '-':'+';
   fprintf(stderr, ">[%d]linkNodes: [%c%s %d,%d ~ %c%s %d,%d] %d %d %c\n",
              lineno, sa, a->id, a->ovx5, a->ovx3, sb, b->id, b->ovx5, b->ovx3,
                       d5, d3, sc);
   }
 #endif  
 return true;
}


bool isET(char* name) {
 return (startsWith(name, "np|") ||
        startsWith(name, "et|") ||
        startsWith(name, "egad|") ||
        startsWith(name, "preegad|")
        );
 }

bool seq_filter(char* s1, char* s2) {
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
 }
 
 
/* addPair() should slowly build the t-clusters and update
   the links for both n1 and n2
*/
bool addPair(char* s1, char* s2, bool minus, 
             int qlen, int q5, int q3, int hlen, int h5, int h3) {
 CNode* n1=all.Find(s1);
 CNode* n2=all.Find(s2); 
 //GMessage("Add pair %s %s\n",n1,n2);
 if (n1==NULL) {
   //###### n1 is new
   if (n2==NULL) {//n2 is also new 
        newPair(s1,s2, n1, n2, minus, qlen, q5, q3, hlen, h5, h3);
        return true;
        }
        //this should set d1 and d2 to the each new CNode*
      else {//n1 is new, but n2 already is in a cluster
        // check n2's ovh boundaries, update n2 and n1's ovh boundaries
        n1 = addNode(s1,n2, minus, qlen, q5, q3, hlen, h5, h3);
        return (n1!=NULL);
        }
   }
  else {
   //###### n1 is in a cluster already
   if (n2==NULL) {//n2 is new
        // check n1's ovh boundaries, update n1 and n2 boundary
        n2=addNode(s2,n1, minus, hlen, h5, h3, qlen, q5, q3);
        return (n2!=NULL);
        }
      else {//n2 also has a t-cluster, merge them
        // check & update n1 and n2's ovh boundaries
        return linkNodes(n1,n2, minus, qlen, q5, q3, hlen, h5, h3);
        }
   }
}

void showCl(GCluster* C) {
  GMessage("{%s",C->Get(0)->id);
  for (int i=1;i<C->Count();i++) 
     GMessage(" %s",C->Get(i)->id);
  GMessage("}");
}


#ifdef DBGTRACE
void memlog(int level, const char* msg1, 
         GCluster* C=NULL, 
         const char* msg2=NULL) {
 GString s;
 if (level>0) {
   s.format("%d>", level);
   s.padR(level+2);
   }
 if (msg1!=NULL) s+=msg1;
 if (C!=NULL) {
     GString s2;
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

int readClones(FILE* f) {
//expects space delimited reads, one clone per line
  char* line;
  int count=0;
  while ((line=fgets(inbuf, BUF_LEN-1,f))!=NULL) {
   int l=strlen(line);
   count++;
   if (line[l-1]=='\n') line[l-1]='\0';
   GString linestr(line);
   linestr.startTokenize("\t ");
   GString token;
   CNode* head=NULL;
   while (linestr.nextToken(token)) {
       if (flt_seqRestrict && !seqonlyList.Find(token.chars()))
          continue;
       if (head==NULL) {
         head=new CNode((char *)token.chars());
         tclusters.Add(head->cluster); 
         all.shkAdd(head->id,head);
         }
        else 
         addClone((char *)token.chars(), head);      
     }
   }
 return count;
}
