#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "GArgs.h"
#include "GStr.h"
#include "GHash.hh"
#include "GList.hh"
#include "GCdbYank.h"
#define USAGE "Usage:\n\
 mblasm <mgblast_sortedhits> [-d <fastadb.cidx>] [-R]\n\
   [-o <outfile.ace>]\n\
   <mgblast_sortedhits> is a mgblast -D5 tabulated stream \n\
                        sorted by alignment score\n\
   -d multi-fasta database index (cdbfasta) for all the sequences\n\
   -o the output acefile is written to <outfile.ace> instead of stdout\n\
   -R assumes the largest sequence in the first hit is also the \n\
      reference sequence and that all the other sequences align to it\n"

#define BUF_LEN 8192*4 //maximum input line length
#define ERR_INVALID_PAIR "Invalid input line encountered:\n%s\n"

/*---
 OVHCOV=xx is the maximum percent length of an overhang \n\
           (maximum overhang length as percentage of participant sequence length)\n\
*/

//============ structures:
//define a cluster object as a list of node names (GStrs)

// dumb node sorting function -- by pointer value (!)
int comparePtr(void* p1, void* p2) {
  return (p1>p2) ? -1 : ((p1<p2)?1:0);
  }

class CNode;

class GCluster :public GList<CNode> {
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
  GCluster():GList<CNode>(true,false,true) {
   //default is: sorted by Node pointer, don't free nodes, unique
    }
  GCluster(bool sorted, bool free_elements=true, bool beUnique=false):
     GList<CNode>(sorted,free_elements,beUnique) {
       }
  };

class CNode {
 public:
   char* id; //actual node ID
   int l_ovx; // left overhang extension
   int r_ovx; // right overhang extension
   int l_trim; // left overhang extension
   int r_trim; // right overhang extension
   // positioning:
   CNode* nto; // node which the offset nofs is relative to
   bool nrev; // reverse complement positioning to nto?
   int nofs; // relative offset to nto?
   //
   GCluster* cluster; /* pointer to a list of nodes (the t-cluster
                         containing this sequence) */
   //GCluster links; //list of all nodes linked to this node
   CNode(char* name) {
    //links are sorted, not free, unique
    id=Gstrdup(name);
    l_ovx=0;
    r_ovx=0;
    l_trim=0;
    r_trim=0;
    nofs=0;
    nto=NULL;
    nrev=false;
    cluster=new GCluster;
    cluster->Add(this);
    }
   // create & add to an existing cluster 
   CNode(char* name, CNode* pairnode) {
    id=Gstrdup(name);
    l_ovx=0;
    r_ovx=0;
    l_trim=0;
    r_trim=0;
    nofs=0;
    nto=NULL;
    nrev=false;
    cluster=pairnode->cluster; //same cluster
    cluster->Add(this); //add self to the cluster
    }
   /* no need for these for basic transitive closure:
   CNode(char* name):links(true,false,true) {
    //links are sorted, not free, unique
    id=Gstrdup(name);    
    cluster=new GCluster;
    cluster->Add(this);
    }
   CNode(char* name, CNode* pairnode):links(true,false,true) {
    id=Gstrdup(name);
    cluster=pairnode->cluster; //same cluster
    cluster->Add(this); //add self to the cluster
    }
   // links management:
   int linksCount(GCluster* C) {
    // returns the number of links that this node  has within cluster C    
    int cnt=0;
    for (int i=0;i<links.Count();i++) 
       if (C->IndexOf(links[i])>=0) cnt++;
    return cnt;
    }
   */ 
  //comparison: just pointer based:
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

FILE* outf;

//all clusters are stored in a unique sorted list
GHash<int> excludeList;
GHash<int> seqonlyList;
GList<GCluster> tclusters(true,true,true);

//=========== functions used:
void addPair(char* s1, char* s2, bool minus, int qlen, int q5, int q3, int hlen, int h5, int h3);
int compareCounts(void* p1, void* p2);
int compareID(void* p1, void* p2); //compare node IDs
int getNextValue(char*& p, char* line);
bool seq_filter(char* s1, char* s2); //returns true if the pair can pass through the filter
int readNames(FILE* f, GHash<int>& xhash);
int readClones(FILE* f);

//========================================================
//====================     main      =====================
//========================================================
int main(int argc, char * const argv[]) {
 GArgs args(argc, argv, "hHf:c:o:s:r:x:SEQFLT=PID=COV=OVHANG=OVL=SCORE=DEBUG=");
 int e;
 if ((e=args.isError())>0)
    GError("%s\nInvalid argument: %s\n", USAGE, argv[e]);
 if (args.getOpt('h')!=NULL) GError("%s\n", USAGE);

 GStr infile;
 if (args.startNonOpt()) {
        infile=args.nextNonOpt();
        GMessage("Given file: %s\n",infile.chars());
        }
 asfasta=(args.getOpt('H')==NULL);
 do_debug=(args.getOpt("DEBUG")!=NULL);
 GStr s=args.getOpt('x');
 FILE* fxclude=NULL;
 if (!s.is_empty()) {
   if ((fxclude=fopen(s, "r"))==NULL)
      GError("Cannot open exclusion file '%s'!\n", s.chars());
   int c=readNames(fxclude, excludeList);
   GMessage("Loaded %d sequences to be excluded.\n", c);
   fclose(fxclude);
   flt_Exclude=(c>0);
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
   
 GStr pfile=args.getOpt('f'); //write filtered pairs here
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
 while ((line=fgets(inbuf, BUF_LEN-1,inf))!=NULL) {
        int l=strlen(line);
        if (line[l-1]=='\n') line[l-1]='\0';
        GStr linecpy(line);
        if (strlen(line)<=1) continue;
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
        if (h5>h3) {
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
          scov = (int) floor(((double)(q3-q5)*100)/qlen);
          lcov = (int) floor(((double)(h3-h5)*100)/hitlen);
          }
        else {
          lcov = (int) floor(((double)(q3-q5)*100)/qlen);
          scov = (int) floor(((double)(h3-h5)*100)/hitlen);
          }
          
        if (scov>=minscov && pid>=minpid && overlap>=minovl
            && score>=minscore && ovh_r <= maxovhang && ovh_l <= maxovhang) {
           //GMessage("%s(%d) %d-%d  | %s(%d) %d-%d, pid=%d%%, score=%d, scov=%d%%, lcov=%d%%\n",
           //          line,qlen,q5,q3,tabpos,hitlen, h5,h3,pid, score, scov, lcov);
           if (fpairs!=NULL) fprintf(fpairs, "%s\n", linecpy.chars());
           //if (fpairs!=NULL) 
           //   fprintf(fpairs, "%s\t%s\t%d\t%d\n",line, tabpos, ovh_l, ovh_r);
           addPair(line, tabpos, minus, qlen, q5, q3, hitlen, h5,h3);
           }
      } //while hits
  if (inf!=stdin) fclose(inf);
  if (fpairs!=NULL && fpairs!=stdout) fclose(fpairs);
  //==========
  //
  GStr outfile=args.getOpt('o');
  if (!outfile.is_empty()) {
     outf=fopen(outfile, "w");
     if (outf==NULL)
        GError("Cannot open file %s for writing!\n",outfile.chars());
     }
   else outf=stdout;
  GMessage("Total t-clusters: %d \n", tclusters.Count());
  tclusters.setSorted(&compareCounts);
  if (tclusters.Count()>0) {
     GMessage("Largest cluster has %d nodes\n", tclusters[0]->Count());
     for (int i=0; i<tclusters.Count(); i++) {
        GCluster* c=tclusters[i];
        c->setSorted(&compareID);
        if (asfasta) fprintf(outf,">CL%d\t%d\n", i+1, c->Count());
        CNode* n=c->Get(0);fprintf(outf,"%s",n->id);
        for (int j=1; j<c->Count();j++)
          fprintf(outf,"\t%s",c->Get(j)->id);
        fprintf(outf,"\n");
        }
     fflush(outf);
     }
  all.Clear(); //this frees the nodes
  tclusters.Clear(); /* this frees the clusters, 
     but should not try to free the nodes inside */
  //GMessage("the tclusters list was cleared!\n");
  if (outf!=stdout) fclose(outf);
  GMessage("*** all done ***\n");
  //getc(stdin);
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
 return (int)floor(d);
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
 h=new CNode(id2, q); //h->cluster set to q->cluster 
 all.shkAdd(q->id,q);
 all.shkAdd(h->id,h);
 q->l_trim=q5-1; q->r_trim=qlen-q3;
 h->l_trim=h5-1; h->r_trim=hlen-h3;

 // now we have to update ovhang boundaries:
 if (minus) { //reverse complement alignmnet
  // relative offset to be made positive:
  int qq5=qlen-q3+1;
  if (h5<qq5) { // q=q is leftmost
     h->nofs=qq5-h5;
     h->nto=q;
     h->nrev=true;
     }
   else { // h=h is the leftmost
     q->nofs=h5-qq5;
     q->nto=h;
     q->nrev=true;
     }
   if (h5>qq5) q->r_ovx=h5-qq5;
          else h->l_ovx=qq5-h5;
   int hh3=hlen-h3+1;
   if (hh3>q5) q->l_ovx=hh3-q5;
          else h->r_ovx=q5-hh3;
   }
 else {  //forward strand alignment
  // relative offset to be made positive:
  if (h5<q5) { // q=q is leftmost
     h->nofs=q5-h5;
     h->nto=q;
     h->nrev=false;
     }
   else { // h=h is the leftmost
     q->nofs=h5-q5;
     q->nto=h;
     q->nrev=false;
     }
  if (h5>q5) q->l_ovx = h5-q5;
        else h->l_ovx = q5-h5;
  int hh3=hlen-h3+1;
  int qq3=qlen-q3+1;
  if (hh3>qq3) q->r_ovx = hh3-qq3;
          else h->r_ovx = qq3-hh3;
  }
 //one cluster created: add it to the global list of t-clusters:
 tclusters.Add(q->cluster);
}

//add a new node q to be paired with an existing node h
// must check if l_ovx/r_ovx allow this
CNode* addNode(char* qname, CNode* h, bool minus, int qlen, int q5, int q3, int hlen, int h5, int h3) {
 int q_r_ovx, q_l_ovx, h_l_ovx, h_r_ovx; //putative adjustments
 q_r_ovx=0;q_l_ovx=0; h_l_ovx=0;h_r_ovx=0;
 //virtual h coverage to limit the virtual ovhang
 int hvcov = (int) floor(((double)(h3-h5)*100)/(hlen+h->r_ovx+h->l_ovx));
 int ovh_l, ovh_r;
 if (minus) { //reverse complement alignmnet
  ovh_r=GMIN(q5-1,    hlen+h->r_ovx-h3);
  ovh_l=GMIN(qlen-q3, h->l_ovx+h5-1);
  if (hvcov<minscov || ovh_r>maxovhang || ovh_l>maxovhang) 
     return NULL;  
  int qq5=qlen-q3+1;
  if (h5>qq5) q_r_ovx=h5-qq5;
          else h_l_ovx=qq5-h5;
   int hh3=hlen-h3+1;
   if (hh3>q5) q_l_ovx=hh3-q5;
          else h_r_ovx=q5-hh3;
   
   }
 else {  //forward strand alignment
  ovh_r=GMIN(qlen-q3, hlen+h->r_ovx-h3);
  ovh_l=GMIN(q5-1,    h5+h->l_ovx-1);
  if (hvcov<minscov || ovh_r>maxovhang || ovh_l>maxovhang) 
     return NULL;  
  if (h5>q5) q_l_ovx = h5-q5;
        else h_l_ovx = q5-h5;
  int hh3=hlen-h3+1;
  int qq3=qlen-q3+1;
  if (hh3>qq3) q_r_ovx = hh3-qq3;
          else h_r_ovx = qq3-hh3;
  
  }
 // it passed the virtual overhang tests, store it
 CNode* q=new CNode(qname,h); //this automatically adds q to the cluster of h
 q->l_trim=q5-1; q->r_trim=qlen-q3;
 q->r_ovx=q_r_ovx;
 q->l_ovx=q_l_ovx;
 //positioning is relative to h as reference
 q->nto=h;
 q->nrev = minus;
 q->nofs = minus ? h5-(qlen-q3+1) : h5-q5;
 //adjust h trims only if shorter:
 if (h->l_trim>h5-1) h->l_trim=h5-1; 
 if (h->r_trim>hlen-h3) h->r_trim=hlen-h3;
 if (h->l_ovx<h_l_ovx) h->l_ovx=h_l_ovx;
 if (h->r_ovx<h_r_ovx) h->r_ovx=h_r_ovx;
 all.shkAdd(q->id,q);
 return q;
}
//override for readClones compatibility!
CNode* addClone(char* qname, CNode* h) {
 //use the same cluster
 CNode* q=new CNode(qname,h); //this automatically adds q to the cluster of h
 all.shkAdd(q->id,q);
 return q;
}


//create a link between two existing nodes..
// if they don't violate the virtual overhangs!
CNode* linkNodes(CNode* q, CNode* h, bool minus, 
             int qlen, int q5, int q3, int hlen, int h5, int h3) {
 if (q->cluster==h->cluster) return NULL;
 int q_r_ovx=0, q_l_ovx=0, h_l_ovx=0, h_r_ovx=0; //putative adjustments
 int hvcov = (int) floor(((double)(h3-h5)*100)/(hlen+h->r_ovx+h->l_ovx));
 int qvcov = (int) floor(((double)(q3-q5)*100)/(qlen+q->r_ovx+q->l_ovx));
 int ovh_l, ovh_r;
 if (minus) { //reverse complement alignmnet
  ovh_r=GMIN(q->l_ovx+q5-1,    hlen+h->r_ovx-h3);
  ovh_l=GMIN(qlen+q->r_ovx-q3, h->l_ovx+h5-1);
  if (hvcov<minscov || qvcov<minscov || ovh_r>maxovhang || ovh_l>maxovhang) 
     return NULL;
  int qq5=qlen-q3+1;
  if (h5>qq5) q_r_ovx=h5-qq5;
          else h_l_ovx=qq5-h5;
   int hh3=hlen-h3+1;
   if (hh3>q5) q_l_ovx=hh3-q5;
          else h_r_ovx=q5-hh3;
   
   }
 else {  //forward strand alignment
  ovh_r=GMIN(qlen+q->r_ovx-q3, hlen+h->r_ovx-h3);
  ovh_l=GMIN(q5+q->l_ovx-1,    h5+h->l_ovx-1);
  if (hvcov<minscov || qvcov<minscov || ovh_r>maxovhang || ovh_l>maxovhang) 
     return NULL;  
  if (h5>q5) q_l_ovx = h5-q5;
        else h_l_ovx = q5-h5;
  int hh3=hlen-h3+1;
  int qq3=qlen-q3+1;
  if (hh3>qq3) q_r_ovx = hh3-qq3;
          else h_r_ovx = qq3-hh3;
  
  }
 // it passed the virtual overhang tests, store it

if (q->l_trim>q5-1) q->l_trim=q5-1; 
if (q->r_trim>qlen-q3) q->r_trim=qlen-q3;
if (q->l_ovx<q_l_ovx) q->l_ovx=q_l_ovx;
if (q->r_ovx<q_r_ovx) q->r_ovx=q_r_ovx;

if (h->l_trim>h5-1) h->l_trim=h5-1; 
if (h->r_trim>hlen-h3) h->r_trim=hlen-h3;
if (h->l_ovx<h_l_ovx) h->l_ovx=h_l_ovx;
if (h->r_ovx<h_r_ovx) h->r_ovx=h_r_ovx;


//we have to merge the two clusters
//add to the bigger cluster the sequences from the smaller one
//(big fish swallows the smaller one)
GCluster* src, *dest;
CNode* nbig;
if (q->cluster->Count()>h->cluster->Count()) {
  // q is with the bigger fish
  // h aligns to q
  nbig=q;
  h->nrev=minus;
  h->nto=q;
  h->nofs=minus ?  qlen-q3+1-h5 : q5-h5;
  dest = q->cluster;
  src  = h->cluster;
  }
else {
  nbig=h;
  q->nrev=minus;
  q->nto=h;
  q->nofs= minus ? h5-(qlen-q3+1) : h5-q5;
  dest = h->cluster;
  src  = q->cluster;
  }
 //merge the clusters (in the bigger one)
 //and remove the smaller cluster  
 for (register int i=0; i<src->Count(); i++) {
     CNode* n= src->Get(i);
     dest->Add(n);  /* check duplicates? they should never be there!
          if (dest->Add(n)<0)  {
                   delete n;
                   src->Put(i,NULL);
                   } */
     //each node from smaller cluster should now point to the bigger (current)cluster
     n->cluster=dest;
     }
 src->setFreeItem(false);
 tclusters.Remove(src);
 return nbig;
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
void addPair(char* s1, char* s2, bool minus, 
             int qlen, int q5, int q3, int hlen, int h5, int h3) {
 CNode* n1=all.Find(s1);
 CNode* n2=all.Find(s2); 
 //GMessage("Add pair %s %s\n",n1,n2);
 if (n1==NULL) {
   //###### n1 is new
   if (n2==NULL) //n2 is also new 
        newPair(s1,s2, n1, n2, minus, qlen, q5, q3, hlen, h5, h3);
          //this should set d1 and d2 to the each new CNode*
      else {//n1 is new, but n2 already is in a cluster
        // check n2's ovh boundaries, update n2 and n1's ovh boundaries
        n1 = addNode(s1,n2, minus, qlen, q5, q3, hlen, h5, h3);
        }
   }
  else {
   //###### n1 is in a cluster already
   if (n2==NULL) {//n2 is new
        // check n1's ovh boundaries, update n1 and n2 boundary
        n2=addNode(s2,n1, minus, hlen, h5, h3, qlen, q5, q3);
        }
      else {//n2 also has a t-cluster, merge them
        // check & update n1 and n2's ovh boundaries
        linkNodes(n1,n2, minus, qlen, q5, q3, hlen, h5, h3);
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

int readClones(FILE* f) {
//expects space delimited reads, one clone per line
  char* line;
  int count=0;
  while ((line=fgets(inbuf, BUF_LEN-1,f))!=NULL) {
   int l=strlen(line);
   count++;
   if (line[l-1]=='\n') line[l-1]='\0';
   GStr linestr(line);
   linestr.startTokenize("\t ");
   GStr token;
   CNode* head=NULL;
   while (linestr.nextToken(token)) {
       if (flt_seqRestrict && !seqonlyList.Find(token.chars()))
          continue;
       if (head==NULL) {
         head=new CNode((char *)token.chars());
         tclusters.Add(head->cluster); 
         all.shkAdd(head->id,head);
         /*if (do_debug)
           GMessage("Added head node '%s'\n",head->id);*/
         }
        else 
         addClone((char *)token.chars(), head);      
     }
   }
 return count;
}
