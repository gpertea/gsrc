#include "GBase.h"
#include "GArgs.h"
#include "GStr.h"
#include "GList.hh"
#include "GHash.hh"
//#include "gffmrna.h"
#include "gff.h"
#include <ctype.h>


#define USAGE "Usage:\n\
   wigcov [-f <input.gff>] [-A] <input.wig>\n\
   \n\
   Shows the coverage islands found in <input.wig>, or,\n\
   if the -f option is given, acts as a filter for the given transcripts,\n\
   only showing those transcripts overlapping the coverage islands. \n\
   \n\
   -A option shows such a transcript even if a single exon overlaps the\n\
    coverage islands (default: all exons should overlap)\n\
    "
class CovIsland {
  uint nextp;
 public:
  uint start;
  uint end;
  CovIsland(uint gstart, int initcov, int cspan=1) {
    init(gstart, initcov, cspan);
    nextp=end+1;
    }
  CovIsland(GLineReader& covreader, GStr& seqname, uint maxgap=10);
     //construct a whole island by reading next lines from file
     //and stopping if a 0-coverage gap of at least maxgap is encountered
  void init(uint gstart, int initcov, int cspan=1) {
    start=gstart;
    end=gstart+cspan-1;
    nextp=end+1;
    }
  void extend(int scov, int cspan=1) {
    if (scov>255) scov=255;
    end+=cspan;
    nextp=end+1;
    }
  bool operator==(CovIsland& c){
    return (start==c.start && this->end==c.end);
    }
  bool operator<(CovIsland& c){
    return (start==c.start) ? (end<c.end) : (start<c.start);
    }
};

CovIsland::CovIsland(GLineReader& covreader, GStr& gseqname, uint maxgap) {
  //build a coverage island by parsing the next .wig lines
 static const char msg_ERR_COVPARSE[] = " Error at parsing coverage file, line: %s\n";
 char *line;
 char* fields[5];
 start=0;
 end=0;
 nextp=0;
 gseqname.clear(); //it will update it
 while ((line=covreader.nextLine())!=NULL) {
  GStr l=line; //keep a copy for error messages
  if (l.length()<3) continue;
  int tidx=strsplit(line,fields, 4, '\t');
  if (tidx<4) continue;
  uint cstart=0, cend=0;
  int lcov=0;
  if (!gseqname.is_empty() && gseqname!=fields[0]) {
     //different genomic seq, return here
     nextp=MAXUINT;
     covreader.pushBack(); //leave it for next CovIsland
     break;
     }

  if (!(parseUInt(fields[1],cstart) && parseUInt(fields[2],cend) && parseInt(fields[3],lcov)))
            GError(msg_ERR_COVPARSE, l.chars());
  cend--;
  if (lcov==0) { // no coverage segment
      if (cend-cstart>maxgap) { //gap larger than minimum intron
            if (end) { nextp=cend+1; break; } //close this CovIsland
                else continue; //skip initially uncovered region
            }
      }// cov gap
  if (end) {
    if (cstart>end+1) { //implied gap
       if (cstart-end-1>maxgap) {
           nextp=cstart;
           covreader.pushBack(); //read it next time
           break; //island end
           }
       extend(0, cstart-end+1); //extend through small gap
       } //implied gap
    extend(lcov, cend-cstart+1);
    }
   else { //new start
    init(cstart, lcov, cend-cstart+1);
    if (gseqname.is_empty()) gseqname=fields[0];
    }
  }//while coverage lines
 // -- we just finished putting together a coverage island

}

class GLocus:public GSeg {
 public:
   int gseq_id; //id of underlying genomic sequence
   GList<GffObj> mrnas; //list of transcripts for this locus
   GArray<GSeg> mexons; //list of merged exons in this region
   GLocus(GffObj* mrna=NULL):mrnas(true,false,false), mexons(true,true) {
     //this will NOT free mrnas!
     gseq_id=-1;
     if (mrna!=NULL) {
        start=mrna->exons.First()->start;
        end=mrna->exons.Last()->end;;
        gseq_id=mrna->gseq_id;
        GSeg seg;
        for (int i=0;i<mrna->exons.Count();i++) {
          seg.start=mrna->exons[i]->start;
          seg.end=mrna->exons[i]->end;
          mexons.Add(seg);
          }
        mrnas.Add(mrna);
        mrna->uptr=this;
        }
     }
  void addMerge(GLocus& locus, GffObj* lnkmrna) {
    //add all the elements of the other locus (merging)
    //-- merge mexons
   GArray<int> ovlexons(true,true); //list of locus.mexons indexes overlapping existing mexons
   int i=0; //index of first mexons with a merge
   int j=0; //index current mrna exon
   while (i<mexons.Count() && j<locus.mexons.Count()) {
     uint istart=mexons[i].start;
     uint iend=mexons[i].end;
     uint jstart=locus.mexons[j].start;
     uint jend=locus.mexons[j].end;
     if (iend<jstart) { i++; continue; }
     if (jend<istart) { j++; continue; }
     //if (mexons[i].overlap(jstart, jend)) {
     //exon overlap was found :
       ovlexons.Add(j);
       //extend mexons[i] as needed
       if (jstart<istart) mexons[i].start=jstart;
       if (jend>iend) { //mexons[i] end extend
            mexons[i].end=jend;
            //now this could overlap the next mexon(s), so we have to merge them all
            while (i<mexons.Count()-1 && mexons[i].end>mexons[i+1].start) {
                uint nextend=mexons[i+1].end;
                mexons.Delete(i+1);
                if (nextend>mexons[i].end) {
                     mexons[i].end=nextend;
                     break; //no need to check next mexons
                     }
                } //while next mexons merge
            } // mexons[i] end extend
       //  } //exon overlap
       j++; //check the next locus.mexon
      }
    // -- merge uexons
    //add to exons:
   GSeg seg;
   for (int i=0;i<locus.mexons.Count();i++) {
      seg.start=locus.mexons[i].start;
      seg.end=locus.mexons[i].end;
      if (!ovlexons.Exists(i)) mexons.Add(seg);
      }
    // -- add locus.mrnas
    for (int i=0;i<locus.mrnas.Count();i++) {
       locus.mrnas[i]->uptr=this;
       if (locus.mrnas[i]!=lnkmrna)
            mrnas.Add(locus.mrnas[i]);
       }
    // -- adjust start/end as needed
    if (start>locus.start) start=locus.start;
    if (end<locus.end) end=locus.end;
    }

  bool add_mRNA(GffObj* mrna) {
    if (mrnas.Count()>0 && mrna->gseq_id!=gseq_id) return false; //mrna must be on the same genomic seq
   //check for exon overlap with existing mexons
   //also update uexons and mexons accordingly, if mrna is added
    uint mrna_start=mrna->exons.First()->start;
    uint mrna_end=mrna->exons.Last()->end;
    if (mrna_start>end || start>mrna_end) return false;
    bool hasovl=false;
    int i=0; //index of first mexons with a merge
    int j=0; //index current mrna exon
    GArray<int> ovlexons(true,true); //list of mrna exon indexes overlapping mexons
    while (i<mexons.Count() && j<mrna->exons.Count()) {
      uint istart=mexons[i].start;
      uint iend=mexons[i].end;
      uint jstart=mrna->exons[j]->start;
      uint jend=mrna->exons[j]->end;
      if (iend<jstart) { i++; continue; }
      if (jend<istart) { j++; continue; }
      //exon overlap found if we're here:
      ovlexons.Add(j);
      hasovl=true;
      //extend mexons[i] as needed
      if (jstart<istart) mexons[i].start=jstart;
      if (jend>iend) { //mexon stretch up
             mexons[i].end=jend;
             //now this could overlap the next mexon(s), so we have to merge them all
             while (i<mexons.Count()-1 && mexons[i].end>mexons[i+1].start) {
                 uint nextend=mexons[i+1].end;
                 mexons.Delete(i+1);
                 if (nextend>mexons[i].end) {
                      mexons[i].end=nextend;
                      break; //no need to check next mexons
                      }
                 } //while next mexons merge
             } //possible mexons merge

      j++; //check the next mrna exon
      }//all vs all exon check loop
    if (hasovl) {
     GSeg seg;
     for (int i=0;i<mrna->exons.Count();i++) {
       seg.start=mrna->exons[i]->start;
       seg.end=mrna->exons[i]->end;
       if (!ovlexons.Exists(i)) mexons.Add(seg);
       }
     // add to mrnas
     mrnas.Add(mrna);
     // adjust start/end
     if (start>mrna_start) start=mrna_start;
     if (end<mrna_end) end=mrna_end;
     mrna->uptr=this;
     gseq_id=mrna->gseq_id;
     }
    return hasovl;
   }
};


class GSeqData {
  public:
    GList<GffObj> mrnas;
    GList<GLocus> loci;
    GSeqData():mrnas(true,true,false), loci(true,true,true) {
      }
};

void openfw(FILE* &f, GArgs& args, char opt) {
  GStr s=args.getOpt(opt);
  if (!s.is_empty()) {
      if (s=='-')
       f=stdout;
      else {
       f=fopen(s,"w");
       if (f==NULL) GError("Error creating file: %s\n", s.chars());
       }
     }
}

// --- globals:
GHash<GSeqData> ctghash(true);

void cluster_mRNAs(GList<GffObj>& mrnas, GList<GLocus>& loci) {
//mrnas sorted by start coordinate
//and so are the loci
for (int t=0;t<mrnas.Count();t++) {
  GArray<int> mrgloci(true);
  GffObj* mrna=mrnas[t];
  int lfound=0; //count of parent loci
  for (int l=0;l<loci.Count();l++) {
    if (loci[l]->end<mrna->exons.First()->start) continue;
    if (loci[l]->add_mRNA(mrna)) {
       //a parent locus was found
       lfound++;
       mrgloci.Add(l);
       }
    }//loci loop
  if (lfound==0) {
     //create a locus with only this mRNA
     loci.Add(new GLocus(mrna));
     }
  else if (lfound>1) {
    //more than one loci found parenting this mRNA, merge loci
    //if (lfound>2) GMessage(" merging %d loci \n",lfound);
     for (int l=1;l<lfound;l++) {
          int mlidx=mrgloci[l]-l+1;
          loci[mrgloci[0]]->addMerge(*loci[mlidx], mrna);
          loci.Delete(mlidx);
          }
     }
  }//mrnas loop
}

class COvlCovList: public GList<CovIsland> {
  public:
    GLocus* ovloc;
    COvlCovList():GList<CovIsland>(true,true,false) {
      ovloc=NULL;
      }
    void addCov(CovIsland& cov, GLocus* loc=NULL) {
      //add CovIsland overlapping locus
      if (this->Count()==0 && loc==NULL)
         GError("Error: cannot add CovIsland to COvlCovList without a locus");
      if (loc!=NULL) {
          if (!loc->overlap(cov.start, cov.end))
            GError("Error: adding non-overlapping CovIsland %d-%d to COvlCovList for locus %d-%d!\n",
                    cov.start, cov.end, loc->start, loc->end);
          if (ovloc!=NULL && ovloc!=loc)
            GError("Error: changing locus when adding new overlapping CovIsland\n");
          ovloc=loc;
          }
      this->Add(new CovIsland(cov));
      }
   bool mrnaOvl(GffObj& mrna, bool fullcov) {
     //TODO: check single or all-exon match
     int i=0;//cov index
     int j=0;//exon index
     int* exonovl=NULL; //exonovl[j]=1 if exon j has overlap
     GCALLOC(exonovl, mrna.exons.Count()*sizeof(int));
     bool aovl=false; //at least an ovl
     while (j<mrna.exons.Count() && i<Count()) {
        uint jstart=mrna.exons[j]->start;
        uint jend=mrna.exons[j]->end;
        uint istart=Get(i)->start;
        uint iend=Get(i)->end;
        if (iend<jstart) { i++; continue; } //next covisland
        if (jend<istart) { j++; continue; } //next exon
        //we have overlap!
        exonovl[j]=1;aovl=true;
        if (!fullcov) { GFREE(exonovl); return true; }
        //     if (iend<jend) i++;
        //               else j++;
        j++; //check next exon
        }
     if (!aovl) { GFREE(exonovl); return false; }
     if (fullcov)
        for (j=0;j<mrna.exons.Count();j++) {
            if (exonovl[j]!=1) { GFREE(exonovl); return false; }
            }
     GFREE(exonovl);
     return true;
    }

   void ovlCheck(bool fullcov, FILE* fout) {
      if (Count()==0) return;
      for (int m=0;m<ovloc->mrnas.Count();m++) {
         if (ovloc->mrnas[m]->udata==0 && mrnaOvl(*(ovloc->mrnas[m]), fullcov)) {
             ovloc->mrnas[m]->udata=1;
             ovloc->mrnas[m]->printGtf(fout);
             }
         }//for each mrna in this locus
     // now clear it, we're done with this
      ovloc=NULL;
      Clear();
     }

};


void load_mRNAs(GStr& fname) {
 // reads mRNAs into ctghash and clusters them into loci
 GMessage("Reading mRNAs..");
 FILE* f=fopen(fname.chars(),"r");
 if (f==NULL) GError("Error opening gff file %s\n",fname.chars());
 GffReader* gffr=new GffReader(f,true);
 gffr->readAll();
 for (int k=0;k<gffr->gflst.Count();k++) {
    GffObj* m=gffr->gflst[k];
    const char* gseqname=m->getGSeqName();
    GSeqData* gdata=ctghash.Find(gseqname);
    if (gdata==NULL) {
          gdata=new GSeqData();
          ctghash.shkAdd(gseqname,gdata);
          }
    gdata->mrnas.Add(m);
    }//for each mRNA loaded
  GMessage("  .. %d mRNAs loaded.\n",gffr->gflst.Count());
  delete gffr;
  GMessage("Clustering mRNAs into loci..\n");
  // -- now for each contig, cluster mRNAs into loci
  ctghash.startIterate();
  GSeqData* gdata=NULL;
  int numloci=0;
  while ((gdata=ctghash.NextData())!=NULL) {
       cluster_mRNAs(gdata->mrnas, gdata->loci);
       numloci+=gdata->loci.Count();
       }
  GMessage("  .. %d loci formed.\n",numloci);
  fclose(f);
 }

GList<GLocus>* getCtgLoci(GStr& ctgname) {
 GSeqData* r=ctghash.Find(ctgname.chars());
 return (r==NULL ? NULL : &(r->loci));
 }

#define FRCLOSE(fh) if (fh!=NULL && fh!=stdin) fclose(fh)
#define FWCLOSE(fh) if (fh!=NULL && fh!=stdout) fclose(fh)

int minintron = 12; //used for crossing small "exon" gaps

int main(int argc, char * const argv[]) {
 GArgs args(argc, argv, "Ahf:o:");
 int e;
 if ((e=args.isError())>0)
    GError("%s\nInvalid argument: %s\n", USAGE, argv[e]);
 if (args.getOpt('h')!=NULL) GError("%s\n", USAGE);
 bool mrna_filter=false;
 GStr s=args.getOpt('f');
 if (!s.is_empty()) {
   load_mRNAs(s); //also clusters them into loci
   mrna_filter=true;
   }
 FILE* f_in=NULL;
 FILE* f_out=NULL;
 openfw(f_out, args, 'o');
 if (f_out==NULL) f_out=stdout;
 if (args.startNonOpt()==0) GError("%s\n", USAGE);
 char* infile=args.nextNonOpt();
 f_in=fopen(infile,"r");
 if (f_in==NULL) GError("Cannot open input file %s!\n",infile);
 bool fullOvl=(mrna_filter && args.getOpt('A')==NULL);
 GLineReader covreader(f_in);
 GStr gseqname;

 if (!mrna_filter) { //simply write the coverage islands and quit
   while (!covreader.eof()) {
     CovIsland cov(covreader, gseqname, minintron);
     if (cov.start==0) continue; //no valid coverage parsed, try next line
     fprintf(f_out, "%s\t%d\t%d\n", gseqname.chars(),cov.start, cov.end);
    }//while reading cov islands
   FWCLOSE(f_out);
   FRCLOSE(f_in);
   return 0;
  }

 GStr prevgseq;
 /* mrna_filter method:
  * when switching to a new contig
  * *skip CovIslands until we find one overlapping the next locus on that contig
  * *create a new GList<CovIsland> and add CovIslands while they are in range of current locus
  * * when the next CovIsland falls outside the current locus, do the overlapCheck
  *   between the GList<CovIsland> and every mRNA in the locus
  * * only output those mRNAs having overlaps for EACH exon with any of GList<CovIsland>
  */
 COvlCovList covlst; //CovIslands overlapping current locus
 GList<GLocus>* ctgloci=NULL;
 int locidx=0; //current locus being checked
 while (!covreader.eof()) {
   CovIsland cov(covreader, gseqname, minintron);
   if (cov.start==0) continue;  //no valid coverage parsed, try next line
   if (prevgseq!=gseqname) { //contig change?
      prevgseq=gseqname;
      covlst.ovlCheck(fullOvl, f_out);
      ctgloci=getCtgLoci(gseqname);
      locidx=0;
      }// contig change

   //entry point when advancing locidx, to recheck the current cov
  checkNextLocus:
   if (ctgloci==NULL) {
     continue; //no loci on this contig, try next
     }
   if (locidx>=ctgloci->Count()) {
     //check out what we got so far
     //covlst.ovlCheck(fullOvl, f_out);
     ctgloci=NULL;
     continue;
     }
   //ctgloci !=NULL, there are mRNAs for this contig
   //if (covlst.Count()==0) { //no overlap found yet
   GLocus* curloc=ctgloci->Get(locidx);
   if (curloc->end<cov.start) {
        //locus is behind current cov, has to catch up
        //check what we have so far and move on
        covlst.ovlCheck(fullOvl, f_out);
        locidx++;
        goto checkNextLocus; //same CovIsland
        } //let
    // -- here we have current locus->end >= cov.start
    if (cov.end<curloc->start) {
       //cov behind locus
       //check what we have so far and move on
       covlst.ovlCheck(fullOvl, f_out);
       //try next covisland against the same locus
       continue; //same locus, next covisland if any
       }
    //--here we have overlap of cov with current locus
    covlst.addCov(cov,curloc);
    covlst.ovlCheck(fullOvl, f_out);
    locidx++;
    goto checkNextLocus; //see if there are more loci overlapping this covisland
    }//while reading wig lines
 covlst.ovlCheck(fullOvl,f_out);    
 FWCLOSE(f_out);
 FRCLOSE(f_in);
 //GMessage(".. press Enter to exit ..\n");
 //getchar();
}
