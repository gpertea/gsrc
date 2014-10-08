#include "GArgs.h"
#include "GStr.h"
#include "gtf_tracking.h"

GList<GSeqData> ref_data(true,true,true); //list of reference mRNAs and loci data for each genomic seq
              //each locus will keep track of any superloci which includes it, formed during the analysis

int main(int argc, char * const argv[]) {
  GArgs args(argc, argv, "FGCRPq:r:o:");
  GList<GSeqData>* rnas;
  args.startNonOpt();
  GStr gtfile(args.nextNonOpt());
  if (gtfile.is_empty()) GError("Usage: test_load [-R] [{-F|-G}] [{-P|-o <outgff>}] <gff_input>\nNo input GFF file given!\n");
  gtf_tracking_verbose=true;
  FILE* RNA_file=fopen(gtfile.chars(), "rb");
  if (RNA_file==NULL) GError("Error opening file %s!\n",gtfile.chars());
  //read_transcripts(mRNA_file, rnas, true);
  bool load_as_ref=(args.getOpt('R')!=NULL);
  int qfilter_redundant=1;
  if (args.getOpt('F')) qfilter_redundant=2;
  if (args.getOpt('G')) qfilter_redundant=0;
  if  (load_as_ref) {
       rnas=&ref_data;
       read_mRNAs(RNA_file, *rnas, NULL, 1, -1, gtfile.chars());
       }
      else {
        rnas=new GList<GSeqData>(true,true,true);
        read_mRNAs(RNA_file, *rnas, &ref_data, qfilter_redundant, 0, gtfile.chars());
        }
  fclose(RNA_file);
  int fcount=0;
  int rcount=0;
  int ucount=0;
  for (int i=0;i<rnas->Count();i++) {
     fcount+=rnas->Get(i)->mrnas_f.Count();
     rcount+=rnas->Get(i)->mrnas_r.Count();
     ucount+=rnas->Get(i)->umrnas.Count();
     }
  int mtotal=fcount+rcount+ucount;
  GStr oname(args.getOpt('o'));
  FILE* fout=NULL;
  if (!oname.is_empty()) fout=fopen(oname.chars(), "w");
      else if (args.getOpt('P')) fout=stdout;
  GMessage("Totals:\n%6d on forward strand\n"
     "%6d on reverse strand\n%6d unknown strand\n"
     "--------------------------\n%6d total\n",
     fcount, rcount, ucount, mtotal);
  if (fout) {
    for (int i=0;i<rnas->Count();i++) {
      for (int k=0;k<rnas->Get(i)->mrnas_f.Count();k++) {
        rnas->Get(i)->mrnas_f[k]->printGff(fout);
        }
      for (int k=0;k<rnas->Get(i)->mrnas_r.Count();k++) {
        rnas->Get(i)->mrnas_r[k]->printGff(fout);
        }
      for (int k=0;k<rnas->Get(i)->umrnas.Count();k++) {
        rnas->Get(i)->umrnas[k]->printGff(fout);
        }
      }
    if (fout!=stdout) fclose(fout);
    }
  if  (!load_as_ref) delete rnas;
  }
