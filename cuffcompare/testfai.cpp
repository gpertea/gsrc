/*
 * testfai.cpp
 *
 *  Created on: Sep 7, 2010
 *      Author: gpertea
 */

#include "GArgs.h"
#include "GFaSeqGet.h"
#include "GFastaIndex.h"
#include "GStr.h"
#include <ctype.h>

#define USAGE "Usage:\n\
  testfai <multi-fasta> [-F] [-q <seqname>] [-r <range>] [-C] [-o <out.fasta>]\n\
  \n\
  Create the index <multi-fasta>.fai if not available and retrieve the\n\
  fasta sequence for <seqname>, optionally only pulling range <range>\n\
  (reverse complemented if -C is given)\n\
  \n\
  Options:\n\
   -F : force rebuilding of the fasta index in the current directory\n\
  "
void openfwrite(FILE* &f, GArgs& args, char opt) {
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

int main(int argc, char * const argv[]) {
  GArgs args(argc, argv, "hFCq:r:o:");
  int e;
  if ((e=args.isError())>0)
    GError("%s\nInvalid argument: %s\n", USAGE, argv[e]);
    if (args.getOpt('h')!=NULL){
      GMessage("%s\n", USAGE);
                exit(1);
      }
  args.startNonOpt();
  GStr fadb(args.nextNonOpt());
  if (fadb.is_empty()) GError("%s Error: multi-fasta file expected!\n",USAGE);
  GStr fname(fadb);
  fname.append(".fai");
  bool createLocal=(args.getOpt('F')!=NULL);
  const char* idxname=(createLocal)? NULL : fname.chars();
  GFastaIndex faidx(fadb.chars(), idxname);
  //also tried to load the index if exists in the current directory
  GStr fnamecwd(fname); //name in current directory (without path)
  int ip=-1;
  if ((ip=fnamecwd.rindex(CHPATHSEP))>=0) {
    fnamecwd.cut(0,ip+1);
    }
  if (!createLocal) { //look for existing indexes to load
    //try the same directory as the fasta file first
    if (!faidx.hasIndex() and fileExists(fnamecwd.chars())>1) { //try current working directory next
        faidx.loadIndex(fnamecwd.chars());
       }
    if (!faidx.hasIndex()) {//could not load any index data
       //try to create it in the same directory as the fasta file
       GMessage("No fasta index found. Rebuilding..\n");
       faidx.buildIndex();
       if (faidx.getCount()==0) GError("Error: no fasta records to be indexed!\n");
       GMessage("Fasta index rebuilt.\n");
       //check if we can create a file there
       FILE* fcreate=fopen(fname.chars(), "w");
       if (fcreate==NULL)
         GMessage("Warning: cannot create fasta index %s! (permissions?)\n", fname.chars());
       else {
         fclose(fcreate);
         if (faidx.storeIndex(fname.chars())<faidx.getCount())
           GMessage("Warning: error writing the index file %s!\n",fname.chars());
         } //creating index file in the same directory as fasta file
       }//trying to create the index file
    }

  if (createLocal || !faidx.hasIndex()) {
    //simply rebuild the index in the current directory and use it:
    //remove directories in path, if any
    if (faidx.getCount()==0) {
          faidx.buildIndex();
          if (faidx.getCount()==0) GError("Error: no fasta records to be indexed!\n");
          }
    if (faidx.storeIndex(fnamecwd.chars())<faidx.getCount())
        GMessage("Warning: error writing the index file %s!\n",fnamecwd.chars());
    }

  GStr qry(args.getOpt('q'));
  if (qry.is_empty()) exit(0);
  GFastaRec* farec=faidx.getRecord(qry.chars());
  if (farec==NULL) {
      GMessage("Error: couldn't find fasta record for '%s'!\n",qry.chars());
      exit(1);
      }
  GFaSeqGet faseq(fadb.chars(),farec->seqlen, farec->fpos, farec->line_len, farec->line_blen);
  //TODO: read these from -r option
  uint qstart=0;
  uint qend=0; //farec->seqlen
  bool revCompl=(args.getOpt('C')!=NULL);
  char* s=args.getOpt('r');
  if (s!=NULL) {
     char *p=s;
     while (isdigit(*p)) p++;
     if (*p=='-') {
                sscanf(s,"%u-%u",&qstart, &qend);
                if (qstart==0 || qend==0)
                      GError("Error parsing sequence range: %s\n",s);
                }
       else if (*p==':') {
                int qlen=0;
                sscanf(s,"%u:%d", &qstart, &qlen);
                if (qstart==0 || qlen==0)
                     GError("Error parsing sequence range: %s\n",s);
                qend=qstart+qlen-1;
                }
       else if (*p=='.') {
               sscanf(s,"%u..%u",&qstart, &qend);
               if (qstart==0 || qend==0)
               GError("Error parsing sequence range: %s\n",s);
               }
     }
  if (qstart==0) qstart=1;
  if (qend==0) qend=farec->seqlen;
  // call faseq.loadall() here if multiple ranges are to be extracted all
  // over this genomic sequence
  char* subseq=faseq.copyRange(qstart, qend, revCompl, true);
  FILE* f_out=NULL;
  openfwrite(f_out, args, 'o');
  if (f_out==NULL) f_out=stdout;
  writeFasta(f_out, qry.chars(), NULL, subseq, 70, qend-qstart+1);
  GFREE(subseq);
}

