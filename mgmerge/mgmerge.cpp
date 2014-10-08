//un/comment this for tracing options
//#define DBGTRACE 1
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "GArgs.h"
#include "GStr.h"
#include "GList.hh"

#define usage "Usage: mgmerge [-C 'cmd'][-M][-o <fbasename> [-s <zsize>] [-V]] files...\n\
 \n\
 Multi-way merge of sorted [and compressed] mgblast output files. Unless -M option\n\
 is given, the input files must be already sorted in this order: \n\
    increasingly by column 11 (E-value)\n\
    then decreasingly by 10ths and then 9th columns (score, pid)\n\
 Optionally, the sorted output stream can be redirected\n\
 to a series of compressed files (if -o option is provided).\n\
 \n\
-C'<filtercmd>' also pipe sorted lines into command '<filtercmd>' and use\n\
                its output instead for -o option\n\
-o<fbasename>   redirect output to gzip to compress it to series of \n\
                sequential files named <fbasename>_NNNN.gz\n\
-M              simply merge the input, with no comparison/sorting\n\
-b              for -o option, use bzip2 with maximum compression instead\n\
                of gzip; output files will have a .bz2 extension instead of .gz\n\
-s<zsize>       only for -o option: specify the approximate maximum size, \n\
                in MegaBytes, of each of the <fbasename>_NNNN.gz files \n\
                (default: 500)\n\
-V              verbose\n"


#define BUF_LEN 8192 //maximum input line length 
#define ERR_INVALID_PAIR "Invalid input line encountered:\n%s\n"

//============ global vars and structures:
/*bool is_numeric=false;
bool is_descending=false;*/
bool verbose=false;
bool use_bzip2=false;
bool nosort=false; // -M option
//int field_no=-1; /* whole line */
//char delimiter[2]="\t";
GStr fbasename;
GStr midpipecmd; //-C command filter
int zmaxsize=500;
int partno=0;
GStr outcmd;

FILE* fout=NULL; /* this could be a popen handle */

class CLine {
public:
 int fileno;
 GStr line;
 int score;
 float pid;
 double evalue;
 bool operator==(CLine& l){
     return (evalue==l.evalue && score==l.score && pid==l.pid);
     }
   bool operator<(CLine& l){ // this is the "better" line
      if (evalue==l.evalue) { 
         if (score==l.score) return (pid<l.pid);
                        else return (score<l.score);
         }
        else return (evalue>l.evalue);
     }
     
   CLine(int fno, const char* s) { //constructor
     GStr field;
     char c=0;
     line=s;
     fileno=fno;
     pid=0;
     score=0;
     evalue=-1;
     line.startTokenize("\t", tkFullString);
     while (c<=10 && line.nextToken(field)) {
       c++;
       switch (c) {
         case 9: pid=field.asReal();break;
         case 10: score=field.asInt();break;
         case 11: evalue=field.asReal();break;
         } //switch
       } //while tokenize
     if (score==0 || pid==0 || evalue<0) {
       GError("Error at parsing input line:\n%s\n",s);
       }
     }
};

GList<CLine> lines(true,true,false);

FILE** srcfiles;
char** fnames;

char* readLine(char* buf, int bufsize, int fno) {
 char* line=fgets(buf, bufsize-1, srcfiles[fno]); 
 if (line==NULL) return NULL;
 int l=strlen(line);
 if (line[l-1]!='\n')
   GError("Input line too long from file '%s':\n%s\n", 
      line, fnames[fno]);
 line[l-1]='\0';     
 return line;
 }

void writeLine(const char* line, int& lno) {
 struct stat filestat;
 GStr partname;
 fprintf(fout, "%s\n", line);
 lno++;
 if (partno>0) { /* we write into gzip's pipe */
   if ( lno%500 != 0 ) return; //only check every 500 lines
   if (use_bzip2)
       partname.format("%s_%03d.bz2", fbasename.chars(), partno);
     else
       partname.format("%s_%03d.gz", fbasename.chars(), partno);
   //fflush(stdout);sync();
   if ( stat(partname.chars(), &filestat)==0 ) {
      if (filestat.st_size>zmaxsize) {
        if (pclose(fout)==-1)
             GError("Error closing compressed writing pipe to '%s'\n",
                 partname.chars());         
        partno++;
        lno=0;
        if (use_bzip2) {
             partname.format("%s_%03d.bz2", fbasename.chars(), partno);
             outcmd.format("%sbzip2 -c9 > %s", midpipecmd.chars(),
                                                           partname.chars());
             }
         else {
             partname.format("%s_%03d.gz", fbasename.chars(), partno);
             outcmd.format("%sgzip -f > %s", midpipecmd.chars(),
                                                    partname.chars());
             }
         if ((fout=popen(outcmd.chars(), "w"))==NULL) {
           perror("popen()");
           GError("Error popen() to '%s'\n", outcmd.chars());
           }
         }
        }
     else if (lno>10000)
        GError("Error getting file size for '%s'\n",partname.chars());
    
   }
}

//========================================================
//====================     main      =====================
//========================================================
int main(int argc, char * const argv[]) {
 char inbuf[BUF_LEN]; // incoming buffer for sequence lines.
 //GArgs args(argc, argv, "hVbrno:s:f:d:");
 GArgs args(argc, argv, "hVMbs:C:o:d:");
 int e,i;
 if ((e=args.isError())>0)
    GError("%s\nInvalid argument: %s\n", usage, argv[e]);
 if (args.getOpt('h')!=NULL) GError("%s\n", usage);


 use_bzip2=(args.getOpt('b')!=NULL);
 //is_numeric=(args.getOpt('n')!=NULL);
 //is_descending=(args.getOpt('r')!=NULL);
 verbose=(args.getOpt('V')!=NULL);
 nosort=(args.getOpt('M')!=NULL);
 //field_no=-1; /* whole line */
 //GStr s=args.getOpt('f');
 //if (!s.is_empty()) {
 //   field_no=s.asInt()-1;
 //   if (field_no<0) GError("Invalid field number provided.\n");
 //   }
 midpipecmd=args.getOpt('C');
 if (!midpipecmd.is_empty()) midpipecmd.append("| ");
 fbasename=args.getOpt('o');
 GStr s=args.getOpt('s');
 if (!s.is_empty()) {
   if (fbasename.is_empty())
     GError("%s\n Option -s cannot be given without the -o option\n", usage);
   zmaxsize=s.asInt();
   if (zmaxsize<=0 || zmaxsize>2000)
    GError("%s\nInvalid compressed file size: %d MB (must be a pozitive number<2000)\n", 
            usage, zmaxsize);
   }
 zmaxsize = zmaxsize*1024*1024;
 int numfiles=args.startNonOpt();   
 if (numfiles==0)
    GError("No input files to process.\n");
 GMALLOC(srcfiles, sizeof(FILE*) * numfiles);
 GMALLOC(fnames, sizeof(char*) * numfiles);
 GStr cmd;
 //open all files
 for (i=0;i<numfiles;i++) {
  fnames[i]=(char*)args.nextNonOpt();
  GStr fname=fnames[i];
  //GMessage("Opening file '%s'\n", fname.chars());
  int dotidx=fname.rindex('.');
  if (dotidx>=0) dotidx-=fname.length();
  GStr ending=fname.substr(dotidx);
  ending.lower();
  if (ending==".z" || ending==".gz") {
    cmd.format("gzip -cd '%s'", fname.chars());
    if ((srcfiles[i]=popen(cmd.chars(), "r"))==NULL) {
        perror("popen()");
        GError("Error at popen('%s')\n", cmd.chars());
        }
    }
   else 
     if (ending==".bz" || ending==".bz2") {
      cmd.format("bzip2 -cd '%s'", fname.chars());
      if ((srcfiles[i]=popen(cmd.chars(), "r"))==NULL) {
         perror("popen()");
         GError("Error at popen('%s')\n", cmd.chars());
         }
      }
     else { //uncompressed files assumed
       if ((srcfiles[i]=fopen(fname.chars(), "r"))==NULL)
        GError("Error at fopen %s \n", fname.chars());
       }
  }


int wlno=0;
int last;
int lno=0;

// prepare/open output file/pipe
if (!fbasename.is_empty()) {
  partno++; //the first part is always 001
  if (use_bzip2)
      outcmd.format("%sbzip2 -c9 > %s_%03d.bz2", midpipecmd.chars(), fbasename.chars(), partno);
    else 
      outcmd.format("%sgzip -f > %s_%03d.gz", midpipecmd.chars(), fbasename.chars(), partno);
  if ((fout=popen(outcmd.chars(), "w"))==NULL) {
       perror("popen()");
       GError("Error popen() to '%s'\n", outcmd.chars());  
       }
  }
 else fout=stdout;


//special case: just stream all input to the output
if (nosort) {
 for (i=0;i<numfiles; i++) {
   char* line;
   while ((line=readLine(inbuf, BUF_LEN-1, i))!=NULL) {
    wlno++;
    writeLine(line, lno);
    } //for each line
  } //for each input file
 goto END_MERGE;
 }

//- init: read the first lines from each file
for (i=0;i<numfiles; i++) {
 char* line=readLine(inbuf, BUF_LEN-1, i);
 //GMessage("%d: read: '%s'\n", i, line);
 if (line!=NULL)
   lines.Add(new CLine(i, line));
 }
/*
for (i=0;i<lines.Count();i++) {
 GMessage("%s\n", lines[i]->line.chars());
 }
GMessage("-----------------------\n"); 
*/

while ((last=lines.Count())>0) {
  last--;
  CLine* best=lines[last]; // best hit, print it
  wlno++;
  writeLine(best->line.chars(), lno);
  //read next one from this file
  int fno=best->fileno;
  char* line=readLine(inbuf, BUF_LEN-1, best->fileno);
  //GMessage("Removing line:'%s'\n", best->line.chars());
  lines.Delete(last); //never read from this file handle
  if (line!=NULL) {
    lines.Add(new CLine(fno, line));
    }
  }
END_MERGE:
if (verbose)
 GMessage("%d lines processed.\n", wlno);
if (fout!=stdout) {
  if (pclose(fout)==-1)
         GError("Error closing compressed writing pipe to the last part (fbasename: '%s')\n",
                 fbasename.chars());
  if (verbose)               
   GMessage("Written %d compressed parts (base name: '%s')\n",
        partno, fbasename.chars());  
  }

//closing time:
for (i=0;i<numfiles; i++) {
  GStr fname=fnames[i];
  GStr ending=fname.substr(-2);
  if (ending=="gz" || ending==".Z" || ending=="z2" || ending=="bz") 
     pclose(srcfiles[i]);
   else 
     fclose(srcfiles[i]);  
  }
 
 GFREE(srcfiles);
 GFREE(fnames);
} //-- end main()
