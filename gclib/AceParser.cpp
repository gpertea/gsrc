#include "AceParser.h"
#include <ctype.h>

bool AceParser::open() { //also checks if it looks like a valid ACE file
   if (f==stdin) return true;
   if ((f=fopen(fname,"rb"))==NULL) return false;
   char first3[4];
   first3[3]='\0';
   if (fread(first3, 1,3, f)==3) {
     f_pos=3;
     return (strcmp(first3, "AS ")==0 || strcmp(first3, "CO ")==0);
     }
   return false;
   }


LytSeqInfo* AceParser::addSeq(char* s, LytCtgData* ctg) {
 LytSeqInfo* seq;
 if (!startsWith(s, "AF ", 3))
       return NULL;
 s+=3;
 char* p=strchrs(s," \t");
 if (p==NULL) return NULL;
 p++;
 char c;
 int offs;
 if (sscanf(p,"%c %d", &c, &offs)!=2) return NULL;
 p--;
 *p='\0';
 if ((seq=seqinfo.Find(s))!=NULL) {
   GMessage("Sequence '%s' already found for contig '%s (%d nt)'\n"
      " so it cannot be added for contig '%s (%d)'\n",
     s, seq->contig->name, seq->contig->len,
     ctg->name, ctg->len);
   return NULL;
   }
 seq = new LytSeqInfo(s, ctg, offs, (c=='C') ? 1 : 0);
 seqinfo.shkAdd(seq->name, seq);
 ctg->seqs.Add(seq);
 return seq;
}


bool AceParser::loadContig(int ctgidx, fnLytSeq* seqfn, bool re_pos) {

    bool forgetCtg = false;
    if (ctgidx>=contigs.Count())
      GError("LayoutParser: invalid contig index '%d'\n", ctgidx);
    LytCtgData* ctgdata=contigs[ctgidx];
    if (re_pos && currentContig!=NULL) { //free previously loaded contig data
      currentContig->seqs.Clear();       // unless it was a parse() call
      seqinfo.Clear();
      }
    currentContig=ctgdata;
    int ctg_numSeqs=ctgdata->numseqs;

    if (re_pos) {
       seek(ctgdata->fpos); //position right where the contig definition starts
       char *r = linebuf->getLine(f,f_pos);
       if (r==NULL) return false;
       }

    if (seqfn!=NULL) { //process the contig sequence!
       char* ctgseq=readSeq();
       forgetCtg=(*seqfn)(numContigs, ctgdata, NULL, ctgseq);
       GFREE(ctgseq); //obviously the caller should have made a copy
       }
    //now look for all the component sequences
    if (fskipTo("AF ")<0) {
       GMessage("AceParser: error finding sequence offsets (AF)"
                " for contig '%s' (%d)\n", ctgdata->name, ctgdata->len);
       return false;
       }
    int numseqs=0;
    while (startsWith(linebuf->chars(), "AF ",3)) {
      if (addSeq(linebuf->chars(), ctgdata)==NULL) {
         GMessage("AceParser: error parsing AF entry:\n%s\n",linebuf->chars());
         return false;
         }
       numseqs++;
      //read next line:
      linebuf->getLine(f,f_pos);
      }
    if (numseqs!=ctg_numSeqs) {
      GMessage("Invalid number of AF entries found (%d) for contig '%s' "
         "(length %d, numseqs %d)\n", numseqs,
                ctgdata->name, ctgdata->len, ctg_numSeqs);
      return false;
      }
    //now read each sequence entry
    off_t seqpos=fskipTo("RD ");
    numseqs=0; //count again, now the RD entries
    if (seqpos<0) {
         GMessage("AceParser: error locating first RD entry for contig '%s'\n",
             ctgdata->name);
         return false;
         }
    //int numseqs=0;
    //reading the actual component sequence details
    while (startsWith(linebuf->chars(), "RD ",3)) {
      char* s=linebuf->chars()+3;
      char* p=strchrs(s, " \t");
      LytSeqInfo* seq;
      if (p==NULL) {
          GMessage("AceParser: Error parsing RD header line:\n%s\n", linebuf->chars());
          return false;
          }
      *p='\0';
      if ((seq=seqinfo.Find(s))==NULL) {
          GMessage("AceParser: unknown RD encountered: '%s'\n", s);
          return false;
          }
      p++; //now p is in linebuf after the RD name
      seq->fpos=seqpos;
      int len;
      if (sscanf(p, "%d", &len)!=1) {
          GMessage("AceParser: cannot parse RD length for '%s'\n", s);
          return false;
          }
      seq->setLength(len);
      //read the sequence data here if a callback fn was given:
      char* sseq=NULL;
      if (seqfn!=NULL)
          sseq=readSeq(seq); //read full sequence here
      if (fskipTo("QA ")<0) {
           GMessage("AceParser: Error finding QA entry for read %s! (fpos=%llu)\n", seq->name, (unsigned long long)f_pos);
           return false;
           }
      //parse QA entry:
      int tmpa, tmpb;
      if (sscanf(linebuf->chars()+3, "%d %d %d %d", &tmpa, &tmpb, &seq->left,&seq->right)!=4 ||
             seq->left<=0 || seq->right<=0) {
           GMessage("AceParser: Error parsing QA entry.\n");
           return false;
           }
      /*
      if (fskipTo("DS")<0) {
           GMessage("AceParser: Error closing RD entry ('DS' not found).\n");
           return false;
           }
           */
      seqpos=getFilePos()+1;
      bool forgetSeq=false;
      if (seqfn!=NULL) {
          forgetSeq=(*seqfn)(numContigs, ctgdata, seq, sseq);
          GFREE(sseq);
          }
      if (forgetSeq) { //parsing the whole stream -- aceconv)
            ctg_numSeqs--;
            seqinfo.Remove(seq->name);
            ctgdata->seqs.RemovePtr(seq);
            }
      numseqs++;
      if (numseqs<ctgdata->numseqs)
        seqpos=fskipTo("RD ", "CO "); //more sequences left to read
      }
    if (numseqs!=ctgdata->numseqs) {
      GMessage("Error: Invalid number of RD entries found (%d) for contig '%s' "
         "(length %d, numseqs %d)\n", numseqs,
                ctgdata->name, ctgdata->len, ctg_numSeqs);
      return false;
      }
    if (forgetCtg) {
      ctgIDs.Remove(ctgdata->name);
      ctgdata->seqs.Clear();
      seqinfo.Clear();
      contigs.RemovePtr(ctgdata);
      }
 return true;
}

//read only the contig headers from the file
//also checks for duplicate seqnames (just in case)
bool AceParser::parseContigs() {
  if (f!=stdin)  seek(0);
  off_t ctgpos;
  numContigs=0;
  while ((ctgpos=fskipTo("CO "))>=0) {
    numContigs++;
    LytCtgData* ctgdata=new LytCtgData(ctgpos);
    char* p=ctgdata->readName(linebuf->chars()+3, ctgIDs);
    if (p==NULL) {
       GMessage("AceParser: error parsing contig name!\n");
       return false;
       }
    int ctglen, ctg_numSeqs, numseqs;
    //p must be after contig name within linebuf!
    if (sscanf(p, "%d %d", &ctglen, &numseqs)!=2) {
      GMessage("Error parsing contig len and seq count at:\n%s\n",
         linebuf->chars());
      }
    ctg_numSeqs=numseqs;
    ctgdata->len=ctglen;
    ctgdata->numseqs=numseqs;
    ctgdata->offs = 1;
    contigs.Add(ctgdata);
    //ctgpos=fskipTo("CO ");
    }
  if (ctgpos==-2) return false; //parsing failed: line too long ?!?
  contigs.setSorted(true);
  return true;
}


bool AceParser::parse(fnLytSeq* seqfn) {
  //read all seqs and their positions from the file
  //also checks for duplicate seqnames (just in case)
  if (f!=stdin)  seek(0);
  ctgIDs.Clear();
  //GHash<int> ctgIDs; //contig IDs, to make them unique!
  //
  off_t ctgpos;
  numContigs=0;
  while ((ctgpos=fskipTo("CO "))>=0) {
    numContigs++;
    LytCtgData* ctgdata=new LytCtgData(ctgpos);
    char* p=ctgdata->readName(linebuf->chars()+3, ctgIDs);
    if (p==NULL) {
       GMessage("AceParser: error parsing contig name!\n");
       return false;
       }
	int ctglen, numseqs;
	//p must be after contig name within linebuf!
	if (sscanf(p, "%d %d", &ctglen, &numseqs)!=2) {
	  GMessage("Error parsing contig len and seq count at:\n%s\n",
		 linebuf->chars());
	  }
	ctgdata->len=ctglen;
	ctgdata->numseqs = numseqs;
	ctgdata->offs = 1;
	int ctgidx=contigs.Add(ctgdata);
	loadContig(ctgidx, seqfn, false);
	} //while contigs
  if (ctgpos==-2) return false; //parsing failed: line too long ?!?
  contigs.setSorted(true);
  return true;
  }


 void add_intron(char*& buf,int accrd, int& rgpos,
        const char* numstr, LytSeqInfo* seqinfo, char lspl, char rspl) {
	int n=atoi(numstr);
	if (n<=0) return;
	if (seqinfo!=NULL) {
	  char sL = 0;
      char sR = 0;
      if (lspl>0) sL=(lspl=='['?'S':'s');
      if (rspl>0) sR=(rspl==']'?'S':'s');
	  seqinfo->addInterSeg(seqinfo->offs+rgpos-1,seqinfo->offs+rgpos+n,
						   0,0,sL,sR,accrd);
	  }
	rgpos+=n;
}

char* AceParser::readSeq(LytSeqInfo* sqinfo) {
//assumes the next line is where a sequence starts!
//stops at the next empty line encountered
char* buf;
static char rlenbuf[12]={0,0,0,0,0,0,0,0,0,0,0,0}; //buffer for parsing the gap length
int rlenbufacc=0; //how many digits accumulated in rlenbuf so far
int buflen=512;
GMALLOC(buf, buflen); //this MUST be freed by the caller
buf[0]='\0';
int accrd=0; //accumulated read length so far -- excludes interseg gaps!
int rgpos=0; //accumulated offset including interseg gaps!
char *r = linebuf->getLine(f,f_pos);
int linelen=linebuf->length();
char r_splice=0, l_splice=0;
while (linelen>0) {
   if (r==NULL) {
	  GMessage("AceParser: error reading sequence data\n");
	  return NULL;
	  }
   //-- pass the line content for accumulation
   int i=0;
   while (r[i]) {
	 while (r[i] && isdigit(r[i])) {
		 rlenbuf[rlenbufacc]=r[i];
		 rlenbufacc++;
		 i++;
		 }
	 if (r[i]==0) break; //end of line reached already
	 //now r[i] is surely a non-digit char
	 if (rlenbufacc>0) { //have we just had a number before?
		rlenbuf[rlenbufacc]=0;
		if (r[i]=='=' || r[i]=='-') {
		     i++;
		     if (r[i]==0) break;
		   } else { //check for splice site markers for this introns
		     if (r[i]=='('||r[i]=='[') {
		         l_splice=r[i];
		         i++;
		         if (r[i]==0) break;
		         }
		     if (r[i]==')'||r[i]==']') {
                 r_splice=r[i];
                 i++;
                 if (r[i]==0) break;
		         }
		   }//splice check
		add_intron(buf, accrd, rgpos, rlenbuf, sqinfo, l_splice, r_splice);
		rlenbufacc=0;
		r_splice=0;l_splice=0;
		//i++;//skip the gap character
		}
	 //check for digits here and break the linebuf as needed
	 int bi=i; //start of non-digit run
	 while (r[i] && !isdigit(r[i])) i++;
	 int nl=(i-bi); //length of non-digit run
	 if (nl>0) {
		int si=accrd;
		accrd+=nl;
		rgpos+=nl;
		if (accrd>=buflen-1) {
			buflen=accrd+512;
			GREALLOC(buf,buflen);
			}
		//append these non-digit chars
		for(int b=0;b<nl;b++) {
			buf[si+b]=r[bi+b];
			}
		}//non-digit run
	 } //while line chars

   /*
   //-- append the line to buf
   accrd+=linelen;
   if (accrd>=buflen) {
		buflen+=1024;
		GREALLOC(buf,buflen);
		}
   strcat(buf, r);
   */

   r=linebuf->getLine(f,f_pos);
   linelen=linebuf->length();
   }//while linelen>0

//add the 0-ending
buf[accrd]=0;
return buf;
}

char* AceParser::getSeq(LytSeqInfo* seq) {
 if (f==stdin || seek(seq->fpos)!=0) {
   GMessage("AceParser: error seeking seq '%s'\n", seq->name);
   return NULL;
   }
 //skip the contig header:
 char* r=linebuf->getLine(f,f_pos);
 if (r==NULL || !startsWith(linebuf->chars(), "RD ", 3)) {
	GMessage("AceParser: error seeking seq '%s'\n"
		" (no RD entry found at location %d)\n", seq->name, seq->fpos);
    return NULL;
	}
 return readSeq(seq);
}

char* AceParser::getContigSeq(LytCtgData* ctg) {
 if (f==stdin || seek(ctg->fpos)!=0) {
   GMessage("AceParser: error seeking contig '%s'\n", ctg->name);
   return NULL;
   }
 //skip the contig header:
 char* r=linebuf->getLine(f,f_pos);
 if (r==NULL || !startsWith(linebuf->chars(), "CO ", 3)) {
    GMessage("AceParser: error seeking contig '%s'\n"
        " (no CO entry found at location %d)\n", ctg->name, ctg->fpos);
    return NULL;
    }
 return readSeq();
}


