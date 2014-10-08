#include "GBam.h"
#include <ctype.h>
#include "kstring.h"

//auxiliary functions for low level BAM record creation
uint8_t* realloc_bdata(bam1_t *b, int size) {
  if (b->m_data < size) {
        b->m_data = size;
        kroundup32(b->m_data);
        b->data = (uint8_t*)realloc(b->data, b->m_data);
        }
  if (b->data_len<size) b->data_len=size;
  return b->data;
}

uint8_t* dupalloc_bdata(bam1_t *b, int size) {
  //same as realloc_bdata, but does not free previous data
  //but returns it instead
  //it ALWAYS duplicates data
  b->m_data = size;
  kroundup32(b->m_data);
  uint8_t* odata=b->data;
  b->data = (uint8_t*)malloc(b->m_data);
  memcpy((void*)b->data, (void*)odata, b->data_len);
  b->data_len=size;
  return odata; //user must FREE this after
}

//from bam_import.c
extern unsigned short bam_char2flag_table[];

GBamRecord::GBamRecord(const char* qname, int32_t gseq_tid,
                 int pos, bool reverse, const char* qseq,
                 const char* cigar, const char* quals):exons(1) {
   novel=true;
   bam_header=NULL;
   b=bam_init1();
   if (pos<=0 || gseq_tid<0) {
               b->core.pos=-1; //unmapped
               b->core.flag |= BAM_FUNMAP;
               gseq_tid=-1;
               }
          else b->core.pos=pos-1; //BAM is 0-based
   b->core.tid=gseq_tid;
   b->core.qual=255;
   b->core.mtid=-1;
   b->core.mpos=-1;
   int l_qseq=strlen(qseq);
   //this may not be accurate, setting CIGAR is the correct way
   //b->core.bin = bam_reg2bin(b->core.pos, b->core.pos+l_qseq-1);
   b->core.l_qname=strlen(qname)+1; //includes the \0 at the end
   memcpy(realloc_bdata(b, b->core.l_qname), qname, b->core.l_qname);
   set_cigar(cigar); //this will also set core.bin
   add_sequence(qseq, l_qseq);
   add_quals(quals); //quals must be given as Phred33
   if (reverse) { b->core.flag |= BAM_FREVERSE ; }
   }

GBamRecord::GBamRecord(const char* qname, int32_t flags, int32_t g_tid,
             int pos, int map_qual, const char* cigar, int32_t mg_tid, int mate_pos,
             int insert_size, const char* qseq, const char* quals,
             GVec<char*>* aux_strings):exons(1)  {
  novel=true;
  bam_header=NULL;
  b=bam_init1();
  b->core.tid=g_tid;
  b->core.pos = (pos<=0) ? -1 : pos-1; //BAM is 0-based
  b->core.qual=map_qual;
  int l_qseq=strlen(qseq);
  b->core.l_qname=strlen(qname)+1; //includes the \0 at the end
  memcpy(realloc_bdata(b, b->core.l_qname), qname, b->core.l_qname);
  set_cigar(cigar); //this will also set core.bin
  add_sequence(qseq, l_qseq);
  add_quals(quals); //quals must be given as Phred33
  set_flags(flags);
  set_mdata(mg_tid, (int32_t)(mate_pos-1), (int32_t)insert_size);
  if (aux_strings!=NULL) {
    for (int i=0;i<aux_strings->Count();i++) {
       add_aux(aux_strings->Get(i));
       }
    }
}

 void GBamRecord::set_cigar(const char* cigar) {
   //requires b->core.pos and b->core.flag to have been set properly PRIOR to this call
   int doff=b->core.l_qname;
   uint8_t* after_cigar=NULL;
   int after_cigar_len=0;
   uint8_t* prev_bdata=NULL;
   if (b->data_len>doff) {
      //cigar string already allocated, replace it
      int d=b->core.l_qname + b->core.n_cigar * 4;//offset of after-cigar data
      after_cigar=b->data+d;
      after_cigar_len=b->data_len-d;
      }
   const char *s;
   char *t;
   int i, op;
   long x;
   b->core.n_cigar = 0;
   if (cigar != NULL && strcmp(cigar, "*") != 0) {
        for (s = cigar; *s; ++s) {
            if (isalpha(*s)) b->core.n_cigar++;
            else if (!isdigit(*s)) {
                 GError("Error: invalid CIGAR character (%s)\n",cigar);
                 }
            }
        if (after_cigar_len>0) { //replace/insert into existing full data
             prev_bdata=dupalloc_bdata(b, doff + b->core.n_cigar * 4 + after_cigar_len);
             memcpy((void*)(b->data+doff+b->core.n_cigar*4),(void*)after_cigar, after_cigar_len);
             free(prev_bdata);
             }
           else {
             realloc_bdata(b, doff + b->core.n_cigar * 4);
             }
        for (i = 0, s = cigar; i != b->core.n_cigar; ++i) {
            x = strtol(s, &t, 10);
            op = toupper(*t);
            if (op == 'M' || op == '=' || op == 'X') op = BAM_CMATCH;
            else if (op == 'I') op = BAM_CINS;
            else if (op == 'D') op = BAM_CDEL;
            else if (op == 'N') op = BAM_CREF_SKIP;
            else if (op == 'S') op = BAM_CSOFT_CLIP;
            else if (op == 'H') op = BAM_CHARD_CLIP;
            else if (op == 'P') op = BAM_CPAD;
            else GError("Error: invalid CIGAR operation (%s)\n",cigar);
            s = t + 1;
            bam1_cigar(b)[i] = x << BAM_CIGAR_SHIFT | op;
        }
        if (*s) GError("Error: unmatched CIGAR operation (%s)\n",cigar);
        b->core.bin = bam_reg2bin(b->core.pos, bam_calend(&b->core, bam1_cigar(b)));
    } else {//no CIGAR string given
        if (!(b->core.flag&BAM_FUNMAP)) {
            GMessage("Warning: mapped sequence without CIGAR (%s)\n", (char*)b->data);
            b->core.flag |= BAM_FUNMAP;
        }
        b->core.bin = bam_reg2bin(b->core.pos, b->core.pos + 1);
    }
   setupCoordinates();
   } //set_cigar()

 void GBamRecord::add_sequence(const char* qseq, int slen) {
   //must be called AFTER set_cigar (cannot replace existing sequence for now)
   if (qseq==NULL) return; //should we ever care about this?
   if (slen<0) slen=strlen(qseq);
   int doff = b->core.l_qname + b->core.n_cigar * 4;
   if (strcmp(qseq, "*")!=0) {
       b->core.l_qseq=slen;
       if (b->core.n_cigar && b->core.l_qseq != (int32_t)bam_cigar2qlen(&b->core, bam1_cigar(b)))
           GError("Error: CIGAR and sequence length are inconsistent!(%s)\n",
                  qseq);
       uint8_t* p = (uint8_t*)realloc_bdata(b, doff + (b->core.l_qseq+1)/2 + b->core.l_qseq) + doff;
       //also allocated quals memory
       memset(p, 0, (b->core.l_qseq+1)/2);
       for (int i = 0; i < b->core.l_qseq; ++i)
           p[i/2] |= bam_nt16_table[(int)qseq[i]] << 4*(1-i%2);
       } else b->core.l_qseq = 0;
   }

 void GBamRecord::add_quals(const char* quals) {
   //requires core.l_qseq already set
   //and must be called AFTER add_sequence(), which also allocates the memory for quals
   uint8_t* p = b->data+(b->core.l_qname + b->core.n_cigar * 4 + (b->core.l_qseq+1)/2);
   if (quals==NULL || strcmp(quals, "*") == 0) {
      for (int i=0;i < b->core.l_qseq; i++) p[i] = 0xff;
      return;
      }
   for (int i=0;i < b->core.l_qseq; i++) p[i] = quals[i]-33;
   }

 void GBamRecord::add_aux(const char* str) {
     //requires: being called AFTER add_quals()
     int strl=strlen(str);
     //int doff = b->core.l_qname + b->core.n_cigar*4 + (b->core.l_qseq+1)/2 + b->core.l_qseq + b->l_aux;
     //int doff0=doff;
     if (strl < 6 || str[2] != ':' || str[4] != ':')
         parse_error("missing colon in auxiliary data");
     tag[0] = str[0]; tag[1] = str[1];
     uint8_t atype = str[3];
     uint8_t* adata=abuf;
     int alen=0;
     if (atype == 'A' || atype == 'a' || atype == 'c' || atype == 'C') { // c and C for backward compatibility
         atype='A';
         alen=1;
         adata=(uint8_t*)&str[5];
         }
      else if (atype == 'I' || atype == 'i') {

		 long long x=strtoll(str+5, NULL, 10); //(long long)atoll(str + 5);
		 //long x=(long)atol(str + 5);
		 if (x < 0) {
             if (x >= -127) {
                 atype='c';
                 abuf[0] =  (int8_t)x;
                 alen=1;
                 }
             else if (x >= -32767) {
                 atype = 's';
                 *(int16_t*)abuf = (int16_t)x;
                 alen=2;
                 }
             else {
                 atype='i';
                 *(int32_t*)abuf = (int32_t)x;
                 alen=4;
                 if (x < -2147483648ll)
                     GMessage("Parse warning: integer %lld is out of range.", x);
                 }
             } else { //x >=0
             if (x <= 255) {
                 atype = 'C';
                 abuf[0] = (uint8_t)x;
                 alen=1;
                 }
             else if (x <= 65535) {
                 atype='S';
                 *(uint16_t*)abuf = (uint16_t)x;
                 alen=2;
                 }
             else {
                 atype='I';
                 *(uint32_t*)abuf = (uint32_t)x;
                 alen=4;
                 if (x > 4294967295ll)
                     GMessage("Parse warning: integer %lld is out of range.", x);
                 }
             }
         } //integer type
         else if (atype == 'f') {
             *(float*)abuf = (float)atof(str + 5);
             alen = sizeof(float);
             }
         else if (atype == 'd') { //?
             *(float*)abuf = (float)atof(str + 9);
             alen=8;
             }
         else if (atype == 'Z' || atype == 'H') {
             if (atype == 'H') { // check whether the hex string is valid
                 if ((strl - 5) % 2 == 1) parse_error("length of the hex string not even");
                 for (int i = 0; i < strl - 5; ++i) {
                     int c = toupper(str[5 + i]);
                     if (!((c >= '0' && c <= '9') || (c >= 'A' && c <= 'F')))
                         parse_error("invalid hex character");
                     }
                 }
             memcpy(abuf, str + 5, strl - 5);
             abuf[strl-5] = 0;
             alen=strl-4;
             } else parse_error("unrecognized aux type");
  this->add_aux(tag, atype, alen, adata);
  }//add_aux()

void interpret_CIGAR(char cop, int cl, int aln_start) {
// NB: Closed ranges
// gpos = current genomic position (will end up as right coordinate on the genome)
// rpos = read position (will end up as the length of the read)
// cop = CIGAR operation, cl = operation length
int rpos = 0;
int gpos = aln_start;
int num_mismatches=0; //NM tag value = edit distance
switch (cop) {
 case BAM_CDIFF:  // X
      num_mismatches+=cl;
 case BAM_CMATCH: // M
      //have to actually check for mismatches: num_mismatches+=count_mismatches;
 case BAM_CEQUAL: // = 
      //printf("[%d-%d]", gpos, gpos + cl - 1);
      gpos+=cl;
      rpos+=cl;
      break;
 case BAM_CPAD:
      // printf("[%d-%d]", pos, pos + cl - 1); // Spans positions, No Coverage
      gpos+=cl;
      break;

 case BAM_CHARD_CLIP:
      // printf("[%d]", pos);  // No coverage 
      // gpos is not advanced by this operation
      break;

 case BAM_CSOFT_CLIP: // S
      //soft clipped bases, present in SEQ
      rpos+=cl;
      break;
 case BAM_CINS: // I
      // No Coverage
      // adds cl bases "throughput" but not genomic position "coverage" (gap in genomic seq)
      // should also add cl to the number of "mismatches" (unaligned read bases)
      num_mismatches+=cl;
      // How you handle this is application dependent 
      // gpos is not advanced by this operation 
      rpos+=cl;
      break;

 case BAM_CDEL: // D
      //deletion in reference sequence relative to the read (gap in read sequence)
      // printf("[%d-%d]", pos, pos + cl - 1);  
      // Spans positions
      num_mismatches+=cl;
      gpos += cl;
      break;

 case BAM_CREF_SKIP: // N
      // intron
      //special skip operation, not contributing to "edit distance", 
      // printf("[%d-%d]", pos, pos + cl - 1); // Spans positions, No Coverage 
      //   so num_mismatches is not updated
      gpos+=cl;
      break;

 default:
      fprintf(stderr, "Unhandled cigar_op %d:%d\n", cop, cl);
      //printf("?");
  }
} // interpret_CIGAR(), just a reference of CIGAR operations interpretation

 void GBamRecord::setupCoordinates() {
   uint32_t *cigar = bam1_cigar(b);
   const bam1_core_t *c = &b->core;
   if (c->flag & BAM_FUNMAP) return; /* skip unmapped reads */
   int l=0;
   start=c->pos+1; //genomic start coordinate, 1-based (BAM core.pos is 0-based)
   int exstart=c->pos;
   for (int i = 0; i < c->n_cigar; ++i) {
        int op = cigar[i]&0xf;
        if (op == BAM_CMATCH || op==BAM_CEQUAL ||
            op == BAM_CDIFF || op == BAM_CDEL) {
               l += cigar[i]>>4;
               }
           else if (op == BAM_CREF_SKIP) { //N
               //intron starts
               //exon ends here
               GSeg exon(exstart+1,c->pos+l);
               exons.Add(exon);
               l += cigar[i]>>4;
               exstart=c->pos+l;
               }
        }
   GSeg exon(exstart+1,c->pos+l);
   exons.Add(exon);
   end=c->pos+l; //genomic start coordinate
 }
 /*
 int GBamRecord::find_tag(const char tag[2], uint8_t* & s, char& tag_type) {
   //position s at the beginning of tag "data" (after the type char)
   //returns the length of tag data, and tag type in tag_type
   uint8_t* aux_start=bam1_aux(b);
   s=aux_start;
   while (s < aux_start + b->l_aux - 1) {
     char key[2];
     key[0] = (char)s[0]; key[1] = (char)s[1];
     s += 2; tag_type = (char)*s; ++s;
     int inc=0;
     int m_inc=0;       //for 'B' case
     uint8_t sub_type=0; // --,,--
     switch (tag_type) {
       case 'A':
       case 'C':
       case 'c':
         inc=1;
         break;
       case 'S':
       case 's':
         inc=2;
         break;
       case 'I':
       case 'i':
       case 'f':
         inc=4;
         break;
       case 'd':
         inc=8;
         break;
       case 'B':
           sub_type = *(s+1);
           int32_t n;
           memcpy(&n, s+1, 4);
           inc += 5;
           //kputc(type, &str); kputc(':', &str); kputc(sub_type, &str);
           m_inc=0;
           if (sub_type == 'c' || sub_type == 'C')
                        { m_inc=1; }
           else if (sub_type == 's' || sub_type == 'S')
                         { m_inc = 2; }
           else if ('i' == sub_type || 'I' == sub_type || 'f' == sub_type)
                         { m_inc = 4; }
           if (m_inc==0)
              GError("Error: invalid 'B' array subtype (%c)!\n",sub_type);
           inc += m_inc*n;
           break;
       case 'H':
       case 'Z':
         while (*(s+inc)) ++inc;
          ++inc; // to skip the terminating \0
         break;
       } //switch (tag_type)
     if (tag[0]==key[0] && tag[1]==key[1]) {
        return inc;
        }
     s+=inc;
     }//while aux data
   return 0;
   }
*/
 uint8_t* GBamRecord::find_tag(const char tag[2]) {
   return bam_aux_get(this->b, tag);
 }

 char GBamRecord::tag_char(const char tag[2]) { //retrieve tag data as single char
   uint8_t* s=find_tag(tag);
   if (s) return ( bam_aux2A(s) );
   return 0;
  }

 int GBamRecord::tag_int(const char tag[2]) { //get the numeric value of tag
   uint8_t *s=find_tag(tag);
   if (s) return ( bam_aux2i(s) );
   return 0;
   }

 char* GBamRecord::tag_str(const char tag[2]) { //return string value for a tag
   uint8_t *s=find_tag(tag);
   if (s) return ( bam_aux2Z(s) );
   return NULL;
   }

 char GBamRecord::spliceStrand() { // '+', '-' from the XS tag, or 0 if no XS tag
   char c=tag_char("XS");
   if (c) return c;
     else return '.';
   }

 char* GBamRecord::sequence() { //user must free this after use
   char *s = (char*)bam1_seq(b);
   char* qseq=NULL;
   GMALLOC(qseq, (b->core.l_qseq+1));
   int i;
   for (i=0;i<(b->core.l_qseq);i++) {
     int8_t v = bam1_seqi(s,i);
     qseq[i] = bam_nt16_rev_table[v];
     }
   qseq[i] = 0;
   return qseq;
   }

 char* GBamRecord::qualities() {//user must free this after use
   char *qual  = (char*)bam1_qual(b);
   char* qv=NULL;
   GMALLOC(qv, (b->core.l_qseq+1) );
   int i;
   for(i=0;i<(b->core.l_qseq);i++) {
     qv[i]=qual[i]+33;
     }
   qv[i]=0;
   return qv;
   }

 char* GBamRecord::cigar() { //returns text version of the CIGAR string; must be freed by user
   kstring_t str;
   str.l = str.m = 0; str.s = 0;
   if (b->core.n_cigar == 0) kputc('*', &str);
    else {
      for (int i = 0; i < b->core.n_cigar; ++i) {
         kputw(bam1_cigar(b)[i]>>BAM_CIGAR_SHIFT, &str);
         kputc("MIDNSHP=X"[bam1_cigar(b)[i]&BAM_CIGAR_MASK], &str);
         }
      }
   return str.s;
   }
