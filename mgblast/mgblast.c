#define MGBLAST_OPTS 1

#include <objsset.h>
#include <blast.h>
#include <txalign.h>
#include <mblast.h>
#include <xmlblast.h>
#include <sqnutils.h>
#include <blfmtutl.h>
#include "mgblast_util.h"

#define DEFLINE_BUF 255


/* Used by the callback function. */
FILE *global_fp=NULL;
/*
    Callback to print out ticks, in UNIX only due to file systems
    portability issues.
*/
static int LIBCALLBACK
dummy_callback(Int4 sequence_number, Int4 number_of_positive_hits)
{
   return 0;
}

static int LIBCALLBACK
tick_callback(Int4 sequence_number, Int4 number_of_positive_hits)

{

#ifdef OS_UNIX
   
   fprintf(global_fp, "%s", ".");
   fflush(global_fp);
#endif
    return 0;
}

enum {
   MBLAST_ENDPOINTS = 0,
   MBLAST_SEGMENTS,
   MBLAST_ALIGNMENTS,
   MBLAST_ALIGN_INFO,
   #ifdef MGBLAST_OPTS
    MGBLAST_FLTHITS, /* tab output */
    MGBLAST_HITGAPS, /* tab output with gap info */
   #endif
   MBLAST_DELAYED_TRACEBACK
};

/* 
{Geo:} MGBLAST_FLTHITS (-D4) will print and filter compact output by max_overhang, min_overlap :

 q_name, q_len, q_n5, q_n3, hit_name, hit_len, hit_n5, hit_n3, pid, score, p-value, strand 

 for MBLAST_HITGAPS option we have the same output but with 2 extra fields
 
 <qgaps> <hitgaps>
 
 where each of them has the gap coordinate and gap length, like this:
 
 gap1pos_gap1len[+gap1len],gap2pos[+gap2len],...
 
*/

int startsWith(char* s, const char* prefix) {
 int i=0;
 if (prefix==NULL || s==NULL) return 0;
 while (prefix[i]!='\0' && prefix[i]==s[i]) i++;
 return (prefix[i]=='\0');
 };

int isET(char* name) {
 char buf[7];
 int i;
 for (i=0;i<6 && name[i]!='\0';i++) 
    buf[i]=tolower(name[i]);
 buf[i]='\0';   
 return (startsWith(buf, "np|") ||
        startsWith(buf, "et|") || 
        startsWith(buf, "egad|") ||
        startsWith(buf, "full|") ||
        startsWith(buf, "mrna|")
        ) ? 1 : 0; 
 };


static int LIBCALLBACK
MegaBlastPrintEndpoints(VoidPtr ptr)
{
   BlastSearchBlkPtr search = (BlastSearchBlkPtr) ptr;
   CharPtr subject_descr;
   SeqIdPtr sip, query_id;
   CharPtr query_buffer, title;
   CharPtr subject_buffer;
   Int4 query_length, q_start, q_end, q_shift=0, s_shift=0;
   Int4 subject_end;
   Int4 hsp_index;
   Boolean numeric_sip_type = FALSE;
   BLAST_HSPPtr hsp; 
   Int2 context;
   Char context_sign;
   Int4 subject_gi, score;
   FILE *fp = (FILE *) search->output;

   if (search->current_hitlist == NULL || search->current_hitlist->hspcnt <= 0) {
      search->subject_info = BLASTSubjectInfoDestruct(search->subject_info);
      return 0;
   }

   if (search->rdfp)
      readdb_get_descriptor(search->rdfp, search->subject_id, &sip,
                            &subject_descr);
   else 
      sip = SeqIdSetDup(search->subject_info->sip);
   
   if (sip->choice != SEQID_GENERAL ||
       StringCmp(((DbtagPtr)sip->data.ptrvalue)->db, "BL_ORD_ID")) {
      if (search->pbp->mb_params->full_seqids) {
         subject_buffer = (CharPtr) Malloc(BUFFER_LENGTH + 1);
         SeqIdWrite(sip, subject_buffer, PRINTID_FASTA_LONG, BUFFER_LENGTH);
      } else
         numeric_sip_type = GetAccessionFromSeqId(SeqIdFindBest(sip, SEQID_GI), 
                                                  &subject_gi, &subject_buffer);
   } else {
      DbtagPtr db_tag = (DbtagPtr) sip->data.ptrvalue;
      if (db_tag->db && 
          (!StringCmp(db_tag->db, "THC") || 
           !StringICmp(db_tag->db, "TI")) && 
          db_tag->tag->id != 0) {
         subject_buffer = (CharPtr) Malloc(16);
         sprintf(subject_buffer, "%ld", (long) db_tag->tag->id);
      } else {
         subject_buffer = StringTokMT(subject_descr, " \t", &subject_descr);
         subject_descr = subject_buffer;
      }
   }

   search->current_hitlist->hspcnt_max = search->current_hitlist->hspcnt;

   /* Only for the two sequences case, get offset shift if subject 
      is a subsequence */
   if (!search->rdfp && search->query_slp->next) {
       s_shift = SeqLocStart(search->query_slp->next);
       subject_end = SeqLocStop(search->query_slp->next);
   } else {
      s_shift = 0;
      subject_end = 
         readdb_get_sequence_length(search->rdfp, search->subject_id);
   }
   /* Get offset shift if query is a subsequence */
   q_shift = SeqLocStart(search->query_slp);

   for (hsp_index=0; hsp_index<search->current_hitlist->hspcnt; hsp_index++) {
      hsp = search->current_hitlist->hsp_array[hsp_index];
      if (hsp==NULL || (search->pbp->cutoff_e > 0 && 
      hsp->evalue > search->pbp->cutoff_e)) 
     continue;
      
      /* Correct query context is already found in BlastGetNonSumStatsEvalue */
      context = hsp->context; 
      query_id = search->qid_array[context/2];


      if (query_id == NULL) /* Bad hsp, something wrong */
     continue; 
      hsp->context = context & 1;      
      query_length = search->query_context_offsets[context+1] -
         search->query_context_offsets[context] - 1;
      hsp->subject.end = hsp->subject.offset + hsp->subject.length;

      if (hsp->context) {
     hsp->query.end = query_length - hsp->query.offset;
     hsp->query.offset = 
        hsp->query.end - hsp->query.length + 1;
     context_sign = '-'; 
      } else {
     hsp->query.end = (++hsp->query.offset) + hsp->query.length - 1;
         if (hsp->query.end > query_length) {
            hsp->subject.end -= (hsp->query.end - query_length);
            hsp->query.end = query_length;
         }
     context_sign = '+';  
      }
      
      if (hsp->subject.end > subject_end) {
         hsp->query.end -= (hsp->subject.end - subject_end);
         hsp->subject.end = subject_end;
      }
      hsp->subject.offset++;
      
      query_buffer = NULL;
      if (query_id->choice == SEQID_LOCAL && 
          search->pbp->mb_params->full_seqids) {
         BioseqPtr query_bsp = BioseqLockById(query_id);
         title = StringSave(BioseqGetTitle(query_bsp));
         if (title)
            query_buffer = StringTokMT(title, " ", &title);
         else {
            Int4 query_gi;
            GetAccessionFromSeqId(query_bsp->id, &query_gi,
                                  &query_buffer);
         }  
         BioseqUnlock(query_bsp);
      } else {
         query_buffer = (CharPtr) Malloc(BUFFER_LENGTH + 1);
         if (!search->pbp->mb_params->full_seqids)
            SeqIdWrite(query_id, query_buffer, PRINTID_TEXTID_ACCESSION,
                       BUFFER_LENGTH);
         else 
            SeqIdWrite(query_id, query_buffer, PRINTID_FASTA_LONG,
                    BUFFER_LENGTH);
      }

      if (search->pbp->gap_open==0 && search->pbp->gap_extend==0)
     score = ((hsp->subject.length + hsp->query.length)*
           search->sbp->reward / 2 - hsp->score) / 
        (search->sbp->reward - search->sbp->penalty);
      else 
     score = hsp->score;

      if (context_sign == '+') {
     q_start = hsp->query.offset;
     q_end = hsp->query.end;
      } else {
     q_start = hsp->query.end;
     q_end = hsp->query.offset;
      }

      /* Adjust offsets if query is a subsequence, only for first query */
      if (context < 2) {
          q_start += q_shift;
          q_end += q_shift;
      }

      hsp->subject.offset += s_shift;
      hsp->subject.end += s_shift;

      if (numeric_sip_type)
     fprintf(fp, "'%ld'=='%c%s' (%d %d %d %d) %d\n", (long) subject_gi, 
         context_sign, query_buffer, hsp->subject.offset, q_start, 
         hsp->subject.end, q_end, score);
      else 
     fprintf(fp, "'%s'=='%c%s' (%d %d %d %d) %d\n", 
         subject_buffer, context_sign, query_buffer, 
         hsp->subject.offset, q_start, 
         hsp->subject.end, q_end, score);
      MemFree(query_buffer);
   }
   if (!numeric_sip_type && subject_buffer != subject_descr)
      MemFree(subject_buffer);
   MemFree(subject_descr);
   sip = SeqIdSetFree(sip);
   return 0;
}

static int LIBCALLBACK
MegaBlastPrintSegments(VoidPtr ptr)
{
   BlastSearchBlkPtr search = (BlastSearchBlkPtr) ptr;
   ReadDBFILEPtr rdfp = search->rdfp;
   BLAST_HSPPtr hsp; 
   Int4 i, subject_gi;
   Int2 context;
   CharPtr query_buffer, title;
   SeqIdPtr sip, query_id; 
   Int4 hsp_index, score;
   Uint1Ptr query_seq, subject_seq = NULL;
   FloatHi perc_ident;
   Char strand;
   GapXEditScriptPtr esp;
   Int4 q_start, q_end, s_start, s_end, query_length, numseg;
   Int4 q_off, num_ident, align_length, total_ident, q_shift=0, s_shift=0;
   Int4Ptr length, start;
   Uint1Ptr strands;
   CharPtr subject_descr, subject_buffer, buffer;
   Char tmp_buffer[BUFFER_LENGTH];
   Int4 buffer_size, max_buffer_size = LARGE_BUFFER_LENGTH;
   Boolean numeric_sip_type = FALSE;
   FILE *fp = (FILE *) search->output;

   if (search->current_hitlist == NULL || search->current_hitlist->hspcnt <= 0) {
      search->subject_info = BLASTSubjectInfoDestruct(search->subject_info);
      return 0;
   }

   subject_seq = search->subject->sequence_start + 1;


   if (rdfp)
      readdb_get_descriptor(rdfp, search->subject_id, &sip, &subject_descr);
   else 
      sip = SeqIdSetDup(search->subject_info->sip);

   if (sip->choice != SEQID_GENERAL ||
       StringCmp(((DbtagPtr)sip->data.ptrvalue)->db, "BL_ORD_ID")) {
      if (search->pbp->mb_params->full_seqids) { 
         subject_buffer = (CharPtr) Malloc(BUFFER_LENGTH + 1);
         SeqIdWrite(sip, subject_buffer, PRINTID_FASTA_LONG, BUFFER_LENGTH);
      } else
         numeric_sip_type = GetAccessionFromSeqId(SeqIdFindBest(sip, SEQID_GI), 
                                                  &subject_gi, &subject_buffer);
   } else {
      subject_buffer = StringTokMT(subject_descr, " \t", &subject_descr);
      subject_descr = subject_buffer;
   }

   buffer = (CharPtr) Malloc(LARGE_BUFFER_LENGTH);

   /* Only for the two sequences case, get offset shift if subject 
      is a subsequence */
   if (!rdfp && search->query_slp->next)
       s_shift = SeqLocStart(search->query_slp->next);
   /* Get offset shift if query is a subsequence */
   q_shift = SeqLocStart(search->query_slp);
   
   for (hsp_index=0; hsp_index<search->current_hitlist->hspcnt; hsp_index++) {
      hsp = search->current_hitlist->hsp_array[hsp_index];
      if (hsp==NULL || (search->pbp->cutoff_e > 0 && 
                        hsp->evalue > search->pbp->cutoff_e)) {
     continue;
      }
      context = hsp->context;
      query_id = search->qid_array[context/2];

     
      if (query_id == NULL) /* Bad hsp, something wrong */
     continue; 
      hsp->context = context & 1;

      if (search->pbp->gap_open==0 && search->pbp->gap_extend==0)
     score = ((hsp->subject.length + hsp->query.length)*
           search->sbp->reward / 2 - hsp->score) / 
        (search->sbp->reward - search->sbp->penalty);
      else 
     score = hsp->score;

      query_length = search->query_context_offsets[context+1] -
         search->query_context_offsets[context] - 1;

      q_off = hsp->query.offset;
      if (hsp->context) {
     strand = '-'; 
     hsp->query.end = query_length - hsp->query.offset;
     hsp->query.offset = 
        hsp->query.end - hsp->query.length;
      } else {
     strand = '+';  
     hsp->query.end = hsp->query.offset + hsp->query.length;
      }

      if (strand == '+') {
     q_start = hsp->query.offset + 1;
     q_end = hsp->query.end;
      } else {
     q_start = hsp->query.end;
     q_end = hsp->query.offset + 1;
      }
      s_start = hsp->subject.offset + 1;
      s_end = hsp->subject.offset + hsp->subject.length;

      /* Adjust offsets if query is a subsequence, only for first query */
      if (context < 2) {
          q_start += q_shift;
          q_end += q_shift;
      }

      s_start += s_shift;
      s_end += s_shift;

      if (query_id->choice == SEQID_LOCAL && 
          search->pbp->mb_params->full_seqids) {
         BioseqPtr query_bsp = BioseqLockById(query_id);
         title = StringSave(BioseqGetTitle(query_bsp));
         if (title)
            query_buffer = StringTokMT(title, " ", &title);
         else {
            Int4 query_gi;
            GetAccessionFromSeqId(query_bsp->id, &query_gi,
                                  &query_buffer);
         }  
         BioseqUnlock(query_bsp);
      } else {
         query_buffer = (CharPtr) Malloc(BUFFER_LENGTH + 1);
         if (!search->pbp->mb_params->full_seqids)
            SeqIdWrite(query_id, query_buffer, PRINTID_TEXTID_ACCESSION,
                       BUFFER_LENGTH);
         else 
            SeqIdWrite(query_id, query_buffer, PRINTID_FASTA_LONG,
                    BUFFER_LENGTH);
      }

      if (numeric_sip_type)
     sprintf(buffer, "\n#'>%ld'=='%c%s' (%d %d %d %d) %d\na {\n  s %d\n  b %d %d\n  e %d %d\n", 
          (long) subject_gi, strand, query_buffer, 
          s_start, q_start, s_end, q_end, score, score, 
          s_start, q_start, s_end, q_end);
      else 
     sprintf(buffer, "\n#'>%s'=='%c%s' (%d %d %d %d) %d\na {\n  s %d\n  b %d %d\n  e %d %d\n", 
          subject_buffer, strand, query_buffer, 
          s_start, q_start, s_end, q_end, score, score, 
          s_start, q_start, s_end, q_end);
      buffer_size = StringLen(buffer);

      query_seq = search->context[context].query->sequence;

      esp = hsp->gap_info->esp;
        
      for (numseg=0; esp; esp = esp->next, numseg++);

      GXECollectDataForSeqalign(hsp->gap_info, hsp->gap_info->esp, numseg,
                &start, &length, &strands, 
                &q_off, &hsp->subject.offset);

      if (start[0] < 0) {
         length[0] += start[0];
         start[1] -= start[0];
         start[0] = 0;
      } 
      if (start[2*(numseg-1)] + length[numseg-1] > query_length) 
         length[numseg-1] = query_length - start[2*(numseg-1)];
      
      total_ident = 0;
      align_length = 0;
      for (i=0; i<numseg; i++) {
         align_length += length[i];
     if (strand == '+') {
        q_start = start[2*i] + 1;
        q_end = q_start + length[i] - 1;
     } else {
        q_start = query_length - start[2*i];
        q_end = q_start - length[i] + 1;
     }
     if (start[2*i] != -1 && start[2*i+1] != -1) {
        num_ident = MegaBlastGetNumIdentical(query_seq, subject_seq, 
                                                 start[2*i], start[2*i+1], 
                                                 length[i], FALSE);
            perc_ident = (FloatHi) num_ident / length[i] * 100;
            total_ident += num_ident;
        sprintf(tmp_buffer, "  l %d %d %d %d (%.0f)\n", start[2*i+1]+1, 
            q_start, start[2*i+1]+length[i],
            q_end, perc_ident);     
        if ((buffer_size += StringLen(tmp_buffer)) > max_buffer_size - 2) {
           max_buffer_size *= 2;
           buffer = (CharPtr) Realloc(buffer, max_buffer_size);
        }
        StringCat(buffer, tmp_buffer);
     }
      }
      if (100*total_ident >= 
          align_length*search->pbp->mb_params->perc_identity) {
        StringCat(buffer, "}");
        fprintf(fp, "%s\n", buffer);
      }
      MemFree(start);
      MemFree(length);
      MemFree(strands);
      MemFree(query_buffer);
   } /* End loop on hsp's */
   if (!numeric_sip_type && subject_buffer != subject_descr)
      MemFree(subject_buffer);
   MemFree(subject_descr);
   MemFree(buffer);
   sip = SeqIdSetFree(sip);
   fflush(fp);
   return 1;
}

/* Breaks up a location like "2000 3000" into two integers 
   that are returned.

   If location is NULL then the integers are set to 0.
*/

static Boolean
Megablast_GetLoc(char* location, Int4* start, Int4* end)
{
        CharPtr delimiters = " ,;";

        if (start == NULL || end == NULL)
           return FALSE;

        *start = 0;
        *end = 0;

        if (location == NULL)
           return TRUE;

        *start =  atoi(StringTokMT(location, delimiters, &location));
        *end = atoi(location);

        return TRUE;
}

typedef enum {
ARG_DB = 0,
ARG_QUERY,
ARG_EVALUE,
ARG_FORMAT,
ARG_OUT,
ARG_FILTER,
ARG_XDROP,
ARG_SHOWGIS,
ARG_MISMATCH,
ARG_MATCH,
ARG_DESCRIPTIONS,
ARG_ALIGNMENTS,
ARG_OUTTYPE,
ARG_THREADS,
ARG_ASNOUT,
ARG_BELIEVEQUERY,
ARG_MAXQUERY,
ARG_WORDSIZE,
ARG_DBSIZE,
ARG_SEARCHSP,
ARG_MAXPOS,
ARG_STRAND,
ARG_HTML,
ARG_GILIST,
ARG_GAPOPEN,
ARG_GAPEXT,
ARG_MINSCORE,
ARG_MASKEDQUERY,
ARG_FULLID,
ARG_LCASE,
ARG_LOGINFO,
ARG_PERC_IDENT,
ARG_QUERYLOC,
ARG_WINDOW,
ARG_XDROP_UNGAPPED,
ARG_XDROP_FINAL,
ARG_TEMPL_LEN,
ARG_EVERYBASE,
ARG_DYNAMIC,
ARG_TEMPL_TYPE,
ARG_MAXHSP,
ARG_FORCE_OLD
} BlastArguments;

#define DO_NOT_SUPPRESS_BLAST_OP

#define NUMARG (sizeof(myargs)/sizeof(myargs[0]))

static Args myargs [] = {
  { "Database", 
    "nr", NULL, NULL, FALSE, 'd', ARG_STRING, 0.0, 0, NULL},       /* ARG_DB */
  { "Query File", 
    NULL, NULL, NULL, FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL}, /* ARG_QUERY */
  { "Expectation value", 
    "10.0", NULL, NULL, FALSE, 'e', ARG_FLOAT, 0.0, 0, NULL},/* ARG_EVALUE */
  { "alignment view options:\n0 = pairwise,\n1 = query-anchored showing identities,"
         "\n2 = query-anchored no identities,\n3 = flat query-anchored, show identities,"
         "\n4 = flat query-anchored, no identities,\n5 = query-anchored no identities and blunt ends,"
         "\n6 = flat query-anchored, no identities and blunt ends,"
         "\n7 = XML Blast output,\n8 = tabular, \n9 tabular with comment lines,"
         "\n10 ASN, text\n11 ASN, binary\n12 mgblast tab (default)\n13 mgblast tab with gapinfo", 
        "12", "0", "13", FALSE, 'm', ARG_INT, 0.0, 0, NULL},       /* ARG_FORMAT */
  { "BLAST report Output File", 
    "stdout", NULL, NULL, TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},/* ARG_OUT */
  { "Filter query sequence",
        "T", NULL, NULL, FALSE, 'F', ARG_STRING, 0.0, 0, NULL},    /* ARG_FILTER */
  { "X dropoff value for gapped alignment (in bits)",
    "20", NULL, NULL, FALSE, 'X', ARG_INT, 0.0, 0, NULL},      /* ARG_XDROP */
  { "Show GI's in deflines",
        "F", NULL, NULL, FALSE, 'I', ARG_BOOLEAN, 0.0, 0, NULL},   /* ARG_SHOWGIS */
  { "Penalty for a nucleotide mismatch",
    "-3", NULL, NULL, FALSE, 'q', ARG_INT, 0.0, 0, NULL},      /* ARG_MISMATCH */
  { "Reward for a nucleotide match",
    "1", NULL, NULL, FALSE, 'r', ARG_INT, 0.0, 0, NULL},       /* ARG_MATCH */
  { "Number of database sequences to show one-line descriptions for (V)",
        "500", NULL, NULL, FALSE, 'v', ARG_INT, 0.0, 0, NULL},     /* ARG_DESCRIPTIONS */
  { "Number of database sequence to show alignments for (B)",
        "250", NULL, NULL, FALSE, 'b', ARG_INT, 0.0, 0, NULL},     /* ARG_ALIGNMENTS */
  { "Type of output:\n0 - alignment endpoints and score\n1 - all ungapped segments endpoints"
     "\n2 - traditional BLAST output\n3 - tab-delimited one line format"
     "\n4 - filtered tabulated hits (default)\n5 - tabulated hits with gap info"
     /* "\n6 - incremental text ASN.1\n7 - incremental binary ASN.1" */
       ,
        "4", NULL, NULL, FALSE, 'D', ARG_INT, 0.0, 0, NULL},       /* ARG_OUTTYPE */
  { "Number of processors to use",
        "1", NULL, NULL, FALSE, 'a', ARG_INT, 0.0, 0, NULL},       /* ARG_THREADS */
  { "ASN.1 SeqAlign file; must be used in conjunction with -D2 option", 
    NULL, NULL, NULL, TRUE, 'O', ARG_FILE_OUT, 0.0, 0, NULL},  /* ARG_ASNOUT */
  { "Believe the query defline",
        "F", NULL, NULL, TRUE, 'J', ARG_BOOLEAN, 0.0, 0, NULL},/* ARG_BELIEVEQUERY */
  { "Maximal total length of queries for a single search",
        "5000000", NULL, NULL, FALSE, 'M', ARG_INT, 0.0, 0, NULL},/* ARG_MAXQUERY */
  { "Word size (length of best perfect match)", 
        "28", NULL, NULL, FALSE, 'W', ARG_INT, 0.0, 0, NULL},      /* ARG_WORDSIZE */
  { "Effective length of the database (use zero for the real size)", 
        "0", NULL, NULL, FALSE, 'z', ARG_FLOAT, 0.0, 0, NULL},     /* ARG_DBSIZE */
  { "Effective length of the search space (use zero for the real size)", 
        "0", NULL, NULL, FALSE, 'Y', ARG_FLOAT, 0.0, 0, NULL},     /* ARG_SEARCHSP */
  { "Maximal number of positions for a hash value (set to 0 to ignore)",
        "0", NULL, NULL, FALSE, 'P', ARG_INT, 0.0, 0, NULL},       /* ARG_MAXPOS */
  { "Query strands to search against database: 3 is both, 1 is top, 2 is bottom",
        "3", NULL, NULL, FALSE, 'S', ARG_INT, 0.0, 0, NULL},       /* ARG_STRAND */
  { "Produce HTML output",
        "F", NULL, NULL, FALSE, 'T', ARG_BOOLEAN, 0.0, 0, NULL},   /* ARG_HTML */
  { "Restrict search of database to list of GI's",
    NULL, NULL, NULL, TRUE, 'l', ARG_STRING, 0.0, 0, NULL},    /* ARG_GILIST */
  { "Cost to open a gap (-1 invokes default behavior)",          /* ARG_GAPOPEN */
        "-1", NULL, NULL, FALSE, 'G', ARG_INT, 0.0, 0, NULL},
  { "Cost to extend a gap (-1 invokes default behavior)",        /* ARG_GAPEXT */
        "-1", NULL, NULL, FALSE, 'E', ARG_INT, 0.0, 0, NULL},
  { "Minimal hit score to report (0 for default behavior)",        /* ARG_MINSCORE */
        "0", NULL, NULL, FALSE, 's', ARG_INT, 0.0, 0, NULL},
  { "Masked query output, must be used in conjunction with -D 2 option", 
    NULL, NULL, NULL, TRUE, 'Q', ARG_FILE_OUT, 0.0, 0, NULL},  /* ARG_MASKEDQUERY */
  { "Show full IDs in the output (default - only GIs or accessions)",
        "F", NULL, NULL, FALSE, 'f', ARG_BOOLEAN, 0.0, 0, NULL},   /* ARG_FULLID */
  {"Use lower case filtering of FASTA sequence",                   
        "F", NULL, NULL, TRUE, 'U', ARG_BOOLEAN, 0.0, 0, NULL},    /* ARG_LCASE */
  {"Report the log information at the end of output",
        "F", NULL, NULL, TRUE, 'R', ARG_BOOLEAN, 0.0, 0, NULL},    /* ARG_LOGINFO */
  {"Identity percentage cut-off", 
        "0", NULL, NULL, FALSE, 'p', ARG_FLOAT, 0.0, 0, NULL},     /* ARG_PERC_IDENT */
  { "Location on query sequence",                             
        NULL, NULL, NULL, TRUE, 'L', ARG_STRING, 0.0, 0, NULL},    /* ARG_QUERYLOC */
  { "Multiple Hits window size; default is 0 (i.e. single-hit extensions) or 40 for discontiguous template (negative number overrides this)",  
        "0", NULL, NULL, FALSE, 'A', ARG_INT, 0.0, 0, NULL},      /* ARG_WINDOW  */
  { "X dropoff value for ungapped extension",
    "10", NULL, NULL, FALSE, 'y', ARG_INT, 0.0, 0, NULL},      /* ARG_XDROP_UNGAPPED */
  { "X dropoff value for dynamic programming gapped extension",
    "50", NULL, NULL, FALSE, 'Z', ARG_INT, 0.0, 0, NULL},       /* ARG_XDROP_FINAL */
#ifdef DO_NOT_SUPPRESS_BLAST_OP
  { "Length of a discontiguous word template (contiguous word if 0)",
    "0", NULL, NULL, FALSE, 't', ARG_INT, 0.0, 0, NULL},       /* ARG_TEMPL_LEN */
  {"Make discontiguous megablast generate words for every base of the database (mandatory with the current BLAST engine)",
        "T", NULL, NULL, TRUE, 'g', ARG_BOOLEAN, 0.0, 0, NULL},    /* ARG_EVERYBASE */
  {"Use non-greedy (dynamic programming) extension for affine gap scores",
        "F", NULL, NULL, TRUE, 'n', ARG_BOOLEAN, 0.0, 0, NULL},    /* ARG_DYNAMIC */
  { "Type of a discontiguous word template (0 - coding, 1 - optimal, 2 - two simultaneous",
    "0", NULL, NULL, FALSE, 'N', ARG_INT, 0.0, 0, NULL},       /* ARG_TEMPL_TYPE */
  { "Maximal number of HSPs to save per database sequence (0 = unlimited)",
    "0", NULL, NULL, FALSE, 'H', ARG_INT, 0.0, 0, NULL},        /* ARG_MAXHSP */
#endif
#if MB_ALLOW_NEW
  {"Force use of the legacy BLAST engine",
        "F", NULL, NULL, TRUE, 'V', ARG_BOOLEAN, 0.0, 0, NULL}    /* ARG_FORCE_OLD */
#endif

};

#define MAX_NUM_QUERIES 16383 /* == 1/2 INT2_MAX */

static Int2 Main_old (void)
 
{
   AsnIoPtr aip, xml_aip = NULL;
   BioseqPtr query_bsp, PNTR query_bsp_array;
   BioSourcePtr source;
   BLAST_MatrixPtr matrix;
   BLAST_OptionsBlkPtr options;
   BLAST_KarlinBlkPtr ka_params=NULL, ka_params_gap=NULL;
   BlastPruneSapStructPtr prune;
   Boolean db_is_na, query_is_na, show_gi, believe_query=FALSE;
   Boolean html=FALSE;
   CharPtr params_buffer=NULL;
   Int4 number_of_descriptions, number_of_alignments;
   SeqAlignPtr  seqalign, PNTR seqalign_array;
   SeqAnnotPtr seqannot;
   SeqEntryPtr PNTR sepp;
   TxDfDbInfoPtr dbinfo=NULL, dbinfo_head;
   Uint1 align_type, align_view, out_type;
   Uint4 align_options, print_options;
   ValNodePtr mask_loc, mask_loc_start, next_mask_loc;
   ValNodePtr vnp, other_returns, error_returns;
   
   CharPtr blast_program, blast_database, blast_inputfile, blast_outputfile;
   FILE *infp, *outfp, *mqfp=NULL;
   Int4 index, num_bsps, total_length, total_processed = 0;
   Int2 ctr = 1;
   Char prefix[2];
   SeqLocPtr last_mask, mask_slp;
   Boolean done, hits_found;
   Boolean lcase_masking;
   MBXmlPtr mbxp = NULL;
   Boolean traditional_formatting;

    blast_program = "blastn";
    blast_database = myargs [ARG_DB].strvalue;
    blast_inputfile = myargs [ARG_QUERY].strvalue;
    blast_outputfile = myargs [ARG_OUT].strvalue;
    if (myargs[ARG_HTML].intvalue)
        html = TRUE;

    if ((infp = FileOpen(blast_inputfile, "r")) == NULL) {
       ErrPostEx(SEV_FATAL, 1, 0, "mgblast: Unable to open input file %s\n", blast_inputfile);
       return (1);
    }

    align_view = (Int1) myargs[ARG_FORMAT].intvalue;
    /* Geo mod: 
      -- replaced myargs[ARG_OUTTYPE].intvalue with out_type from now on
    */
    out_type=(Int1) myargs[ARG_OUTTYPE].intvalue;
    if (out_type==MGBLAST_FLTHITS || out_type==MGBLAST_HITGAPS) {
      align_view = 12 + (out_type-MGBLAST_FLTHITS ); 
      out_type=MBLAST_ALIGNMENTS;
      //Attention: 12 MUST be the -m mgblast tab option for MGBLAST_FLTHITS format
      // and MGBLAST_HITGAPS = MGBLAST_FLTHITS+1
       if (align_view>12) { // this is MGBLAST_HITGAPS output
            gap_Info=TRUE;
            if (dbgaps_buf==NULL)
                  dbgaps_buf=(CharPtr) Malloc(dbgaps_bufsize + 1);
            if (qgaps_buf==NULL) 
                qgaps_buf=(CharPtr) Malloc(qgaps_bufsize + 1);
            }
      }

    outfp = NULL;

    traditional_formatting = 
        (out_type == MBLAST_ALIGNMENTS ||
         out_type == MBLAST_DELAYED_TRACEBACK);

    if ((!traditional_formatting ||
            (align_view != 7 && align_view != 10 && align_view != 11)) && 
            blast_outputfile != NULL) {
       if ((outfp = FileOpen(blast_outputfile, "w")) == NULL) {
          ErrPostEx(SEV_FATAL, 1, 0, "blast: Unable to open output file %s\n", blast_outputfile);
          return (1);
       }
    }

    //align_type = BlastGetTypes(blast_program, &query_is_na, &db_is_na);
    align_type = BlastGetTypes(blast_program, &query_is_na, &db_is_na);
    /*
    if (!traditional_formatting)
        believe_query = TRUE;
    else
        believe_query = (Boolean) myargs[ARG_BELIEVEQUERY].intvalue;
    */
    //Geo mod: 
    believe_query=FALSE;
    //If ASN.1 output is requested and believe_query is not set to TRUE,
    //   exit with an error.    
    if (!believe_query && (myargs[ARG_ASNOUT].strvalue ||
                           align_view == 10 || align_view == 11)) {
        ErrPostEx(SEV_FATAL, 1, 0, 
                  "-J option must be TRUE to produce ASN.1 output; before "
                  "changing -J to TRUE please also ensure that all query "
                  "sequence identifiers are unique");
        return -1;
    }
        
    options = BLASTOptionNewEx(blast_program, TRUE, TRUE);
    if (options == NULL)
        return 3;

    options->do_sum_stats = FALSE;
    options->is_neighboring = FALSE;
        options->expect_value  = (Nlm_FloatHi) myargs [ARG_EVALUE].floatvalue;
    number_of_descriptions = myargs[ARG_DESCRIPTIONS].intvalue;    
    number_of_alignments = myargs[ARG_ALIGNMENTS].intvalue;    
    options->hitlist_size = MAX(number_of_descriptions, number_of_alignments);

    if (myargs[ARG_XDROP].intvalue != 0)
           options->gap_x_dropoff = myargs[ARG_XDROP].intvalue;
    if (myargs[ARG_XDROP_UNGAPPED].intvalue != 0)
           options->dropoff_2nd_pass = myargs[ARG_XDROP_UNGAPPED].intvalue;
        if (myargs[ARG_XDROP_FINAL].intvalue != 0)
           options->gap_x_dropoff_final = myargs[ARG_XDROP_FINAL].intvalue;

    if (StringICmp(myargs[ARG_FILTER].strvalue, "T") == 0)
       options->filter_string = StringSave("D");
    else
       options->filter_string = StringSave(myargs[ARG_FILTER].strvalue);
    
    show_gi = (Boolean) myargs[ARG_SHOWGIS].intvalue;
    options->penalty = myargs[ARG_MISMATCH].intvalue;
    options->reward = myargs[ARG_MATCH].intvalue;
        if (myargs[ARG_GAPOPEN].intvalue >= 0)
        options->gap_open = myargs[ARG_GAPOPEN].intvalue;
        if (myargs[ARG_GAPEXT].intvalue >= 0)
        options->gap_extend = myargs[ARG_GAPEXT].intvalue;

    if (options->gap_open == 0 && options->reward % 2 == 0 && 
        options->gap_extend == options->reward / 2 - options->penalty)
       /* This is the default value */
    options->gap_extend = 0;

    options->genetic_code = 1;
    options->db_genetic_code = 1; /* Default; it's not needed here anyway */
    options->number_of_cpus = myargs[ARG_THREADS].intvalue;
    if (myargs[ARG_WORDSIZE].intvalue != 0)
           options->wordsize = myargs[ARG_WORDSIZE].intvalue;
        if (myargs[ARG_MINSCORE].intvalue == 0)
           options->cutoff_s2 = options->wordsize*options->reward;
        else 
           options->cutoff_s2 = myargs[ARG_MINSCORE].intvalue;

        options->db_length = (Int8) myargs[ARG_DBSIZE].floatvalue;
        options->searchsp_eff = (Nlm_FloatHi) myargs[ARG_SEARCHSP].floatvalue;

    options->perform_culling = FALSE;
    /* Kludge */
    options->block_width  = myargs[ARG_MAXPOS].intvalue;

    options->strand_option = myargs[ARG_STRAND].intvalue;
        options->window_size = myargs[ARG_WINDOW].intvalue;
#ifdef DO_NOT_SUPPRESS_BLAST_OP        
        options->mb_template_length = myargs[ARG_TEMPL_LEN].intvalue;
        if (myargs[ARG_TEMPL_LEN].intvalue != 0)
            options->mb_one_base_step = (Boolean) myargs[ARG_EVERYBASE].intvalue;
        options->mb_disc_type = myargs[ARG_TEMPL_TYPE].intvalue;
#endif
        lcase_masking = (Boolean) myargs[ARG_LCASE].intvalue;
        /* Allow dynamic programming gapped extension only with affine 
           gap scores */
        if (options->gap_open != 0 || options->gap_extend != 0)
           options->mb_use_dyn_prog = (Boolean) myargs[ARG_DYNAMIC].intvalue;

        print_options = 0;
        align_options = 0;
        align_options += TXALIGN_COMPRESS;
        align_options += TXALIGN_END_NUM;
        if (show_gi) {
       align_options += TXALIGN_SHOW_GI;
       print_options += TXALIGN_SHOW_GI;
        }
            
        if (align_view) {
       align_options += TXALIGN_MASTER;
       if (align_view == 1 || align_view == 3)
          align_options += TXALIGN_MISMATCH;
       if (align_view == 3 || align_view == 4 || align_view == 6)
          align_options += TXALIGN_FLAT_INS;
       if (align_view == 5 || align_view == 6)
          align_options += TXALIGN_BLUNT_END;
        } else {
       align_options += TXALIGN_MATRIX_VAL;
       align_options += TXALIGN_SHOW_QS;
    }

    if (html) {
       align_options += TXALIGN_HTML;
       print_options += TXALIGN_HTML;
    }

    if (myargs[ARG_GILIST].strvalue)
       options->gifile = StringSave(myargs[ARG_GILIST].strvalue);
   
    if (out_type == MBLAST_ENDPOINTS)
      options->no_traceback = 1;
   else if (out_type == MBLAST_DELAYED_TRACEBACK)
       options->no_traceback = 2;
    else
       options->no_traceback = 0;

    options->megablast_full_deflines = (Boolean) myargs[ARG_FULLID].intvalue;
    options->perc_identity = (FloatLo) myargs[ARG_PERC_IDENT].floatvalue;
    options->hsp_num_max = myargs[ARG_MAXHSP].intvalue;

    if (!believe_query)
           options->megablast_full_deflines = TRUE;
        /*if (options->megablast_full_deflines)
          believe_query = FALSE;*/

    query_bsp_array = (BioseqPtr PNTR) MemNew((MAX_NUM_QUERIES+1)*sizeof(BioseqPtr));
    sepp = (SeqEntryPtr PNTR) MemNew(MAX_NUM_QUERIES*sizeof(SeqEntryPtr));

    StrCpy(prefix, "");

    global_fp = outfp;
        options->output = outfp;

    if (traditional_formatting) {
       if (align_view < 7) {
              if (html) {
                 fprintf(outfp, "<HTML>\n<TITLE>MEGABLAST Search Results</TITLE>\n");
                 fprintf(outfp, "<BODY BGCOLOR=\"#FFFFFF\" LINK=\"#0000FF\" "
                         "VLINK=\"#660099\" ALINK=\"#660099\">\n");
                 fprintf(outfp, "<PRE>\n");
              }
              init_buff_ex(90);
              BlastPrintVersionInfo("mgblast", html, outfp);
              fprintf(outfp, "\n");
              MegaBlastPrintReference(html, 90, outfp);
              fprintf(outfp, "\n");
              
              if(!PrintDbInformation(blast_database, !db_is_na, 70, outfp, html))
                 return 1;
              
              free_buff();
    
#ifdef OS_UNIX
              fprintf(global_fp, "%s", "Searching");
#endif
           }
    }
    
        aip = NULL;
        if (myargs[ARG_ASNOUT].strvalue != NULL) {
           if ((aip = AsnIoOpen (myargs[ARG_ASNOUT].strvalue,"w")) == NULL) {
              ErrPostEx(SEV_FATAL, 1, 0, "blast: Unable to open output file %s\n", myargs[ARG_ASNOUT].strvalue);
              return 1;
           }
        }
        else if (align_view == 10 || align_view == 11)
        {
            const char* mode = (align_view == 10) ? "w" : "wb";
            if ((aip = AsnIoOpen (blast_outputfile, (char*) mode)) == NULL) {
                    ErrPostEx(SEV_FATAL, 1, 0, "blast: Unable to open output file %s\n", blast_outputfile);
                    return 1;
            }
        }


        if (align_view == 7) {
           xml_aip = AsnIoOpen(blast_outputfile, "wx");
        }

        if (myargs[ARG_QUERYLOC].strvalue) {       
            Int4 start, end;
            Megablast_GetLoc(myargs[ARG_QUERYLOC].strvalue, &start, &end);
            options->required_start = start - 1;
            options->required_end = end -1;
        }

    done = FALSE;
    while (!done) {
       num_bsps = 0;
       total_length = 0;
       done = TRUE;
       SeqMgrHoldIndexing(TRUE);
       mask_slp = last_mask = NULL;
   
       while ((sepp[num_bsps]=FastaToSeqEntryForDb(infp, query_is_na, NULL,
                               believe_query, prefix, &ctr, 
                               &mask_slp)) != NULL) {
              if (!lcase_masking) /* Lower case ignored */
                 mask_slp = SeqLocFree(mask_slp);
         if (mask_slp) {
           if (!last_mask)
              options->query_lcase_mask = last_mask = mask_slp;
           else {
              last_mask->next = mask_slp;
              last_mask = last_mask->next;
              }
           mask_slp = NULL;
           }
          query_bsp = NULL;
         SeqEntryExplore(sepp[num_bsps], &query_bsp, FindNuc);
         //debug:
         /*
         char query_buffer[255];
         SeqIdWrite(query_bsp->id, query_buffer, PRINTID_FASTA_LONG, BUFFER_LENGTH);
         fprintf(stderr, "===> query_buf=%s\n", query_buffer);
         */
         if (query_bsp == NULL) {
           ErrPostEx(SEV_FATAL, 1, 0, "Unable to obtain bioseq\n");
           return 2;
          }
          
          source = BioSourceNew();
          source->org = OrgRefNew();
          source->org->orgname = OrgNameNew();
          source->org->orgname->gcode = options->genetic_code;
          ValNodeAddPointer(&(query_bsp->descr), Seq_descr_source, source);
          
          query_bsp_array[num_bsps++] = query_bsp;
          
          total_length += query_bsp->length;
          if (total_length > myargs[ARG_MAXQUERY].intvalue || 
          num_bsps >= MAX_NUM_QUERIES) {
         done = FALSE;
         break;
          }
       }

           if (num_bsps == 0)
               break;

       SeqMgrHoldIndexing(FALSE);
       other_returns = NULL;
       error_returns = NULL;
       
       if (out_type==MBLAST_ENDPOINTS) 
          seqalign_array = BioseqMegaBlastEngine(query_bsp_array, blast_program,
                             blast_database, options,
                             &other_returns, &error_returns,
                             dummy_callback, NULL, NULL, 0, 
                             MegaBlastPrintEndpoints);
       else if (out_type==MBLAST_SEGMENTS) 
          seqalign_array = BioseqMegaBlastEngine(query_bsp_array, blast_program,
                             blast_database, options,
                             &other_returns, &error_returns,
                             dummy_callback, NULL, NULL, 0,
                             MegaBlastPrintSegments);
       else if (out_type==MBLAST_ALIGN_INFO) {
              /* -- Geo mod: do not print header
              PrintTabularOutputHeader(blast_database, 
                                       (num_bsps==1) ? query_bsp_array[0] : NULL,
                                       NULL, "megablast", 0, believe_query,
                                       global_fp);*/
          seqalign_array = BioseqMegaBlastEngine(query_bsp_array, blast_program,
                             blast_database, options,
                             &other_returns, &error_returns,
                             dummy_callback, NULL, NULL, 0,
                                MegaBlastPrintAlignInfo);
       } else if (out_type==MBLAST_ALIGNMENTS) {
          seqalign_array = BioseqMegaBlastEngine(query_bsp_array, blast_program,
                  blast_database, options, &other_returns, 
                                  &error_returns, align_view < 7 ? tick_callback : NULL,
                                  NULL, NULL, 0, NULL);
          }
       
#ifdef OS_UNIX
       fflush(global_fp);
#endif

       if (error_returns) {
             BlastErrorPrint(error_returns);
              for (vnp = error_returns; vnp; vnp = vnp->next) {
                 BlastDestroyErrorMessage((BlastErrorMsgPtr)vnp->data.ptrvalue);
              }
              ValNodeFree(error_returns);
           }
              
              
       if (traditional_formatting) {
          dbinfo = NULL;
          ka_params = NULL;
          ka_params_gap = NULL;
          params_buffer = NULL;
          mask_loc = NULL;
          matrix = NULL;
          for (vnp=other_returns; vnp; vnp = vnp->next) {
           switch (vnp->choice) {
           case TXDBINFO:
              dbinfo = vnp->data.ptrvalue;
              break;
           case TXKABLK_NOGAP:
              ka_params = vnp->data.ptrvalue;
              break;
           case TXKABLK_GAP:
              ka_params_gap = vnp->data.ptrvalue;
              break;
           case TXPARAMETERS:
              params_buffer = vnp->data.ptrvalue;
              break;
           case TXMATRIX:
              matrix = vnp->data.ptrvalue;
              break;
           case SEQLOC_MASKING_NOTSET:
           case SEQLOC_MASKING_PLUS1:
           case SEQLOC_MASKING_PLUS2:
           case SEQLOC_MASKING_PLUS3:
           case SEQLOC_MASKING_MINUS1:
           case SEQLOC_MASKING_MINUS2:
           case SEQLOC_MASKING_MINUS3:
              ValNodeAddPointer(&mask_loc, vnp->choice, vnp->data.ptrvalue);
              break;
           default:
              break;
           }
          }    
          
#ifdef OS_UNIX
        if(align_view < 7) {
                 fprintf(global_fp, "%s\n", " done");
              }
#endif
          
        if (myargs[ARG_MASKEDQUERY].strvalue) {
                 if ((mqfp = FileOpen(myargs[ARG_MASKEDQUERY].strvalue, "w")) == NULL)
                    ErrPostEx(SEV_WARNING, 1, 0, "Unable to open file %s for masked query\n",
                              myargs[ARG_MASKEDQUERY].strvalue);
              }

        hits_found = FALSE;

        mask_loc_start = next_mask_loc = mask_loc;
        mask_loc = NULL;

        if (align_view == 7) {
           mbxp = PSIXmlInit(xml_aip, "megablast", blast_database, 
                             options, query_bsp_array[0], 0);
           }

        if (seqalign_array) { //results returned back for processing
             ReadDBBioseqFetchEnable ("megablast", blast_database, db_is_na, TRUE);
             for (index=0; index<num_bsps; index++) {
                    seqalign = seqalign_array[index];
                    if (next_mask_loc && 
                        SeqIdComp(SeqLocId((SeqLocPtr)next_mask_loc->data.ptrvalue), 
                                  query_bsp_array[index]->id) == SIC_YES) {
                       mask_loc = (SeqLocPtr) 
                       MemDup(next_mask_loc, sizeof(SeqLoc));
                       next_mask_loc = next_mask_loc->next;
                       mask_loc->next = NULL;
                    }
                    if (mqfp) {
                       /* convert mask locations from all sources into
                          a single seqloc */
                       mask_slp = NULL;
                       if (mask_loc) 
                          mask_slp = blastMergeFilterLocs(mask_slp, 
                              (SeqLocPtr)mask_loc->data.ptrvalue,
                              FALSE, 0, 0);
                       PrintMaskedSequence(query_bsp_array[index], mask_slp,
                                           mqfp, 50, lcase_masking);
                       SeqLocSetFree(mask_slp);
                       }
                    if (seqalign==NULL) {
                       mask_loc = MemFree(mask_loc);
                       continue;
                    }
                    hits_found = TRUE;
                    if (align_view < 7) {
                       init_buff_ex(70);
                       AcknowledgeBlastQuery(query_bsp_array[index], 70, outfp, 
                                             believe_query, html);
                       free_buff();
                       }
                    if (align_view == 8 || align_view == 9) {
                       if (align_view == 9)
                          PrintTabularOutputHeader(blast_database, 
                             query_bsp_array[index], NULL, blast_program, 0,
                             believe_query, global_fp);
                       /* debug:
                       char qbuf[512];
                       strcpy(qbuf, BioseqGetTitle(query_bsp_array[index]));
                       fprintf(stderr, "---> Here: query title=%s\n", qbuf);
                       */
                       BlastPrintTabulatedResults(seqalign, 
                           query_bsp_array[index], NULL, number_of_alignments,
                            blast_program, !options->gapped_calculation, 
                            believe_query, 0, 0, 
                            global_fp, (align_view == 9));
                            

                       ObjMgrFreeCache(0);

                       SeqAlignSetFree(seqalign);
                       mask_loc = MemFree(mask_loc);
                       continue;
                    } 
                       //Geo mod:   
                   else if (align_view>=12)  {
                        MGBlastPrintTab(seqalign, 
                            query_bsp_array[index], number_of_alignments,
                            !options->gapped_calculation, 
                            global_fp);
                        ObjMgrFreeCache(0);

                        SeqAlignSetFree(seqalign);
                        mask_loc = MemFree(mask_loc);
                        continue;
                        }
                    else if(align_view == 7) {
                       IterationPtr iterp;

                       iterp = BXMLBuildOneQueryIteration(seqalign, 
                                  NULL, FALSE, 
                                  !options->gapped_calculation, index, 
                                  NULL, query_bsp_array[index], mask_loc);
                       IterationAsnWrite(iterp, mbxp->aip, mbxp->atp);
                       AsnIoFlush(mbxp->aip);
                       IterationFree(iterp);
                       SeqAlignSetFree(seqalign);
                       mask_loc = MemFree(mask_loc);
                       continue;
                    }
                    seqannot = SeqAnnotNew();
                    seqannot->type = 2;
                    AddAlignInfoToSeqAnnot(seqannot, align_type);
                    seqannot->data = seqalign;
                    if (aip) {
                       SeqAnnotAsnWrite((SeqAnnotPtr) seqannot, aip, NULL);
                       AsnIoReset(aip);
                    }
                    if (outfp) { /* Uncacheing causes problems with ordinal nos. vs. gi's. */
                       prune = BlastPruneHitsFromSeqAlign(seqalign, number_of_descriptions, NULL);
                       ObjMgrSetHold();
                       init_buff_ex(85);
                       PrintDefLinesFromSeqAlign(prune->sap, 80,
                                                 outfp, print_options, FIRST_PASS, NULL);
                       free_buff();
                       
                       prune = BlastPruneHitsFromSeqAlign(seqalign, number_of_alignments, prune);
                       seqannot->data = prune->sap;
                       if (align_view != 0)
                          ShowTextAlignFromAnnot(seqannot, 60, outfp, NULL,
                                                 NULL, align_options, NULL, 
                                                 mask_loc, NULL);
                       else
                          ShowTextAlignFromAnnot(seqannot, 60, outfp, NULL, NULL, align_options, NULL, mask_loc, FormatScoreFunc);
                       seqannot->data = seqalign;
                       prune = BlastPruneSapStructDestruct(prune);
                       ObjMgrClearHold();
                       ObjMgrFreeCache(0);
                    }
                    seqannot = SeqAnnotFree(seqannot);
                    mask_loc = MemFree(mask_loc);
                 } /* End loop on seqaligns for different queries */
                 ReadDBBioseqFetchDisable();
              } 

              if (mbxp != NULL) {
                 MBXmlClose(mbxp, other_returns, !options->gapped_calculation);
              }

              if (mqfp)
                 FileClose(mqfp);

              if (!hits_found && align_view < 7)
                 fprintf(outfp, "\n\n ***** No hits found ******\n\n");

              matrix = BLAST_MatrixDestruct(matrix);
          
              if(html) 
                 fprintf(outfp, "<PRE>\n");
              init_buff_ex(85);
              dbinfo_head = dbinfo;
              if(align_view < 7) {
                 while (dbinfo) {
                    PrintDbReport(dbinfo, 70, outfp);
                    dbinfo = dbinfo->next;
                 }
              }
              dbinfo_head = TxDfDbInfoDestruct(dbinfo_head);
              
              if (ka_params) {
                 if(align_view < 7)
                    PrintKAParameters(ka_params->Lambda, ka_params->K, ka_params->H, 70, outfp, FALSE);
                 MemFree(ka_params);
              }
              if (ka_params_gap) {
                 if(align_view < 7)
                    PrintKAParameters(ka_params_gap->Lambda, ka_params_gap->K, ka_params_gap->H, 70, outfp, TRUE);
                 MemFree(ka_params_gap);
              }
              if(align_view < 7)
                 PrintTildeSepLines(params_buffer, 70, outfp);
              MemFree(params_buffer);
              free_buff();
              mask_loc = mask_loc_start;
              while (mask_loc) {
                 SeqLocSetFree(mask_loc->data.ptrvalue);
                 mask_loc = mask_loc->next;
              }
              ValNodeFree(mask_loc_start);
       } else { //not traditional formatting
          /* Just destruct all other_returns parts */
          for (vnp=other_returns; vnp; vnp = vnp->next) {
         switch (vnp->choice) {
         case TXDBINFO:
            TxDfDbInfoDestruct(vnp->data.ptrvalue);
            break;
         case TXKABLK_NOGAP:
         case TXKABLK_GAP:
         case TXPARAMETERS:
            MemFree(vnp->data.ptrvalue);
            break;
         case TXMATRIX:
            BLAST_MatrixDestruct(vnp->data.ptrvalue);
            break;
         case SEQLOC_MASKING_NOTSET:
         case SEQLOC_MASKING_PLUS1:
         case SEQLOC_MASKING_PLUS2:
         case SEQLOC_MASKING_PLUS3:
         case SEQLOC_MASKING_MINUS1:
         case SEQLOC_MASKING_MINUS2:
         case SEQLOC_MASKING_MINUS3:
                    mask_loc = vnp->data.ptrvalue;
                    SeqLocSetFree(mask_loc);
         default:
            break;
         }
          }
       }
       other_returns = ValNodeFree(other_returns);
       MemFree(seqalign_array);
           options->query_lcase_mask = 
              SeqLocSetFree(options->query_lcase_mask);

       /* Freeing SeqEntries can be very expensive, do this only if 
          this is not the last iteration of search */
       if (!done) { 
          for (index=0; index<num_bsps; index++) {
         sepp[index] = SeqEntryFree(sepp[index]);
         query_bsp_array[index] = NULL;
          }       
           }
           total_processed += num_bsps;
    } /* End of loop on complete searches */
        
        aip = AsnIoClose(aip);

        /*if (align_view == 7)
          xml_aip = AsnIoClose(xml_aip);*/

        if (align_view < 7 && html) 
           fprintf(outfp, "</PRE>\n</BODY>\n</HTML>\n");
        if (align_view < 7 && myargs[ARG_LOGINFO].intvalue)
           fprintf(outfp, "Mega BLAST run finished, processed %d queries\n",
                   total_processed);
    MemFree(query_bsp_array);
    MemFree(sepp);
    MemFree(qgaps_buf);
    MemFree(dbgaps_buf);
    options = BLASTOptionDelete(options);
    FileClose(infp);
        FileClose(outfp);
    
    return 0;
}


Int2 Nlm_Main(void) {
    Boolean use_new_engine=FALSE;
    char buf[256] = { '\0' };

    StringCpy(buf, "mgblast ");
    StringNCat(buf, BlastGetVersionNumber(), sizeof(buf)-StringLen(buf)-1);
    if (! GetArgs (buf, NUMARG, myargs))
       return (1);

    UseLocalAsnloadDataAndErrMsg ();

    if (! SeqEntryLoad())
        return 1;

    ErrSetMessageLevel(SEV_WARNING);
    /*
    if (myargs[ARG_FORCE_OLD].intvalue == 0 &&
                  myargs[ARG_OUTTYPE].intvalue > 1 &&
                      myargs[ARG_GILIST].strvalue == NULL)
          use_new_engine = readdb_use_new_blast(myargs[ARG_DB].strvalue);
    if (myargs[ARG_MASKEDQUERY].strvalue)
        use_new_engine = FALSE;
    */
    
    /*if (use_new_engine)
        return Main_new();
    else */
        return Main_old();
}
