#include <objsset.h>
#include <blast.h>
#include <txalign.h>
#include <mblast.h>
#include <xmlblast.h>
#include <sqnutils.h>
#include <blfmtutl.h>
#include "mgblast_util.h"

int max_overhang=INT_MAX;
int min_overlap=0;
int slice_clustering=FALSE;
int gap_Info=FALSE;
int db_skipto=0;
int appendName=0; /* for -D4, appends query or subject defline */
CharPtr dbgaps_buf=NULL;
CharPtr qgaps_buf=NULL;
Int4 dbgaps_buf_used=0;
Int4 dbgaps_bufsize=1024;
Int4 qgaps_buf_used=0;
Int4 qgaps_bufsize=1024;
/*-- how many query slicing iterations? */
int qread_base = 0;
/* extern int max_num_queries = 4000; */

static Char tmp_buf[1024]; /* temporary text buffer for preparing gap info */

/*
 MGBlastPrintTab() based on BlastPrintTabularResults() in blfmtutl.c
*/

#define INTSWAP(a,b) a ^= b ^= a ^= b


/* utility function for gap buffers: */

gaplistAdd(CharPtr* buf, Int4* buf_used, Int4* buf_cap, CharPtr addstr) {
  int addlen=strlen(addstr);
  if (*buf_used+addlen+1>=*buf_cap) {
     *buf_cap+=64;
     *buf=(CharPtr) Realloc(*buf, *buf_cap);
     }
  if ((*buf)[0]!='\0') {
     strcat(*buf, ",");
     *buf_used++;
     }
  strcat(*buf, addstr);
  *buf_used+=addlen;
}

/* BlastBioseqGetNumIdentical() 
   lifted from blfmtutl.c, needed by MGBlastPrintTab() */

static Int4
BlastBioseqGetNumIdentical(BioseqPtr q_bsp, BioseqPtr s_bsp, Int4 q_start,
                     Int4 s_start, Int4 length,
                     Uint1 q_strand, Uint1 s_strand)
{
   SeqLocPtr q_slp, s_slp;
   SeqPortPtr q_spp, s_spp;
   Int4 i, ident = 0;
   Uint1 q_res, s_res;

   if (!q_bsp || !s_bsp)
      return 0;

   q_slp = SeqLocIntNew(q_start, q_start+length-1, q_strand, q_bsp->id);
   s_slp = SeqLocIntNew(s_start, s_start+length-1, s_strand, s_bsp->id);
   if (ISA_na(q_bsp->mol))
      q_spp = SeqPortNewByLoc(q_slp, Seq_code_ncbi4na);
   else
      q_spp = SeqPortNewByLoc(q_slp, Seq_code_ncbistdaa);
   if (ISA_na(s_bsp->mol))
      s_spp = SeqPortNewByLoc(s_slp, Seq_code_ncbi4na);
   else
      s_spp = SeqPortNewByLoc(s_slp, Seq_code_ncbistdaa);

   for (i=0; i<length; i++) {
      while ((q_res = SeqPortGetResidue(q_spp)) != SEQPORT_EOF &&
             !IS_residue(q_res));
      while ((s_res = SeqPortGetResidue(s_spp)) != SEQPORT_EOF &&
             !IS_residue(s_res));
      if (q_res == SEQPORT_EOF || s_res == SEQPORT_EOF)
         break;
      else if (q_res == s_res)
         ident++;
   }

   SeqLocFree(q_slp);
   SeqLocFree(s_slp);
   SeqPortFree(q_spp);
   SeqPortFree(s_spp);

   return ident;
}

/* modified from ScoreAndEvalueToBuffers() in txalign.c */
void formatScores (FloatHi bit_score, FloatHi evalue,
                        CharPtr bit_score_buf, CharPtr evalue_buf)
{
 unsigned int int_bitscore=(unsigned int) bit_score;
#ifdef OS_MAC
   if (evalue < 1.0e-280) {
      sprintf(evalue_buf, "0.0");
   } else if (evalue < 1.0e-99) {
      sprintf(evalue_buf, "%2.0Le", evalue);
      /* if (format_options & TX_KNOCK_OFF_ALLOWED)
         (*evalue_buf)++;  Knock off digit. */
   } else if (evalue < 0.0009) {
      sprintf(evalue_buf, "%3.0Le", evalue);
   } else if (evalue < 0.1) {
      sprintf(evalue_buf, "%4.3Lf", evalue);
   } else if (evalue < 1.0) {
      sprintf(evalue_buf, "%3.2Lf", evalue);
   } else if (evalue < 10.0) {
      sprintf(evalue_buf, "%2.1Lf", evalue);
   } else {
      sprintf(evalue_buf, "%5.0Lf", evalue);
   }
   /* if (bit_score > 9999)
      sprintf(bit_score_buf, "%4.3Le", bit_score);
   else if (bit_score > 99.9)
      sprintf(bit_score_buf, "%4.0ld", (long)bit_score);
   else // %4.1Lf is bad on 68K Mac, so cast to long
      sprintf(bit_score_buf, "%4.0ld", (long)bit_score);
      */
#else
   if (evalue < 1.0e-280) {
      sprintf(evalue_buf, "0.0");
   } else if (evalue < 1.0e-99) {
      sprintf(evalue_buf, "%2.0le", evalue);
      /*
      if (format_options & TX_KNOCK_OFF_ALLOWED)
         (*evalue_buf)++;  Knock off digit. */
   } else if (evalue < 0.0009) {
      sprintf(evalue_buf, "%3.0le", evalue);
   } else if (evalue < 0.1) {
      sprintf(evalue_buf, "%4.3lf", evalue);
   } else if (evalue < 1.0) {
      sprintf(evalue_buf, "%3.2lf", evalue);
   } else if (evalue < 10.0) {
      sprintf(evalue_buf, "%2.1lf", evalue);
   } else {
      sprintf(evalue_buf, "%5.0lf", evalue);
   }
   /*
   if (bit_score > 9999)
      sprintf(bit_score_buf, "%4.3le", bit_score);
   else if (bit_score > 99.9)
      sprintf(bit_score_buf, "%4.0ld", (long)bit_score);
   else if (format_options & TX_INTEGER_BIT_SCORE)
      sprintf(bit_score_buf, "%4.0lf", bit_score);
   else
      sprintf(bit_score_buf, "%4.1lf", bit_score);
  */
#endif
  sprintf(bit_score_buf, "%d", int_bitscore);
}


/* 
{Geo:} this should print this tab-delimited output:

 q_name, q_len, q_n5, q_n3, hit_name, hit_len, hit_n5, hit_n3, pid, score, p-value, strand 

 for MBLAST_HITGAPS option we have the same output but with 2 extra fields added:
 
 <qgaps> <hitgaps>
 
 ..where each of them has the gap coordinate and gap length, like this:
 
 gap1pos_gap1len[+gap1len],gap2pos[+gap2len],...
 
*/

void MGBlastPrintTab(SeqAlignPtr seqalign, BioseqPtr query_bsp,
        Int4 num_alignments, Boolean is_ungapped, 
        FILE *fp)
{
   // -- un-needed parameters:
    SeqLocPtr query_slp=NULL;
    CharPtr blast_program=NULL; //"megablast"
    Boolean is_ooframe=FALSE;
    Boolean believe_query=FALSE;
    Int4 q_shift=0;
    Int4 s_shift=0;
    Boolean print_query_info=FALSE;
    int *num_formatted=NULL;
   // ---- they should not be used below:
   SeqAlignPtr sap, sap_tmp = NULL;
   FloatHi perc_ident, bit_score, evalue;
   Int4 numseg, num_gap_opens, num_mismatches, num_ident, score;
   Int4 number, align_length, index, i, i2, j;
   Int4 q_start, q_end, s_start, s_end, q_len, s_len;
   Char bit_score_buff[10];
   Char eval_buff[16];
   Boolean is_translated;
   SeqIdPtr query_id, old_query_id = NULL, subject_id, old_subject_id = NULL;
   BioseqPtr subject_bsp=NULL;
   Char query_buffer[BUFFER_LENGTH+1], subject_buffer[BUFFER_LENGTH+1];
   char aln_strand='+';
   DenseSegPtr dsp;
   StdSegPtr ssp = NULL;
   DenseDiagPtr ddp = NULL;
   AlignSumPtr asp = NULL;
   CharPtr defline, title;
   SeqLocPtr slp;
   Int4 alignments_count;
   Int4 objmgr_count = 0;
   /* Geo add: */
   Int4 dbgaplen, qgaplen;
   Int4 qseg_start, qseg_end, dbseg_start, dbseg_end, qseg_prev_end, dbseg_prev_end;


   /*is_translated = (StringCmp(blast_program, "blastn") &&
                    StringCmp(blast_program, "blastp"));*/
   is_translated = FALSE;
   if (is_translated) {
      asp = MemNew(sizeof(AlignSum));
      asp->matrix = load_default_matrix();
      asp->is_aa = TRUE;
      asp->ooframe = is_ooframe;
   }

   if (is_ungapped)
      sap_tmp = SeqAlignNew();

   slp = query_slp;
   q_len=0;
   if (query_bsp) {
      query_id = query_bsp->id;
      q_len=BioseqGetLen(query_bsp);
      }

   /* Evalue buffer is dynamically allocated to avoid compiler warnings 
      in calls to ScoreAndEvalueToBuffers. */
   //eval_buff = Malloc(10);
   s_len=0;
   for (sap = seqalign; sap; sap = sap->next) {
      if (query_slp)
         query_id = TxGetQueryIdFromSeqAlign(sap);
      if (SeqIdComp(query_id, old_query_id) != SIC_YES) {
         if (old_query_id && num_formatted)
            (*num_formatted)++;
         alignments_count = num_alignments;
         /* New query: find the corresponding SeqLoc */
         while (slp && SeqIdComp(query_id, SeqLocId(slp)) != SIC_YES)
            slp = slp->next;
         if (slp != NULL) {
            query_id = old_query_id = SeqLocId(slp);
            /* Print new query information */
            /*
            if (print_query_info)
               PrintTabularOutputHeader(NULL, NULL, slp, NULL, 0, 
                                        believe_query, fp); */
         } else if (query_bsp)
            old_query_id = query_bsp->id;
         defline = (CharPtr) Malloc(BUFFER_LENGTH+1);
         SeqIdWrite(query_id, defline, PRINTID_FASTA_LONG, BUFFER_LENGTH);
         if (StringNCmp(defline, "lcl|", 4))
            StringCpy(query_buffer, defline);
         else if (!believe_query) {
            if (slp) {
               BioseqUnlock(query_bsp);
               query_bsp = BioseqLockById(query_id);
            }
            if ((title = StringSave(BioseqGetTitle(query_bsp))) != NULL) {
               defline = MemFree(defline);
               defline = StringTokMT(title, " ", &title);
               StringNCpy_0(query_buffer, defline, BUFFER_LENGTH);
               defline = MemFree(defline);
            } else
               StringCpy(query_buffer, defline+4);
            defline = MemFree(defline);
         } else
            StringCpy(query_buffer, defline+4);
      } else
         query_id = old_query_id;      

      subject_id = TxGetSubjectIdFromSeqAlign(sap);
      if (SeqIdComp(subject_id, old_subject_id) != SIC_YES) {
         /* New subject sequence has been found in the seqalign list */
         if (--alignments_count < 0)
            continue;
         BioseqUnlock(subject_bsp);

         /* object manager cache is limited in size */
         if (++objmgr_count > 8000) {
            objmgr_count = 0;
            ObjMgrFreeCache(OBJ_MAX);
         }

         subject_bsp = BioseqLockById(subject_id);
         if (!subject_bsp || !subject_bsp->id)
            continue;
         s_len=BioseqGetLen(subject_bsp);
         if (subject_bsp->id->choice != SEQID_GENERAL ||
             StringCmp(((DbtagPtr)subject_id->data.ptrvalue)->db, "BL_ORD_ID")) {
            SeqIdPtr use_this_gi_id = GetUseThisGi(sap); 
            defline = (CharPtr) Malloc(BUFFER_LENGTH+1);
            if (use_this_gi_id) {
                BlastDefLinePtr bdlp, actual_bdlp;
                bdlp=FDGetDeflineAsnFromBioseq(subject_bsp);
                actual_bdlp=getBlastDefLineForSeqId(bdlp, use_this_gi_id);
                
                SeqIdWrite(actual_bdlp->seqid, defline, PRINTID_FASTA_LONG, BUFFER_LENGTH);
                BlastDefLineSetFree(bdlp);
            } else {
                SeqIdWrite(subject_bsp->id, defline, PRINTID_FASTA_LONG, BUFFER_LENGTH);
            }
            if (StringNCmp(defline, "lcl|", 4))
               StringCpy(subject_buffer, defline);
            else
               StringCpy(subject_buffer, defline+4);
         } else {
            defline = StringSave(BioseqGetTitle(subject_bsp));
            defline = StringTokMT(defline, " \t", &title);
            StringCpy(subject_buffer, defline);
         }
         defline = MemFree(defline);
      }
      
      perc_ident = 0;
      align_length = 0;
      num_gap_opens = 0;
      num_mismatches = 0;

      GetScoreAndEvalue(sap, &score, &bit_score, &evalue, &number);

      //ScoreAndEvalueToBuffers()
       formatScores(bit_score, evalue, bit_score_buff, eval_buff);

      /* Loop on segments within this seqalign (in ungapped case) */
      while (TRUE) {
         if (sap->segtype == SAS_DENSEG) { /* Geo: this is the typical megablast case */
            Boolean get_num_ident = TRUE;
            dsp = (DenseSegPtr) sap->segs;
            numseg = dsp->numseg;
            /* Query Bioseq is needed for calculating number of identities.
               NB: even if number of identities is already filled in the 
               seqalign score list, that is not enough here, because we need to
               know number of identities in each segment in order to calculate
               number of mismatches correctly. */
            if (!query_bsp) {
               query_bsp = BioseqLockById(query_id);
            }
            qseg_prev_end=-1;
            dbseg_prev_end=-1;

            aln_strand= (dsp->strands[0] == dsp->strands[1]) ? '+' : '-';
            if (gap_Info) {
              dbgaps_buf[0]='\0';dbgaps_buf_used=1;
              qgaps_buf[0]='\0';qgaps_buf_used=1;
              for (i=0; i<numseg; i++) {
                 align_length += dsp->lens[i];
                 i2=i<<1;
                 if (dsp->starts[i2] != -1 && dsp->starts[i2+1] != -1) {
                   if (get_num_ident) {
                      num_ident = BlastBioseqGetNumIdentical(query_bsp, subject_bsp,
                                     dsp->starts[i2], dsp->starts[i2+1],
                                     dsp->lens[i], dsp->strands[i2],
                                     dsp->strands[i2+1]);
                      perc_ident += num_ident;
                      num_mismatches += dsp->lens[i] - num_ident;
                      }
                   /* Geo: keep track of gap information here */
                   qseg_start = dsp->starts[i2] + 1;
                   qseg_end = qseg_start + dsp->lens[i] - 1;

                   if (aln_strand == '+') {
                          /* qseg_start = dsp->starts[i2] + 1;
                          qseg_end = qseg_start + dsp->lens[i] - 1; */
                          dbseg_start = dsp->starts[i2+1] + 1;
                          dbseg_end = dbseg_start + dsp->lens[i] - 1;
                          } else {
                          /* qseg_start = q_len - dsp->starts[i2];
                          qseg_end = qseg_start - dsp->lens[i] + 1; */
                          dbseg_start = s_len - dsp->starts[i2+1];
                          dbseg_end = dbseg_start - dsp->lens[i] + 1;
                          }
                   /* if (qseg_prev_end<0) {
                          qseg_prev_end=qseg_end;
                          //dbseg_prev_end=dsp->starts[i2+1]+dsp->lens[i];
                          dbseg_prev_end=dbseg_end;
                          continue; //end of gap chain?
                          } */
                   if (dbseg_prev_end<0) {
                        dbseg_prev_end=dbseg_end;
                        //qseg_prev_end=dsp->starts[i2]+dsp->lens[i];
                        qseg_prev_end=qseg_end;
                        continue; //end of gap chain?
                        }
                   dbgaplen=abs(qseg_start-qseg_prev_end)-1;
                     //qgaplen=dsp->starts[i2+1]-dbseg_prev_end;
                   qgaplen=abs(dbseg_start-dbseg_prev_end)-1;
                   if (qgaplen>0) {/* gap in query */
                         /* fprintf(stderr, "gaplen=%d in query\n", qgaplen);*/
                         if (aln_strand=='-') {
                            dbseg_prev_end--; //? adjustment needed
                            }
                         if (qgaplen>1) sprintf(tmp_buf, "%d+%d", qseg_prev_end+1, qgaplen);
                                   else sprintf(tmp_buf, "%d",    qseg_prev_end+1);
                         gaplistAdd(&qgaps_buf, &qgaps_buf_used, &qgaps_bufsize, tmp_buf);
                         } /* gap in query */
                   tmp_buf[0]='\0';
                   if (dbgaplen>0) {/* gap in db seq */
                       if (dbgaplen>1) sprintf(tmp_buf, "%d+%d", dbseg_prev_end+1, dbgaplen);
                                  else sprintf(tmp_buf, "%d", dbseg_prev_end+1);
                       gaplistAdd(&dbgaps_buf, &dbgaps_buf_used, &dbgaps_bufsize, tmp_buf);
                       } /* gap in db seq */
                    qseg_prev_end=qseg_end;
                    dbseg_prev_end=dbseg_end;
                  /* gap info requested */
                  } else {
                    if (dsp->starts[i2] < 0) {
                       //TODO: WORKS - check it again!
                       int qcoord=(aln_strand=='+') ? dsp->starts[2*(i+1)]+1 : dsp->starts[2*(i-1)]+1;
                       fprintf(stderr, "gap in Qry  at %d+%d)\n",qcoord, dsp->lens[i]);
                       }
                    else if (dsp->starts[i2+1] < 0) {
                      fprintf(stderr, "gap in Subj at %d (len %d)\n",dsp->starts[2*(i+1)+1]+1, dsp->lens[i]);
                       }
                    num_gap_opens++;
                    }
               } /* for each segment */
              }
            else { /* no gap tracking */
            for (i=0; i<numseg; i++) {
               align_length += dsp->lens[i];
               i2=i<<1;
               if (dsp->starts[i2] != -1 && dsp->starts[i2+1] != -1) {
                 if (get_num_ident) {
                    num_ident = BlastBioseqGetNumIdentical(query_bsp, subject_bsp,
                                   dsp->starts[i2], dsp->starts[i2+1],
                                   dsp->lens[i], dsp->strands[i2],
                                   dsp->strands[i2+1]);
                    perc_ident += num_ident;
                    num_mismatches += dsp->lens[i] - num_ident;
                    }
               } else {
                  num_gap_opens++;
               }
             }
            } // no gap info needed
            perc_ident = perc_ident / align_length * 100;
            /* compute half the sequence offsets (account for
               leading gaps in the alignment) */
            if (dsp->starts[0] == -1) {
                i = 1; j = 0;
            }
            else if (dsp->starts[1] == -1) {
                i = 0; j = 1;
            }
            else {
                i = j = 0;
            }
            if (dsp->strands[0] != dsp->strands[1]) {
               q_end = dsp->starts[2*i] + dsp->lens[i];
               s_end = dsp->starts[2*j+1] + 1;
            } else {
               q_start = dsp->starts[2*i] + 1;
               s_start = dsp->starts[2*j+1] + 1;
            }

            /* compute half the sequence offsets (account for
               trailing gaps in the alignment) */
            if (dsp->starts[2*numseg-2] == -1) {
                i = numseg-1; j = numseg;
            }
            else if (dsp->starts[2*numseg-1] == -1) {
                i = numseg; j = numseg-1;
            }
            else {
                i = j = numseg;
            }
            if (dsp->strands[0] != dsp->strands[1]) {
               q_start = dsp->starts[2*i-2] + 1;
               s_start = dsp->starts[2*j-1] + dsp->lens[j-1];
            } else {
               q_end = dsp->starts[2*i-2] + dsp->lens[i-1];
               s_end = dsp->starts[2*j-1] + dsp->lens[j-1];
            }

         } else if (sap->segtype == SAS_STD) {

            if (!ssp)
               ssp = (StdSegPtr) sap->segs;
            
            if (is_ungapped) {
               sap_tmp->segtype = SAS_STD;
               sap_tmp->segs = ssp;
               GetScoreAndEvalue(sap_tmp, &score, &bit_score, &evalue, &number);
               formatScores(bit_score, evalue, bit_score_buff, eval_buff);
               find_score_in_align(sap_tmp, 1, asp);
            } else
               find_score_in_align(sap, 1, asp);
            
            if (asp->m_frame < 0)
               q_start = SeqLocStop(ssp->loc) + 1;
            else
               q_start = SeqLocStart(ssp->loc) + 1;
            
            if (asp->t_frame < 0)
               s_start = SeqLocStop(ssp->loc->next) + 1;
            else
               s_start = SeqLocStart(ssp->loc->next) + 1;
            
            if (!is_ungapped) {
               for (index=1; ssp->next; index++)
                  ssp = ssp->next;
               num_gap_opens = index / 2;
            } else 
               num_gap_opens = 0;

            if (asp->m_frame < 0)
               q_end = SeqLocStart(ssp->loc) + 1;
            else
               q_end = SeqLocStop(ssp->loc) + 1;
            
            if (asp->t_frame < 0)
               s_end = SeqLocStart(ssp->loc->next) + 1;
            else
               s_end = SeqLocStop(ssp->loc->next) + 1;
            
            align_length = asp->totlen;
            num_mismatches = asp->totlen - asp->gaps - asp->identical;
            perc_ident = ((FloatHi) 100*asp->identical)/ (asp->totlen);
         } else if (sap->segtype == SAS_DENDIAG) {
           if (!ddp)
               ddp = (DenseDiagPtr) sap->segs;
            sap_tmp->segtype = SAS_DENDIAG;
            sap_tmp->segs = ddp;
            GetScoreAndEvalue(sap_tmp, &score, &bit_score, &evalue, &number);
            formatScores(bit_score, evalue, bit_score_buff, eval_buff);
            align_length = ddp->len;
            /*always show plus strand for query*/
            if (ddp->strands[0] == Seq_strand_minus &&
                ddp->strands[1] == Seq_strand_plus) { 
                ddp->strands[0] = Seq_strand_plus;
                ddp->strands[1] = Seq_strand_minus;   
            }
            aln_strand=(ddp->strands[0] == ddp->strands[1]) ? '+' : '-';
            if (ddp->strands[0] == Seq_strand_minus) {
               q_start = ddp->starts[0] + align_length;
               q_end = ddp->starts[0] + 1;
            } else {
               q_start = ddp->starts[0] + 1;
               q_end = ddp->starts[0] + align_length;
            }

            if (ddp->strands[1] == Seq_strand_minus) {
               s_start = ddp->starts[1] + align_length;
               s_end = ddp->starts[1] + 1;
            } else {
               s_start = ddp->starts[1] + 1;
               s_end = ddp->starts[1] + align_length;
            }
            num_gap_opens = 0;
            /* Query Bioseq is needed for calculating number of identities.
               NB: even if number of identities is already filled in the 
               seqalign score list, that is not enough here, because we need to
               know number of identities in each segment in order to calculate
               number of mismatches correctly. */
            if (!query_bsp) {
               query_bsp = BioseqLockById(query_id);
            }

            num_ident = BlastBioseqGetNumIdentical(query_bsp, subject_bsp, 
                           ddp->starts[0], ddp->starts[1], align_length, 
                           ddp->strands[0], ddp->strands[1]);
            num_mismatches = align_length - num_ident;
            perc_ident = ((FloatHi)num_ident) / align_length * 100;
         }
         if (!is_translated) {
            /* Adjust coordinates if query and/or subject is a subsequence */
            q_start += q_shift;
            q_end += q_shift;
            s_start += s_shift;
            s_end += s_shift;
         }
         
         if (perc_ident >= 99.995 && perc_ident < 100.00)
            perc_ident = 99.99;
         
         /* fprintf(fp,
                 "%s\t%s\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n",
                 query_buffer, subject_buffer, perc_ident, align_length, 
                 num_mismatches, num_gap_opens, q_start, 
                 q_end, s_start, s_end, eval_buff, bit_score_buff);
                 */
         if (s_start>s_end) { INTSWAP(s_start,s_end); }
         if ((q_start>q_end && aln_strand=='+') || (q_start<q_end && aln_strand=='-')) {
             INTSWAP(q_start,q_end);
             }
         if (gap_Info)
           fprintf(fp,"%s\t%d\t%d\t%d\t%s\t%d\t%d\t\%d\t%.2f\t%s\t%s\t%c\t%s\t%s\n",
                          query_buffer, q_len, q_start, q_end, subject_buffer, s_len, s_start, s_end,
                          perc_ident, bit_score_buff, eval_buff, aln_strand, qgaps_buf, dbgaps_buf);
         else
           fprintf(fp,"%s\t%d\t%d\t%d\t%s\t%d\t%d\t\%d\t%.2f\t%s\t%s\t%c\n",
                          query_buffer, q_len, q_start, q_end, subject_buffer, s_len, s_start, s_end,
                          perc_ident, bit_score_buff, eval_buff, aln_strand);

         old_subject_id = subject_id;
         if (sap->segtype == SAS_DENSEG)
            break;
         else if (sap->segtype == SAS_DENDIAG) {
            if ((ddp = ddp->next) == NULL)
               break;
         } else if (sap->segtype == SAS_STD) {
            if ((ssp = ssp->next) == NULL)
               break;
         }
      }
   }

   if (is_ungapped)
      sap_tmp = MemFree(sap_tmp);

   if (is_translated) {
      free_default_matrix(asp->matrix);
      MemFree(asp);
   }

   BioseqUnlock(subject_bsp);
   if (query_slp)
      BioseqUnlock(query_bsp);
} /* MGBlastPrintTab -- based on BlastPrintTabularResults from blfmtutl.c */

