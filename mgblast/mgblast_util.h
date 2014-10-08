#ifndef MGBLAST_UTIL_H
#define MGBLAST_UTIL_H

#include <objsset.h>
#include <blast.h>
#include <sqnutils.h>

extern int max_overhang;
extern int max_overhang;
extern int min_overlap;
extern int slice_clustering;
extern int gap_Info;
extern int db_skipto;
extern int appendName; /* for -D4, appends query or subject defline */
extern CharPtr dbgaps_buf;
extern CharPtr qgaps_buf;
extern Int4 dbgaps_buf_used;
extern Int4 dbgaps_bufsize;
extern Int4 qgaps_buf_used;
extern Int4 qgaps_bufsize;

/*-- how many query slicing iterations? */
extern int qread_base;
/* int max_num_queries; */

void MGBlastPrintTab(SeqAlignPtr seqalign, BioseqPtr query_bsp,
        Int4 num_alignments, Boolean is_ungapped, FILE *fp);

#endif
