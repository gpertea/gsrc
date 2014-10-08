28a29
    Char buf[256] = { '\0' };
32d32
<         const char *dummystr;
33a34,35
    BlastOutputPtr boutp = NULL;
    MBXmlPtr mbxp = NULL;
35c37,39
<         if (! GetArgs ("megablast", NUMARG, myargs))
---
     StringCpy(buf, "megablast ");
     StringNCat(buf, BlastGetVersionNumber(), sizeof(buf)-StringLen(buf)-1);
     if (! GetArgs (buf, NUMARG, myargs))
53c57
< 	   ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open input file %s\n", blast_inputfile);
---
 	   ErrPostEx(SEV_FATAL, 1, 0, "blast: Unable to open input file %s\n", blast_inputfile);
56a61
 	align_view = (Int1) myargs[3].intvalue;
58c63
< 	if (blast_outputfile != NULL) {
---
 	if (align_view != 7 && align_view != 10 && align_view != 11 && blast_outputfile != NULL) {
60c65
< 	      ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open output file %s\n", blast_outputfile);
---
 	      ErrPostEx(SEV_FATAL, 1, 0, "blast: Unable to open output file %s\n", blast_outputfile);
65,66d69
< 	align_view = (Int1) myargs[3].intvalue;
< 
82,83c85,86
< 	if (believe_query == FALSE && myargs[14].strvalue) 
< 	   ErrPostEx(SEV_FATAL, 0, 0, "-J option must be TRUE to produce a SeqAlign file");
---
 	if (believe_query == FALSE && (myargs[14].strvalue || align_view == 10 || align_view == 11)) 
 	   ErrPostEx(SEV_FATAL, 1, 0, "-J option must be TRUE to produce a SeqAlign file");
85c88
< 	options = BLASTOptionNew(blast_program, TRUE);
---
 	options = BLASTOptionNewEx(blast_program, TRUE, TRUE);
        
129,130d131
< 	options->cutoff_s = (options->wordsize + 4)*options->reward;
< 

139,140d139
<         lcase_masking = (Boolean) myargs[28].intvalue;
< /*        
144,145c143,147
<          Allow dynamic programming gapped extension only with affine 
<            gap scores 
---
         options->mb_disc_type = myargs[38].intvalue;
 #endif
         lcase_masking = (Boolean) myargs[28].intvalue;
         /* Allow dynamic programming gapped extension only with affine 
            gap scores */
149,151d150
<         options->mb_disc_type = myargs[38].intvalue;
< #endif
< */
182d180
< 	options->is_megablast_search = TRUE;
189a188,189
         options->hsp_num_max = myargs[39].intvalue;
 
195c195
< 	query_bsp_array = (BioseqPtr PNTR) MemNew(MAX_NUM_QUERIES*sizeof(BioseqPtr));
---
 	query_bsp_array = (BioseqPtr PNTR) MemNew((MAX_NUM_QUERIES+1)*sizeof(BioseqPtr));
232c232,240
<               ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open output file %s\n", "blastngp.sat");
---
               ErrPostEx(SEV_FATAL, 1, 0, "blast: Unable to open output file %s\n", "blastngp.sat");
               return 1;
            }
         }
     	else if (align_view == 10 || align_view == 11)
     	{
         	const char* mode = (align_view == 10) ? "w" : "wb";
         	if ((aip = AsnIoOpen (blast_outputfile, (char*) mode)) == NULL) {
                 	ErrPostEx(SEV_FATAL, 1, 0, "blast: Unable to open output file %s\n", myargs[20].strvalue);
237,238d244
<         if (align_view == 7)
<            xml_aip = AsnIoOpen(blast_outputfile, "wx");
239a246,248
         if (align_view == 7) {
            xml_aip = AsnIoOpen(blast_outputfile, "wx");
         }
249,259d257
<         /* Geo extra params: */
<         if (myargs[35].intvalue > 0) {
<             options->first_db_seq=myargs[35].intvalue;
<             db_skipto=myargs[35].intvalue;
<             }
<         max_overhang=myargs[33].intvalue;
<         min_overlap=myargs[34].intvalue;
< 
<         if (myargs[36].strvalue[0] == 'T' || myargs[36].strvalue[0] == 't' ||
<                                              myargs[36].strvalue[0] == '1')
<              slice_clustering = TRUE;
261d258
<         /* Geo extra params: */
288c285
< 		 ErrPostEx(SEV_FATAL, 0, 0, "Unable to obtain bioseq\n");
---
 		 ErrPostEx(SEV_FATAL, 1, 0, "Unable to obtain bioseq\n");
315d311
< 	   ReadDBBioseqFetchEnable ("megablast", blast_database, db_is_na, TRUE);
342,362c338
< 	    } 
<            else if (myargs[12].intvalue==MBLAST_FLTHITS) {
< 	      seqalign_array = BioseqMegaBlastEngine(query_bsp_array, blast_program,
< 						     blast_database, options,
< 						     &other_returns, &error_returns,
< 						     dummy_callback, NULL, NULL, 0,
< 						     MegaBlastPrintFltHits);
<               }
< 	   else if (myargs[12].intvalue==MBLAST_HITGAPS) {
<               gap_Info=TRUE;
<               if (dbgaps_buf==NULL)
<                   dbgaps_buf=(CharPtr) Malloc(dbgaps_bufsize + 1);
<               if (qgaps_buf==NULL) 
<                  qgaps_buf=(CharPtr) Malloc(qgaps_bufsize + 1);
< 	      seqalign_array = BioseqMegaBlastEngine(query_bsp_array, blast_program,
< 						     blast_database, options,
< 						     &other_returns, &error_returns,
< 						     dummy_callback, NULL, NULL, 0,
< 						     MegaBlastPrintGapHits);
<               }           
<            else /* if (myargs[12].intvalue==MBLAST_ALIGNMENTS) */
---
 	   } else /* if (myargs[12].intvalue==MBLAST_ALIGNMENTS) */
418,421c394
< 	if (gap_Info) {
<          MemFree(dbgaps_buf);
<          MemFree(qgaps_buf);
<          }      
---
 	      
432a406,411
 
               if (align_view == 7) {
                  mbxp = PSIXmlInit(xml_aip, "megablast", blast_database, 
                                    options, query_bsp_array[0], 0);
               }
 
433a413
 	         ReadDBBioseqFetchEnable ("megablast", blast_database, db_is_na, TRUE);
484,487c464,472
<                        BXMLPrintOutput(xml_aip, seqalign, 
<                                        options, "megablast", blast_database, 
<                                        query_bsp_array[index], other_returns, 0, NULL);
<                        AsnIoReset(xml_aip);
---
                        IterationPtr iterp;
 
                        iterp = BXMLBuildOneQueryIteration(seqalign, 
                                   NULL, FALSE, 
                                   !options->gapped_calculation, index, 
                                   NULL, query_bsp_array[index]);
                        IterationAsnWrite(iterp, mbxp->aip, mbxp->atp);
                        AsnIoFlush(mbxp->aip);
                        IterationFree(iterp);
523a509,513
                  ReadDBBioseqFetchDisable();
               } 
 
               if (mbxp != NULL) {
                  MBXmlClose(mbxp, other_returns, !options->gapped_calculation);
566,567d555
<               
<               ReadDBBioseqFetchDisable();
615,616c603,604
<         if (align_view == 7)
<            xml_aip = AsnIoClose(xml_aip);
---
         /*if (align_view == 7)
           xml_aip = AsnIoClose(xml_aip);*/
626a615
         FileClose(outfp);
