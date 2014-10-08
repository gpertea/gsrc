Int2 Main (void)
 
{
	AsnIoPtr aip, xml_aip;
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
	Uint1 align_type, align_view;
	Uint4 align_options, print_options;
	ValNodePtr mask_loc, mask_loc_start, next_mask_loc;
        ValNodePtr vnp, other_returns, error_returns;

	CharPtr blast_program, blast_database, blast_inputfile, blast_outputfile;
	FILE *infp, *outfp, *mqfp;
	Int4 index, num_bsps, total_length, total_processed = 0;
	Int2 ctr = 1;
	Char prefix[2];
	SeqLocPtr last_mask, mask_slp;
	Boolean done, first_seq = TRUE, hits_found;
	CharPtr masked_query_file;
        const char *dummystr;
        Boolean lcase_masking;

        if (! GetArgs ("megablast", NUMARG, myargs))
	   return (1);

	UseLocalAsnloadDataAndErrMsg ();

	if (! SeqEntryLoad())
		return 1;

	ErrSetMessageLevel(SEV_WARNING);

	blast_program = "blastn";
        blast_database = myargs [0].strvalue;
        blast_inputfile = myargs [1].strvalue;
        blast_outputfile = myargs [4].strvalue;
	if (myargs[21].intvalue)
		html = TRUE;

	if ((infp = FileOpen(blast_inputfile, "r")) == NULL) {
	   ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open input file %s\n", blast_inputfile);
	   return (1);
	}

	outfp = NULL;
	if (blast_outputfile != NULL) {
	   if ((outfp = FileOpen(blast_outputfile, "w")) == NULL) {
	      ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open output file %s\n", blast_outputfile);
	      return (1);
	   }
	}

	align_view = (Int1) myargs[3].intvalue;

	align_type = BlastGetTypes(blast_program, &query_is_na, &db_is_na);
 
        if (myargs[15].strvalue) {
           if (myargs[15].strvalue[0] == 'f' || myargs[15].strvalue[0] == 'F' ||
               myargs[15].strvalue[0] == '0')
              believe_query = FALSE;
           else
              believe_query = TRUE;
        } else if (myargs[12].intvalue == MBLAST_ALIGNMENTS && 
                   !myargs[14].strvalue)
           believe_query = FALSE;
        else
           believe_query = TRUE;
        
        
	if (believe_query == FALSE && myargs[14].strvalue) 
	   ErrPostEx(SEV_FATAL, 0, 0, "-J option must be TRUE to produce a SeqAlign file");

	options = BLASTOptionNew(blast_program, TRUE);
	if (options == NULL)
		return 3;

	options->do_sum_stats = FALSE;
	options->is_neighboring = FALSE;
        options->expect_value  = (Nlm_FloatHi) myargs [2].floatvalue;
	number_of_descriptions = myargs[10].intvalue;	
	number_of_alignments = myargs[11].intvalue;	
	options->hitlist_size = MAX(number_of_descriptions, number_of_alignments);

	if (myargs[6].intvalue != 0)
           options->gap_x_dropoff = myargs[6].intvalue;
	if (myargs[33].intvalue != 0)
           options->dropoff_2nd_pass = myargs[33].intvalue;
        if (myargs[34].intvalue != 0)
           options->gap_x_dropoff_final = myargs[34].intvalue;

	if (StringICmp(myargs[5].strvalue, "T") == 0)
	   options->filter_string = StringSave("D");
	else
	   options->filter_string = StringSave(myargs[5].strvalue);
	
	show_gi = (Boolean) myargs[7].intvalue;
	options->penalty = myargs[8].intvalue;
	options->reward = myargs[9].intvalue;
	options->gap_open = myargs[23].intvalue;
	options->gap_extend = myargs[24].intvalue;

	if (options->gap_open == 0 && options->reward % 2 == 0 && 
	    options->gap_extend == options->reward / 2 - options->penalty)
	   /* This is the default value */
	   options->gap_extend = 0;

	options->genetic_code = 1;
	options->db_genetic_code = 1; /* Default; it's not needed here anyway */
	options->number_of_cpus = myargs[13].intvalue;
	if (myargs[17].intvalue != 0)
           options->wordsize = myargs[17].intvalue;
        if (myargs[25].intvalue == 0)
           options->cutoff_s2 = options->wordsize*options->reward;
        else 
           options->cutoff_s2 = myargs[25].intvalue;

	options->cutoff_s = (options->wordsize + 4)*options->reward;

        options->db_length = (Int8) myargs[18].floatvalue;

	options->perform_culling = FALSE;
	/* Kludge */
        options->block_width  = myargs[19].intvalue;

	options->strand_option = myargs[20].intvalue;
        options->window_size = myargs[32].intvalue;
        lcase_masking = (Boolean) myargs[28].intvalue;
/*        
#ifdef DO_NOT_SUPPRESS_BLAST_OP
        options->mb_template_length = myargs[35].intvalue;
        options->mb_one_base_step = (Boolean) myargs[36].intvalue;
         Allow dynamic programming gapped extension only with affine 
           gap scores 
        if (options->gap_open != 0 || options->gap_extend != 0)
           options->mb_use_dyn_prog = (Boolean) myargs[37].intvalue;

        options->mb_disc_type = myargs[38].intvalue;
#endif
*/
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

	if (myargs[22].strvalue)
	   options->gifile = StringSave(myargs[22].strvalue);
   
	options->is_megablast_search = TRUE;
	if (myargs[12].intvalue == MBLAST_ENDPOINTS)
	   options->no_traceback = TRUE;
	else
	   options->no_traceback = FALSE;

	options->megablast_full_deflines = (Boolean) myargs[27].intvalue;
        options->perc_identity = (FloatLo) myargs[30].floatvalue;
        if (!believe_query)
           options->megablast_full_deflines = TRUE;
        /*if (options->megablast_full_deflines)
          believe_query = FALSE;*/

	query_bsp_array = (BioseqPtr PNTR) MemNew(MAX_NUM_QUERIES*sizeof(BioseqPtr));
	sepp = (SeqEntryPtr PNTR) MemNew(MAX_NUM_QUERIES*sizeof(SeqEntryPtr));

	StrCpy(prefix, "");

	global_fp = outfp;
        if (myargs[12].intvalue != MBLAST_ALIGNMENTS)
           options->output = outfp;

	if (myargs[12].intvalue==MBLAST_ALIGNMENTS) {
	   if (align_view < 7) {
              if (html) {
                 fprintf(outfp, "<HTML>\n<TITLE>MEGABLAST Search Results</TITLE>\n");
                 fprintf(outfp, "<BODY BGCOLOR=\"#FFFFFF\" LINK=\"#0000FF\" "
                         "VLINK=\"#660099\" ALINK=\"#660099\">\n");
                 fprintf(outfp, "<PRE>\n");
              }
              init_buff_ex(90);
              BlastPrintVersionInfo("megablast", html, outfp);
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
        if (myargs[14].strvalue != NULL) {
           if ((aip = AsnIoOpen (myargs[14].strvalue,"w")) == NULL) {
              ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open output file %s\n", "blastngp.sat");
              return 1;
           }
        }

        if (align_view == 7)
           xml_aip = AsnIoOpen(blast_outputfile, "wx");


        if (myargs[31].strvalue) {      
            CharPtr delimiters = " ,;";
            CharPtr location;
            location = myargs[31].strvalue;
            options->required_start = 
                atoi(StringTokMT(location, delimiters, &location)) - 1;
            options->required_end = atoi(location) - 1;
        }
        /* Geo extra params: */
        if (myargs[35].intvalue > 0) {
            options->first_db_seq=myargs[35].intvalue;
            db_skipto=myargs[35].intvalue;
            }
        max_overhang=myargs[33].intvalue;
        min_overlap=myargs[34].intvalue;

        if (myargs[36].strvalue[0] == 'T' || myargs[36].strvalue[0] == 't' ||
                                             myargs[36].strvalue[0] == '1')
             slice_clustering = TRUE;

        /* Geo extra params: */
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

	      if (query_bsp == NULL) {
		 ErrPostEx(SEV_FATAL, 0, 0, "Unable to obtain bioseq\n");
		 return 2;
	      }
	      
	      source = BioSourceNew();
	      source->org = OrgRefNew();
	      source->org->orgname = OrgNameNew();
	      source->org->orgname->gcode = options->genetic_code;
	      ValNodeAddPointer(&(query_bsp->descr), Seq_descr_source, source);
	      
	      query_bsp_array[num_bsps++] = query_bsp;
	      
	      total_length += query_bsp->length;
	      if (total_length > myargs[16].intvalue || 
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
	   
	   ReadDBBioseqFetchEnable ("megablast", blast_database, db_is_na, TRUE);
#if 0
	   fprintf(stderr, "Process %d queries with total length %ld\n", 
		   num_bsps, total_length);
#endif
	   if (myargs[12].intvalue==MBLAST_ENDPOINTS) 
	      seqalign_array = BioseqMegaBlastEngine(query_bsp_array, blast_program,
						     blast_database, options,
						     &other_returns, &error_returns,
						     dummy_callback, NULL, NULL, 0, 
						     MegaBlastPrintEndpoints);
	   else if (myargs[12].intvalue==MBLAST_SEGMENTS) 
	      seqalign_array = BioseqMegaBlastEngine(query_bsp_array, blast_program,
						     blast_database, options,
						     &other_returns, &error_returns,
						     dummy_callback, NULL, NULL, 0,
						     MegaBlastPrintSegments);
	   else if (myargs[12].intvalue==MBLAST_ALIGN_INFO) {
              PrintTabularOutputHeader(blast_database, 
                                       (num_bsps==1) ? query_bsp_array[0] : NULL,
                                       NULL, "megablast", 0, believe_query,
                                       global_fp);
	      seqalign_array = BioseqMegaBlastEngine(query_bsp_array, blast_program,
						     blast_database, options,
						     &other_returns, &error_returns,
						     dummy_callback, NULL, NULL, 0,
						     MegaBlastPrintAlignInfo);
	    } 
           else if (myargs[12].intvalue==MBLAST_FLTHITS) {
	      seqalign_array = BioseqMegaBlastEngine(query_bsp_array, blast_program,
						     blast_database, options,
						     &other_returns, &error_returns,
						     dummy_callback, NULL, NULL, 0,
						     MegaBlastPrintFltHits);
              }
	   else if (myargs[12].intvalue==MBLAST_HITGAPS) {
              gap_Info=TRUE;
              if (dbgaps_buf==NULL)
                  dbgaps_buf=(CharPtr) Malloc(dbgaps_bufsize + 1);
              if (qgaps_buf==NULL) 
                 qgaps_buf=(CharPtr) Malloc(qgaps_bufsize + 1);
	      seqalign_array = BioseqMegaBlastEngine(query_bsp_array, blast_program,
						     blast_database, options,
						     &other_returns, &error_returns,
						     dummy_callback, NULL, NULL, 0,
						     MegaBlastPrintGapHits);
              }           
           else /* if (myargs[12].intvalue==MBLAST_ALIGNMENTS) */
	      seqalign_array = BioseqMegaBlastEngine(query_bsp_array, blast_program,
				  blast_database, options, &other_returns, 
                                  &error_returns, align_view < 7 ? tick_callback : NULL,
                                  NULL, NULL, 0, NULL);
	   
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
              
              
	   if (myargs[12].intvalue==MBLAST_ALIGNMENTS) {
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
	if (gap_Info) {
         MemFree(dbgaps_buf);
         MemFree(qgaps_buf);
         }      
#ifdef OS_UNIX
              if(align_view < 7) {
                 fprintf(global_fp, "%s\n", " done");
              }
#endif
	      
              masked_query_file = myargs[26].strvalue;
              hits_found = FALSE;

              mask_loc_start = next_mask_loc = mask_loc;
              mask_loc = NULL;
              if (seqalign_array) {
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
                    if (masked_query_file) {
                       mask_slp = MaskSeqLocFromSeqAlign(seqalign);
                       if (mask_loc) 
                          mask_slp = blastMergeFilterLocs(mask_slp, 
                              (SeqLocPtr)mask_loc->data.ptrvalue,
                              FALSE, 0, 0);
                       if (first_seq)
                          mqfp = FileOpen(masked_query_file, "w");
                       PrintMaskedSequence(query_bsp_array[index], mask_slp,
                                           mqfp, 50, lcase_masking);
                       first_seq = FALSE;
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

                       BlastPrintTabulatedResults(seqalign, 
                          query_bsp_array[index], NULL, number_of_alignments,
                          blast_program, !options->gapped_calculation, 
                          believe_query, options->required_start, 0, 
                          global_fp, (align_view == 9));

                       SeqAlignSetFree(seqalign);
                       mask_loc = MemFree(mask_loc);
                       continue;
                    } else if(align_view == 7) {
                       BXMLPrintOutput(xml_aip, seqalign, 
                                       options, "megablast", blast_database, 
                                       query_bsp_array[index], other_returns, 0, NULL);
                       AsnIoReset(xml_aip);
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
              } 

              if (masked_query_file)
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
              
              ReadDBBioseqFetchDisable();
	   } else {
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

        if (align_view == 7)
           xml_aip = AsnIoClose(xml_aip);

        if (align_view < 7 && html) 
           fprintf(outfp, "</PRE>\n</BODY>\n</HTML>\n");
        if (align_view < 7 && myargs[29].intvalue)
           fprintf(outfp, "Mega BLAST run finished, processed %d queries\n",
                   total_processed);
	MemFree(query_bsp_array);
	MemFree(sepp);
	options = BLASTOptionDelete(options);
	FileClose(infp);
	
	return 0;
}
