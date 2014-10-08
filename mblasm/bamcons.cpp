#include "GArgs.h"
#include "GStr.h"
#include "GHash.hh"
#include "GVec.hh"
#include "GList.hh"
#include "GapAssem.h"
#include "GFaSeqGet.h"
#include "GFastaIndex.h"
#include "GBam.h"
#define USAGE "Usage:\n\
 bamcons <file.sorted.bam> [-r <genome.fa>] [-c <clipmax[%]>] \\ \n\
    [-o <outfile.ace>]\n\
   The input BAM file <file.sorted.bam> must be sorted by genomic location.\n\
   \n\
   -r <genome.fa> is a multi-fasta file with all the genomic sequences,\n\
      preferably indexed with 'samtools faidx'\n\
   -c maximum clipping allowed for component sequences when\n\
      mapped onto a reference; if it ends with '%' then maximum\n\
      clipping is estimated for each read as <clipmax> percent\n\
      of its length\n\
   -o the output ACE file is written to <outfile.ace> instead of stdout\n\
   -G do not remove consensus gaps (default is to edit reads\n\
      by removing gap-dominated columns in the MSA)\n"

#define LOG_MSG_CLIPMAX "Overlap between %s and reference %s rejected due to clipmax=%4.2f constraint.\n"
#define LOG_MSG_OVLCLIP "Overlap between %s and %s invalidated by the %dnt clipping of %s at %d' end.\n"
//-------- global variables 
bool debugMode = false;
bool removeConsGaps = false;
bool verbose = false;
int rlineno = 0;

static FILE* outf;
//static GHash<GASeq> seqs(false);
//static GList<GSeqAlign> alns(true, true, false);
// sorted, free element, not unique

// flag for reference status:
static const unsigned char flag_IS_REF = 0;

//this is only for tree/transitive embedding:
static const unsigned char flag_HAS_PARENT = 1;

float clipmax = 0;

//genomic fasta sequence handling
class GFastaHandler {
public:
	char* fastaPath;
	GFastaIndex* faIdx;
	GFastaHandler(const char* fpath = NULL) :
			fastaPath(NULL), faIdx(NULL) {
		if (fpath != NULL && fpath[0] != 0)
			init(fpath);
	}

	void init(const char* fpath) {
		if (fpath == NULL || fpath[0] == 0)
			return;
		if (!fileExists(fpath))
			GError("Error: file/directory %s does not exist!\n", fpath);
		fastaPath = Gstrdup(fpath);
		if (fastaPath != NULL) {
			if (fileExists(fastaPath) > 1) { //exists and it's not a directory
				GStr fainame(fastaPath);
				//the .fai name might have been given directly
				if (fainame.rindex(".fai") == fainame.length() - 4) {
					//.fai index file given directly
					fastaPath[fainame.length() - 4] = 0;
					if (!fileExists(fastaPath))
						GError("Error: cannot find fasta file for index %s !\n", fastaPath);
				} else
					fainame.append(".fai");
				//fainame.append(".fai");
				faIdx = new GFastaIndex(fastaPath, fainame.chars());
				GStr fainamecwd(fainame);
				int ip = -1;
				if ((ip = fainamecwd.rindex('/')) >= 0)
					fainamecwd.cut(0, ip + 1);
				if (!faIdx->hasIndex()) { //could not load index
					//try current directory
					if (fainame != fainamecwd) {
						if (fileExists(fainamecwd.chars()) > 1) {
							faIdx->loadIndex(fainamecwd.chars());
						}
					}
				} //tried to load index
				if (!faIdx->hasIndex()) {
					GMessage(
					    "No fasta index found for %s. Building it now, please wait..\n",
					    fastaPath);
					faIdx->buildIndex();
					if (faIdx->getCount() == 0)
						GError("Error: no fasta records found!\n");
					GMessage("FASTA index rebuilt.\n");
					FILE* fcreate = fopen(fainame.chars(), "w");
					if (fcreate == NULL) {
						GMessage(
						    "Warning: cannot create FASTA index file %s! (permissions?)\n",
						    fainame.chars());
						if (fainame != fainamecwd)
							fcreate = fopen(fainamecwd.chars(), "w");
						if (fcreate == NULL)
							GError("Error: cannot create fasta index %s!\n",
							    fainamecwd.chars());
					}
					if (faIdx->storeIndex(fcreate) < faIdx->getCount())
						GMessage("Warning: error writing the index file!\n");
				} //index created and attempted to store it
			} //multi-fasta
		} //genomic sequence given
	}

	char* fetchSeq(const char* gseqname, int& seqlen) {
		if (fastaPath == NULL || faIdx == NULL || gseqname == NULL)
			return NULL;
		//genomic sequence given
		seqlen=0;
		GFastaRec* farec = faIdx->getRecord(gseqname);
		if (farec != NULL) {
			 GFaSeqGet faseq(fastaPath, farec->seqlen, farec->fpos,
			    farec->line_len, farec->line_blen);
			return faseq.fetchSeq(&seqlen); //allocate a ptr and load the whole sequence, bypassing the buffer
		} else {
			GMessage("Warning: couldn't find fasta record for '%s'!\n", gseqname);
			return NULL;
		}
	}
	~GFastaHandler() {
		GFREE(fastaPath);
		delete faIdx;
	}
};

//--------------------------------
class RefAlign {
	char* linecpy;
	int lineno; //input file line#
	char* gpos;
	char* gpos_onref;
	char* gaps; //char string with gaps in this read
	char* gaps_onref; //char string with gaps in reference
public:
	char* seqname; //aligned read
	int seqlen; //length of this read
	int offset; //offset of this read relative to reference
	int clip5;
	int clip3;
	char reverse; //is this reversed?
	void parseErr(int fldno);
	RefAlign(char* line, int len, int lno);
	~RefAlign();
	int nextRefGap(int& pos);
	int nextSeqGap(int& pos);
};

// bundle data structure, holds all data needed for
// infering contigs from a bundle
enum BundleStatus {
	BUNDLE_STATUS_CLEAR=0, //available for loading/prepping
	BUNDLE_STATUS_LOADING, //being prepared by the main thread (there can be only one)
	BUNDLE_STATUS_READY //ready to be processed, or being processed
};


struct CReadAln:public GSeg {
	//DEBUG ONLY:
	GStr name;
	// 0: strand; 1: NH; 2: pair's no; 3: coords of read; 4: junctions
	char strand; // 1, 0 (unkown), -1 (reverse)
	short int nh;
	//int pair_idx; //mate index alignment in CReadAln list
	//GVec<GSeg> segs; //"exons"
	GStr rname;

	//DEBUG ONLY: (discard rname when no debugging needed)
	CReadAln(char _strand=0, short int _nh=0,
			int rstart=0, int rend=0 /*,  const char* rname=NULL */): GSeg(rstart, rend), //name(rname),
					strand(_strand) { }
};


struct BundleData { //clustered (contiguous) read alignments
 BundleStatus status;
 int idx; //index in the main bundles array
 int start;
 int end;
 GStr refseq;
 GList<GSeqAlign> alns;
 GHash<GASeq> seqs;
 GVec<float> bpcov;
 BundleData():status(BUNDLE_STATUS_CLEAR), idx(0), start(0), end(0),
		 refseq(), alns(false,true), seqs(false), bpcov(1024) { }

 void getReady(int currentstart, int currentend) {
	 start=currentstart;
	 end=currentend;
	 status=BUNDLE_STATUS_READY;
 }

 void Clear() {
	alns.Clear();
	seqs.Clear();
	bpcov.Clear();
	bpcov.setCapacity(1024);
	start=0;
	end=0;
	status=BUNDLE_STATUS_CLEAR;
 }

 ~BundleData() {
	Clear();
 }
};


void processBundle(BundleData* bundle) {
	if (verbose) {
		//printTime(stderr);
		GMessage(">bundle %s:%d-%d(%d) begins processing...\n",
				bundle->refseq.chars(), bundle->start, bundle->end, bundle->alns.Count());
	}
	//TODO: MSA processing and printing here
	for (int i = 0; i < bundle->alns.Count(); i++) {
		GSeqAlign* a = bundle->alns.Get(i);
		//loadAlnSeqs(a, refcdb); //,cdbyank, refcdb); //loading actual sequences for this cluster
		if (debugMode) { //write plain text alignment file
			fprintf(outf, ">Alignment%d (%d)\n", i + 1, a->Count());
			a->print(outf, 'v');
		} else { //write a real ACE file
			//a->buildMSA();
			GStr ctgname;
			ctgname.format("AorContig%d", i + 1);
			a->writeACE(outf, ctgname.chars());
			a->freeMSA(); //free MSA and seq memory
		}
	} // for each PMSA cluster
	// oooooooooo D O N E oooooooooooo
	//alns.Clear();
	//seqs.Clear();

	if (verbose) {
	  //printTime(stderr);
	  GMessage("^bundle %s:%d-%d(%d) done.\n",bundle->refseq.chars(),
	  		bundle->start, bundle->end, bundle->alns.Count());
	}
	bundle->Clear();
}



//void loadAlnSeqs(GSeqAlign* aln, GFastaHandler* refcdb = NULL); //, GCdbYank* cdbynk, GCdbYank* refcdb=NULL);
//void loadRefSeq(GSeqAlign* aln, GFastaHandler* refcdb);
void printDebugAln(FILE* f, GSeqAlign* aln, int num, GFastaHandler* refcdb = NULL); // GCdbYank* cdbynk, GCdbYank* refcdb=NULL);

//returns the end of a space delimited token
void skipSp(char*& p) {
	while (*p == ' ' || *p == '\t' || *p == '\n')
		p++;
}
char* endSpToken(char* str) {
	if (str == NULL || *str == 0)
		return NULL;
	char* p = str;
	for (; *p != 0; p++)
		if (*p == ' ' || *p == '\t' || *p == '\n')
			return p;
	return p; // *p is '\0'
}


//-- prepareMerge checks clipping and even when no clipmax is given,
//   adjusts clipping as appropriately

//========================================================
//====================     main      =====================
//========================================================
int main(int argc, char * const argv[]) {
	//GArgs args(argc, argv, "DGvd:o:c:");
	GArgs args(argc, argv, "DGvr:o:c:");
	GFastaHandler* refcdb = NULL;
	int e;
	if ((e = args.isError()) > 0)
		GError("%s\nInvalid argument: %s\n", USAGE, argv[e]);
	if (args.getOpt('h') != NULL)
		GError("%s\n", USAGE);
	debugMode = (args.getOpt('D') != NULL);
	removeConsGaps = (args.getOpt('G') == NULL);
	verbose = (args.getOpt('v') != NULL);
	if (debugMode)
		verbose = true;

	MSAColumns::removeConsGaps = removeConsGaps;
	GStr infile;
	if (args.startNonOpt()) {
		infile = args.nextNonOpt();
		//GMessage("Given file: %s\n",infile.chars());
	}
	else
		GError("%s\nInput BAM file required!", USAGE);
	//==
	/*
	FILE* inf;
	if (!infile.is_empty()) {
		inf = fopen(infile, "r");
		if (inf == NULL)
			GError("Cannot open input file %s!\n", infile.chars());
	} else
		inf = stdin;
  */
	GStr s = args.getOpt('c');
	if (!s.is_empty()) {
		bool ispercent = (s[-1] == '%');
		if (ispercent)
			s.trimR('%');
		int c = s.asInt();
		if (c <= 0)
			GError("Error: invalid -c <clipmax> (%d) option provided (must be "
					"a positive integer)!\n", c);
		if (ispercent && c > 99)
			GError("Error: invalid percent value (%d) for -c option "
					" (must be an integer between 1 and 99)!\n", c);
		if (ispercent) {
			clipmax = float(c) / 100;
			int maxp = iround(100 * clipmax);
			if (verbose)
				fprintf(stderr, "Percentual max clipping set to %d%%\n", maxp);
		} else {
			clipmax = c;
			if (verbose)
				fprintf(stderr, "Max clipping set to %d bases\n", c);
		}

	} //clipmax option
	GStr outfile = args.getOpt('o');
	if (!outfile.is_empty()) {
		outf = fopen(outfile, "w");
		if (outf == NULL)
			GError("Cannot open file %s for writing!\n", outfile.chars());
	} else
		outf = stdout;
	//************************************************

	/*
	 GCdbYank* refcdb = NULL;
	 dbidx=args.getOpt('r');
	 if (!dbidx.is_empty())
	 refcdb=new GCdbYank(dbidx.chars());
	 */
	s = args.getOpt('r');
	if (!s.is_empty()) {
		refcdb = new GFastaHandler(s.chars());
	}

	GBamReader bamreader(infile.chars());
	GStr lastref;
	int currentstart=0, currentend=0;

	BundleData bundles[1];
	BundleData* bundle = &(bundles[0]);
	int readthr=3;      // read coverage per bundle bp to accept it; otherwise considered noise
	uint bundledist=0;  // reads at what distance should be considered part of separate bundles


	GBamRecord* brec=NULL;
	bool more_alns=true;
	GASeq* refseq=NULL;
	while (more_alns) {
		 bool chr_changed=false;
		 int pos=0;
		 const char* rname=NULL;
		 int nh=1;
		 int hi=0;
		 bool new_bundle=false;
		 delete brec;
		 if ((brec=bamreader.next())!=NULL) {
			 if (brec->isUnmapped()) continue;
			 rname=brec->refName();
			 if (rname==NULL) GError("Error: cannot retrieve target seq name from BAM record!\n");
			 pos=brec->start; //BAM is 0 based, but GBamRecord makes it 1-based
			 chr_changed=(lastref.is_empty() || lastref!=rname);
			 if (chr_changed) {
				 //TODO: load chr sequence here
				 refseq=new GASeq(lastref.chars(), 0, 0);
				 refseq->setFlag(flag_IS_REF);
			 }
			 if (!chr_changed && currentend>0 && pos>currentend+(int)bundledist)
				   new_bundle=true;
		 }
		 else { //no more alignments
			 more_alns=false;
			 new_bundle=true; //fake a new start (end of last bundle)
		 }
		 if (new_bundle || chr_changed) {
			 //hashread.Clear();
			 if (bundle->alns.Count()>0) { // process reads in previous bundle
				// geneno=infer_transcripts(geneno, lastref, $label,\@readlist,$readthr,\@junction,$junctionthr,$mintranscriptlen,\@keepguides);
				// (readthr, junctionthr, mintranscriptlen are globals)
				bundle->getReady(currentstart, currentend);
				//TODO: processBundle(bundle);
			 } //have alignments to process
			 else { //no read alignments in this bundle?
				bundle->Clear();
			 }
			 if (chr_changed) {
				 lastref=rname;
				 currentend=0;
			 }
			 if (!more_alns) {
					if (verbose) {
						//printTime(stderr);
						GMessage(" Done reading alignments.\n");
					}
				 //noMoreBundles();
				 break;
			 }
			currentstart=pos;
			bundle->refseq=lastref;
			bundle->start=currentstart;
			bundle->end=currentend;
		 } //<---- new bundle
		 //TODO:
		 //currentend=processRead(currentstart, currentend, *bundle, *brec);
	} //for each read alignment

	 //cleaning up
	 delete brec;
	 bamreader.bclose();
	 if (verbose) {
	    //printTime(stderr);
	    GMessage(" Done.\n");
	 }

	/*************** now build MSAs and write them
	if (verbose) {
		fprintf(stderr, "\n%d input lines processed.\n", rlineno);
		fprintf(stderr, "Refining and printing %d MSA(s)..\n", alns.Count());
		fflush(stderr);
	}
  */
	//print all the alignments
	fflush(outf);
	if (outf != stdout)
		fclose(outf);
	delete refcdb;
	GMessage("*** all done ***\n");
#ifdef __WIN32__
	//getc(stdin);
#endif
}

//---------------------- RefAlign class
void RefAlign::parseErr(int fldno) {
	fprintf(stderr, "Error parsing input line #%d (field %d):\n%s\n", lineno,
	    fldno, linecpy);
	exit(3);
}

RefAlign::RefAlign(char* line, int len, int lno) {
	lineno = lno;
	linecpy = Gstrdup(line);
	seqname = line;
	char* p = endSpToken(seqname);
	if (p == NULL)
		parseErr(0);
	*p = 0;
	p++;
	skipSp(p);
	if (*p != '-' && *p != '+')
		parseErr(1);
	reverse = 0;
	if (*p == '-')
		reverse = 1;
	p++;
	clip5 = 0;
	clip3 = 0;
	if (!parseInt(p, seqlen))
		parseErr(2);
	if (!parseInt(p, offset))
		parseErr(3);
	offset--;
	skipSp(p);
	if (isdigit(*p)) {  //clipping follows
		if (!parseInt(p, clip5))
			parseErr(4);
		if (!parseInt(p, clip3))
			parseErr(5);
	}
	if (reverse != 0)
		Gswap(clip5, clip3);
	skipSp(p);
	gaps = NULL;
	gaps_onref = NULL;
	gpos = NULL;
	gpos_onref = NULL;
	if (*p == 0)
		return;
	while (*p != 0 && *p != 'R') {
		while (*p != ' ' && *p != '\t' && *p != 0)
			p++;
		skipSp(p);
	}
	if (*p == 0)
		return;
	// it is R
	p++;
	if (*p != ':')
		parseErr(6); //invalid field
	p++;
	if (isdigit(*p)) { //should be a number here
		gaps = p;
		gpos = p;
	} else if (*p != '/')
		parseErr(6);
	while (*p != 0 && *p != ' ' && *p != '/' && *p != '\t')
		p++;
	if (*p == 0)
		return;
	if (*p == '/') {
		*p = 0; //split here
		p++;
		gaps_onref = p;
		gpos_onref = p;
	}
}
RefAlign::~RefAlign() {
	GFREE(linecpy); //that was just for error messages
}

int RefAlign::nextSeqGap(int& pos) {
	int r = 1;
	if (gpos == NULL || *gpos == 0)
		return 0;
	if (!parseInt(gpos, pos))
		parseErr(6);
	if (pos <= 0)
		parseErr(6);
	if (*gpos == '+') { //gap length following
		gpos++;
		if (!parseInt(gpos, r))
			parseErr(6);
		if (r <= 0)
			parseErr(6);
	}
	if (*gpos == ',')
		gpos++;
	return r;
}

int RefAlign::nextRefGap(int& pos) {
	int r = 1;
	if (gpos_onref == NULL || *gpos_onref == 0)
		return 0;
	if (!parseInt(gpos_onref, pos))
		parseErr(6);
	if (pos <= 0)
		parseErr(6);
	if (*gpos_onref == '+') { //gap length following
		gpos_onref++;
		if (!parseInt(gpos_onref, r))
			parseErr(6);
		if (r <= 0)
			parseErr(6);
	}
	if (*gpos_onref == ',')
		gpos_onref++;
	return r;
}

void printDebugAln(FILE* f, GSeqAlign* aln, int num, GFastaHandler* refcdb) { //, GCdbYank* cdbynk, GCdbYank* refcdb) {
	fprintf(f, ">[%d]DebugAlign%d (%d)\n", rlineno + 1, num, aln->Count());
	//loadAlnSeqs(aln, refcdb); //,cdbynk, refcdb);
	aln->print(f, '=');
}

void loadRefSeq(GASeq& s, GFastaHandler* refcdb) { //, GCdbYank* cdbynk, GCdbYank* refcdb) {
	if (s.len == 0 && refcdb) {  //reference seq not loaded yet
		//fetch it from the fasta index
		int seqlen=0;
		 char* seq=refcdb->fetchSeq(s.id, seqlen);
		 if (seq) {
			 s.setSeqPtr(seq, seqlen);
			 s.allupper();
			 s.loadProcessing();
			}
	}
}

