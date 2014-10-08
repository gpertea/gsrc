#!/usr/bin/perl
use strict;

# this version is similar to v6 but it uses the transfrags coming into the child/or parent to decide the continuation path, in a different way than v10: v10 computes overall max_compon based on transfrag values, while this version computes max_compon based on 'in's and then on 'out's (using intrcov, and outtrcov instead of transfrag values) or viceversa (depending if I am extending path for child or for parent). the goal is to model more closely the v4 output, but in a more systematic implementation than v4

my ($bamfile,$label,$guidegff)=@ARGV;

## Global parameters
my $epsilon=0.0000000001;
my $trthr=1;   # transfrag pattern threshold 


### Input parameters ###

my $mintranscriptlen=200; # minimum number for a transcript to be printed
#my $junctionsupport=20; # anchor length for junction to be considered well supported <- consider shorter??
my $junctionsupport=10; # anchor length for junction to be considered well supported <- consider shorter??
#my $junctionthr=3; # number of reads needed to support a particular junction
my $junctionthr=2; # number of reads needed to support a particular junction
my $readthr=3;     # read coverage per bundle bp to accept it; otherwise considered noise
my $bundledist=0; # reads at what distance should be considered part of separate bundles <- this is not addressed everywhere, e.g. in infer_transcripts -> look into this
if(!$label) {  $label="RLINK";}

my $guided=0;
if($guidegff) { $guided=1;}

### End parameters


my %refguides; # 0: strand; 1: coords on genome

if($guided) {

  my %ng;
  my %parent;

  # load gff and sort
  open(F,$guidegff);
  while(<F>) {
    if(/\texon\t/) {
      my ($ref,$seq,$type,$start,$end,$val,$strand,@rest)=split;
      my ($parentname)=/Parent=([^;]+)/;
      if(!defined $parent{$parentname}) {
	if(!defined $ng{$ref}) { $ng{$ref}=0;}
	$parent{$parentname}=$ng{$ref}++;
	if($strand eq "-") { $refguides{$ref}[$parent{$parentname}][0]=0;}
	elsif($strand eq "+") { $refguides{$ref}[$parent{$parentname}][0]=2;}
	else { $refguides{$ref}[$parent{$parentname}][0]=1;}
      }
      push(@{$refguides{$ref}[$parent{$parentname}][1]},[$start,$end]);

      #print STDERR "ref=$ref parentname=$parentname strand=$strand start=$start end=$end\n";
      #print STDERR "ng{$ref}=",$ng{$ref}," computed ng=",scalar(@{$refguides{$ref}})," parent=$parentname\n";
      
      if($ng{$ref} != scalar(@{$refguides{$ref}})) { exit;}

    }
  }
  close(F);

}

open(F,"samtools view $bamfile|");

my $currentstart;
my $currentend;

my @readlist; # 0: strand; 1: NH; 2: pair's no; 3: coords of read; 4: junctions
my @junction; # 0: strand; 1: start; 2: end; 3: nreads; 4: nreads_good;
my %hashread;
my %hashjunction;
my @guides=();
my $ng_start=0;
my $ng_end=-1;
my $ng=0;

my $lastref="";
my $geneno=0;
my $ncluster=0;

while(<F>) { 

    my ($qname,$flag,$rname,$pos,$mapq,$cigar,$rnext,$pnext,@rest)=split;

    my $strand=0;
    my ($sign)=($_ =~ /XS:A:([\+\-])/);
    $sign='' unless $sign;
    if($sign eq "+") { $strand=1;}      # strand is +1 is positive, -1 if negative, and 
    elsif($sign eq "-") { $strand=-1;}  # 0 if unkown (no XS paramater present)


    my ($nh)= ($_ =~ /NH:i:(\d+)/);

    my ($hi)= ($_ =~ /HI:i:(\d+)/);
    $hi='' unless $hi;

    if($rname ne $lastref) { # start a new chromosome

	%hashread=();

	if($#readlist > -1) { # readlist is not empty

	  my @keepguides=();
	  if($ng_end>=$ng_start) { @keepguides=@guides[$ng_start..$ng_end];}
	  $geneno=infer_transcripts($geneno,$lastref,$label,\@readlist,$readthr,\@junction,$junctionthr,$mintranscriptlen,\@keepguides);
	  $ncluster++;
	  #print STDERR "Done cluster $ncluster\n";
	}

	if($guided) {
	  @guides=();
	  @guides = sort { $a->[1][0][0] <=> $b->[1][0][0] } @{$refguides{$rname}}; # check this
	  $ng=scalar(@guides);
	}

	@readlist=();
	@junction=();
	%hashjunction=();
	$lastref=$rname;
	$currentend=0;
	$currentstart=$pos;
    
	if($guided) { # adjust currentstart and currentend based on the guiding transfrags
	  $ng_start=$ng_end+1;
	  while($ng_start<$ng && $guides[$ng_start][1][-1][1] < $pos) { $ng_start++;} # skip guides that have no read coverage

	  if($ng_start<$ng && $guides[$ng_start][1][0][0]<$pos) { 
	    $currentstart=$guides[$ng_start][1][0][0];
	    $currentend=$guides[$ng_start][1][-1][1];
	    $ng_end=$ng_start;
	    while($ng_end+1<$ng && $guides[$ng_end+1][1][0][0]<=$pos) {
	      $ng_end++;
	      if($guides[$ng_end][1][-1][1]>$currentend) { 
		$currentend=$guides[$ng_end][1][-1][1];
	      }
	    }
	  }
	}
    }
    elsif($pos>$currentend+$bundledist) { # read starts after end of current "bundle"; SHOULD I ADD A -1 HERE?
      %hashread=();
      %hashjunction=(); # ??? does this get used in infer_transcripts???

      if($#readlist > -1) { # readlist is not empty

	#  print STDERR "0: junction has ",scalar(@junction)," elements\n";
	    
	my @keepguides=();
	if($ng_end>=$ng_start) { @keepguides=@guides[$ng_start..$ng_end];}
	$geneno=infer_transcripts($geneno,$lastref,$label,\@readlist,$readthr,\@junction,$junctionthr,$mintranscriptlen,\@keepguides);
	$ncluster++;
	#print STDERR "Done cluster $ncluster\n";
      }
      @readlist=();
      @junction=();
      $currentstart=$pos;
          
      if($guided) { # adjust currentstart and currentend based on the guiding transfrags
	$ng_start=$ng_end+1;
	while($ng_start<$ng && $guides[$ng_start][1][-1][1] < $pos) { $ng_start++;} # skip guides that have no read coverage
	
	if($ng_start<$ng && $guides[$ng_start][1][0][0]<$pos) { 
	  $currentstart=$guides[$ng_start][1][0][0];
	  $currentend=$guides[$ng_start][1][-1][1];
	  $ng_end=$ng_start;
	  while($ng_end+1<$ng && $guides[$ng_end+1][1][0][0]<=$pos) {
	    $ng_end++;
	    if($guides[$ng_end][1][-1][1]>$currentend) { 
	      $currentend=$guides[$ng_end][1][-1][1];
	    }
	  }
	}
      }
    }

    $currentend=process_read($currentend,\@readlist,\%hashread,\@junction,\%hashjunction,$junctionsupport,$qname,$flag,$pos,$cigar,$pnext,$strand,$nh,$hi);  # here I might update current end to be at least as big as the start of the read pair in the gragment?? -> maybe not because then I could introduce some false positives with paired reads mapped badly

    if($guided) { # I need to adjust end according to guides
      while($ng_end+1 < $ng && $guides[$ng_end+1][1][0][0]<=$currentend) {
	$ng_end++;
	if($guides[$ng_end][1][-1][1]>$currentend) { 
	  $currentend=$guides[$ng_end][1][-1][1];
	}
      }
    }

}
close(F);

# no more reads to read -> process last bundle
if($#readlist > -1) { # readlist is not empty
  my @keepguides=();
  if($ng_end>=$ng_start) { @keepguides=@guides[$ng_start..$ng_end];}
  $geneno=infer_transcripts($geneno,$lastref,$label,\@readlist,$readthr,\@junction,$junctionthr,$mintranscriptlen,\@keepguides);
  $ncluster++;
  #print STDERR "Done cluster $ncluster\n";
}

sub infer_transcripts {
    my ($geneno,$refname,$label,$readlist,$readthr,$junction,$junctionthr,$mintranscriptlen,$guides)=@_;

    #print STDERR "1: infer_transcripts on $lastref from $currentstart $currentend\n";

    #print STDERR "juntion=$junction\n";
    

    # clean unsupported junctions
    clean_junctions($junction,$junctionthr);

    #print STDERR "Cleaned junctions\n";

    $geneno=build_graphs($geneno,$refname,$label,$readlist,$readthr,$junction,$mintranscriptlen,$guides);

}

sub build_graphs {
    my ($geneno,$refname,$label,$readlist,$readthr,$junction,$mintranscriptlen,$guides)=@_;
    
    my $nreads=scalar(@{$readlist});

    ### form groups on strands
    # all groups below are like this: 0 = negative strand; 1 = unknown strand; 2 = positive strand
    my @group; # 0: start; 1: end; 2: color; 3: group no; 4: sum of covered bp; 5: address of next group 
    my $ngroup=0; # number of groups 
    my @currgroup=(0,0,0);  # current group of each type 
    my @startgroup=(0,0,0); # start group of each type

    my $color=0;    # next color to assign

    my @readgroup;  # remebers groups for each read
    my %equalcolor; # remembers colors for the same bundle
    my %merge;      # remembers merged groups

    for(my $n=0;$n<$nreads;$n++) {

	# first see if read contains deleted junctions -> if it does ignore it	
	my $keep=1; # assume I keep read
	my $i=0;
	while(defined $$readlist[$n][4][$i]) {
	    if(!$$junction[$$readlist[$n][4][$i]][0]) {
		$keep=0; # ignore read if it contains deleted junctions
		$$readlist[$n][1]=0; # set number of hits to 0
		my $np=$$readlist[$n][2];
		if($np>-1) { $$readlist[$np][2]=-1;} # delete pairing for the paired read
		last;
	    }
	    $i++;
	}

	if($keep) { # add read to groups

	    ($color,$ngroup)=add_read_to_group($n,$readlist,$color,\@group,$ngroup,\@currgroup,\@startgroup,\@readgroup,\%equalcolor,\%merge);

	}
    }

=begin

    print STDERR "$ngroup groups created!\n";
    for(my $sno=0;$sno<3;$sno++) {
	print STDERR "Groups on strand $sno:\n";
	my $group=$startgroup[$sno];
	while($group) {
	    print STDERR " gr ",$$group[3],"(",$$group[2],",",$$group[4],"): ",$$group[0],"-",$$group[1];
	    $group=$$group[5];
	}
	print STDERR "\n";
    }

=cut    

    ### form bundles here

    ## first do color assignment
    @currgroup=@startgroup;
    my @prevgroup=(0,0,0);
    
    my %eqposcol;
    my %eqnegcol;
   
    while($currgroup[0] || $currgroup[1] || $currgroup[2]) { # there are still groups to process
	
	my $nextgr=get_min_start(\@currgroup); # gets the index of currgroup with the left most begining
	
	my $grcol = $currgroup[$nextgr][2];    # set smallest color for currgroup[$nextgr]
	#print STDERR "start grcol=$grcol\n";
	while(defined $equalcolor{$grcol}) {
	    $grcol=$equalcolor{$grcol};
	}
	$currgroup[$nextgr][2]=$grcol;

	#print STDERR "nextgr=$nextgr grcol=$grcol current group is at coords: ",$currgroup[$nextgr][0],"-",$currgroup[$nextgr][1],"\n";

	if($nextgr == 0) { # negative strand group
	    if($prevgroup[1] && $currgroup[$nextgr][0]<$prevgroup[1][1]) { # overlaps an unknown strand group
		set_strandcol(\@{$prevgroup[1]},\@{$currgroup[$nextgr]},$grcol,\%eqnegcol,\%equalcolor);
	    }
	}
	elsif($nextgr == 2) { # positive strand group
	    if($prevgroup[1] && $currgroup[$nextgr][0]<$prevgroup[1][1]) { # overlaps an unknown strand group
		set_strandcol(\@{$prevgroup[1]},\@{$currgroup[$nextgr]},$grcol,\%eqposcol,\%equalcolor);
	    }
	}
	else { # unknown strand group
	    if($prevgroup[0] && $currgroup[$nextgr][0]<$prevgroup[0][1]) { # overlaps negative strand group
		set_strandcol(\@{$currgroup[$nextgr]},\@{$prevgroup[0]},$prevgroup[0][2],\%eqnegcol,\%equalcolor);
	    }
	    if($prevgroup[2] && $currgroup[$nextgr][0]<$prevgroup[2][1]) { # overlaps positive strand group
		set_strandcol(\@{$currgroup[$nextgr]},\@{$prevgroup[2]},$prevgroup[2][2],\%eqposcol,\%equalcolor);
	    }
	}

	$prevgroup[$nextgr]=$currgroup[$nextgr];
	$currgroup[$nextgr]=$currgroup[$nextgr][5];

	#if($currgroup[$nextgr]) { print STDERR "currgroup[$nextgr] is now ",$currgroup[$nextgr][0],"-",$currgroup[$nextgr][1],"\n";}
	#else { print STDERR "currgroup[$nextgr] is now 0\n";}
    }

    #print STDERR "Colors assigned!\n";

	    
    ## create bundles
    @currgroup=@startgroup;
    @prevgroup=(0,0,0);

    my @bundle; # 0: len; 1: bp coverage; 2: reference of start node; ### 3: ref to currrent node; 
    my @nbundle=(0,0,0);
    my %bundlecol; # associates a bundle number to a group color
    my %group2bundle; # to retrace reads from group no to bundle

    my @bnode; # bundlenodes (bnode's) are defined like this: 0:start; 1: end; 2: link to next bnode; 3: bnode coverage
    
    while($currgroup[0] || $currgroup[1] || $currgroup[2]) { # there are still groups to process


      #print STDERR "Process groups: ",$currgroup[0],",",$currgroup[1],",",$currgroup[2],"\n";


	my $nextgr=get_min_start(\@currgroup);

	#print STDERR "nextgr=$nextgr\n";


	# get group color; I need to redo this to ensure I equalize all colors -> they could still be hanged by set_strandcol
	my $grcol = $currgroup[$nextgr][2];
	#print STDERR "Start grcol=$grcol\n";

	while(defined $equalcolor{$grcol}) {
	    $grcol=$equalcolor{$grcol};
	}
	$currgroup[$nextgr][2]=$grcol;

	#print STDERR "grcol=$grcol at ",$currgroup[$nextgr][0],"-",$currgroup[$nextgr][1]," where group no is ",$currgroup[$nextgr][3]," and sum=",$currgroup[$nextgr][4],"\n";


	if($nextgr == 0 || $nextgr ==2 || ($nextgr==1 &&(!defined $eqnegcol{$grcol}) && (!defined $eqposcol{$grcol}))) { # negative or positive strand bundle or unstranded bundle 

	    my $bno=$bundlecol{$grcol};

	    #print STDERR "Bundle: bno=$bno\n";



	    if(defined $bno) { # bundle for group has been created before

	      #print STDERR "Add group ",$currgroup[$nextgr][3]," to bundle $bno\n";
	      
		$bnode[$nextgr][$bno]=add_group_to_bundle(\@{$currgroup[$nextgr]},\@{$bundle[$nextgr][$bno]},$bnode[$nextgr][$bno]);

=begin

		print STDERR "Current bundle[$nextgr][$bno] is: ";
	      my $bn=$bundle[$nextgr][$bno][2];
	      while($bn) {
		print STDERR " ",$$bn[0],"-",$$bn[1]," cov=",$$bn[3];
		$bn=$$bn[2];
	      }
	      print STDERR "\n";

=cut

	    }
	    else { # create new bundle
		$bno=create_bundle(\@bundle,\@nbundle,$nextgr,\@{$currgroup[$nextgr]});
		$bnode[$nextgr][$bno]=$bundle[$nextgr][$bno][2];

		#print STDERR "Added group ",$currgroup[$nextgr][3]," to new bundle $bno\n";

		$bnode[$nextgr][$bno]=$bundle[$nextgr][$bno][2];
		$bundlecol{$grcol}=$bno;
	    }
	    my $id = $currgroup[$nextgr][3].":".$nextgr;
	    $group2bundle{$id}=$bnode[$nextgr][$bno]; # gives specific bnode in bundle


	}
	else { # unknown strand

	    if(defined $eqnegcol{$grcol}){
		my $negcol=$eqnegcol{$grcol};
		while(defined $equalcolor{$negcol}) {
		    $negcol=$equalcolor{$negcol};
		}

		my $bno=$bundlecol{$negcol};
		if(defined $bno) { # bundle for group has been created before
		  #print STDERR "Add neutral group ",$currgroup[$nextgr][3]," to neg bundle $bno\n";

		    $bnode[0][$bno]=add_group_to_bundle(\@{$currgroup[$nextgr]},\@{$bundle[0][$bno]},$bnode[0][$bno]);

=begin

		    print STDERR "Current bundle[0][$bno] is: ";
		  my $bn=$bundle[0][$bno][2];
		  while($bn) {
		    print STDERR " ",$$bn[0],"-",$$bn[1]," cov=",$$bn[3];
		    $bn=$$bn[2];
		  }
		  print STDERR "\n";

=cut

		}
		else { # create new bundle
		    $bno=create_bundle(\@bundle,\@nbundle,0,\@{$currgroup[$nextgr]});
		    $bnode[0][$bno]=$bundle[0][$bno][2];

		    #print STDERR "Added neutral group ",$currgroup[$nextgr][3]," to new neg bundle $bno\n";

		    $bundlecol{$negcol}=$bno;
		}
		my $id = $currgroup[$nextgr][3].":0";
		$group2bundle{$id}=$bnode[0][$bno];
	    }

	    if(defined $eqposcol{$grcol}){
		my $poscol=$eqposcol{$grcol};
		while(defined $equalcolor{$poscol}) {
		    $poscol=$equalcolor{$poscol};
		}

		my $bno=$bundlecol{$poscol};
		if(defined $bno) { # bundle for group has been created before
		  #print STDERR "Add neutral group ",$currgroup[$nextgr][3]," to pos bundle $bno\n";

		    $bnode[2][$bno]=add_group_to_bundle(\@{$currgroup[$nextgr]},\@{$bundle[2][$bno]},$bnode[2][$bno]);

=begin

		  print STDERR "Current bundle[2][$bno] is: ";
		  my $bn=$bundle[2][$bno][2];
		  while($bn) {
		    print STDERR " ",$$bn[0],"-",$$bn[1]," cov=",$$bn[3];
		    $bn=$$bn[2];
		  }
		  print STDERR "\n";

=cut

		}
		else { # create new bundle
		    $bno=create_bundle(\@bundle,\@nbundle,2,\@{$currgroup[$nextgr]});
		    $bnode[2][$bno]=$bundle[2][$bno][2];

		    #print STDERR "Added neutral group ",$currgroup[$nextgr][3]," to new pos bundle $bno\n";

		    $bundlecol{$poscol}=$bno;
		}
		my $id = $currgroup[$nextgr][3].":2";
		$group2bundle{$id}=$bnode[2][$bno];
	    }	    
	}

	#print STDERR "Done with groups: ",$currgroup[0],",",$currgroup[1],",",$currgroup[2],"\n";


	$currgroup[$nextgr]=$currgroup[$nextgr][5];

    }
		
    #print STDERR "Bundles created!\n";

=begin

    print STDERR "There are ",$nbundle[1]," unstranded bundles ",$nbundle[0]," negative bundles and ",$nbundle[2]," positive bundles\n";    

    for(my $sno=0;$sno<3;$sno++) {
	print STDERR "Bundles on strand $sno:\n";
	for(my $b=0;$b<$nbundle[$sno];$b++) {
	    my $bnode=$bundle[$sno][$b][2];
	    print STDERR "***Bundle $b starting at $bnode:";
	    while($bnode) {
		print STDERR " ",$$bnode[0],"-",$$bnode[1]," cov=",$$bnode[3]/($$bnode[1]-$$bnode[0]+1);
		$bnode=$$bnode[2];
	    }
	    print STDERR "\n";
	}
    }

=cut

    ### predict transcripts for unstranded bundles here
    for(my $b=0;$b<$nbundle[1];$b++) {

      #print STDERR "For $label.$geneno len is ",$bundle[1][$b][0]," coverage is ",$bundle[1][$b][1]," with coords:",$bundle[1][$b][2][0],"-",$bundle[1][$b][2][1],"\n";

	#my $cov=sprintf ("%.6f",$bundle[1][$b][1]/$bundle[1][$b][0]);
	#if($cov >= $readthr ) { # there might be small transfrags that are worth showing

	if($bundle[1][$b][0]>$mintranscriptlen ) { # there might be small transfrags that are worth showing

	    # bundle might contain multiple fragments of a transcript but since we don't know the complete structure -> print only the pieces that are well represented

	    my $bnode=$bundle[1][$b][2];
	    my $t=1;
	    while($bnode) {
		
		#print STDERR "bnode=$bnode at ",$$bnode[1],"-",$$bnode[0],"\n";

		my $cov=$$bnode[3]/($$bnode[1]-$$bnode[0]+1);
		if($cov>=$readthr) {		
		    if($t==1) { $geneno++;}
		    print "$refname\tRlink\ttranscript\t",$$bnode[0],"\t",$$bnode[1],"\t1000\t.\t.\tgene_id \"",$label,".$geneno\"; transcript_id \"",$label,".$geneno.$t\"; cov \"$cov\";\n";
		    print "$refname\tRlink\texon\t",$$bnode[0],"\t",$$bnode[1],"\t1000\t.\t.\tgene_id \"",$label,".$geneno\"; transcript_id \"",$label,".$geneno.$t\"; exon_number \"1\";cov \"$cov\";\n";
		    $t++;
		}
		$bnode=$$bnode[2];
	    }
	}
	#}
    }
    
    #print STDERR "Done with unstranded bundles\n";

    ### build graphs for stranded bundles here
    
    if($startgroup[0] || $startgroup[2]) { # there are stranded groups to process

	# sort junctions
	my @sjunction;
	@{$sjunction[0]}= sort { $a->[1] <=> $b->[1]} @{$junction};
	@{$sjunction[1]}= sort { $a->[2] <=> $b->[2]} @{$junction};
	
	my %bundle2graph;

	my @ngraph=(0,0,0); # one set of graphs for each strand: negative, unknown (ignored), positive
	my @graphno; # how many nodes are in a certain graph
	my @no2gnode;
	
	# build graph structure
	for(my $sno=0;$sno<3;$sno+=2) { # skip neutral bundles -> those shouldn't have junctions
	    for(my $b=0;$b<$nbundle[$sno];$b++) {
		#if($bundle[$sno][$b][1]/$bundle[$sno][$b][0] >= $readthr ) { # bundle is worth processing: it might be that there are small transfrags from source to sink that are worth processing

		if($bundle[$sno][$b][0] >= $mintranscriptlen ) { # bundle is worth processing: it might be that there are small transfrags from source to sink that are worth processing
		  #print STDERR "Process bundle[$sno][$b] to create graph ",$ngraph[$sno]," on strand $sno\n";

		    $graphno[$sno][$ngraph[$sno]]=create_graph($sno,$ngraph[$sno],\@{$bundle[$sno][$b]},\@sjunction,\%bundle2graph,\@no2gnode); # also I need to remember graph coverages somewhere -> probably in the create_graph procedure
		    $ngraph[$sno]++;
		#}
		}
	    }
	}

	#print STDERR "Done creating graphs\n";
	
	my @transfrag; # for each transfrag t on a strand s, in a graph g, transfrag[s][g][t] gives it's abundance
	my @tr2no;  # for each transfrag pattern p on a strand s, in a graph g, tr2no[s][g]{p} gives it's number t
	my @trnumber; # how many transfrags are on a strand s, in a graph g
	my @no2tr; # # for each transfrag t on a strand s, in a graph g, no2tr[s][g][t] gives its pattern p

	# compute probabilities for stranded bundles 
	for(my $n=0;$n<$nreads;$n++) { # this might be a read that was deleted
	    my $np=$$readlist[$n][2];
	    
	    #print STDERR "Consider read $n with pair $np starting at:",$$readlist[$n][3][0][0],"-",$$readlist[$n][3][0][1],"\n";

	    if($np==-1 || $n<$np) {
		my @f = get_fragment_pattern($readlist,$n,$np,\@readgroup,\%merge,\%group2bundle,\%bundle2graph,\@graphno,\@no2gnode);

		my $nf=scalar(@f); # how many fragment patterns were identified (one or two)
	    
		#print STDERR "nf=$nf\n";

		for(my $i=0;$i<$nf;$i++) { # for each pattern
		    for(my $s=0;$s<3;$s+=2) { # for each strand
			if(defined $f[$i][$s][0]) { # this is a valid pattern

			    if(!defined $tr2no[$s][$f[$i][$s][1]]{$f[$i][$s][0]}) {
				$trnumber[$s][$f[$i][$s][1]]++;
				$tr2no[$s][$f[$i][$s][1]]{$f[$i][$s][0]}=$trnumber[$s][$f[$i][$s][1]]-1;
				$no2tr[$s][$f[$i][$s][1]][$tr2no[$s][$f[$i][$s][1]]{$f[$i][$s][0]}]=$f[$i][$s][0];
			    }
			    $transfrag[$s][$f[$i][$s][1]][$tr2no[$s][$f[$i][$s][1]]{$f[$i][$s][0]}]+=1/$$readlist[$n][1];

			}
		    }
		}
	    }
	} # here I am done with read parsing -> they could be deleted if necessary
	
	#print STDERR "Done read parsing and updating frequencies\n";
	    

	my @compatible; # compatibility matrix

	#print STDERR "Transfrags in graphs are:\n";
	for(my $s=0;$s<3;$s+=2) {
	  for(my $g=0;$g<$ngraph[$s];$g++) {

=begin

            # this was the original code to eliminate transfrags that were below threshold -> now I have to coordinate this with the guides
	    my $update=0;
	    my $t=0;
	    while($t<$trnumber[$s][$g]) {
	      my $nt=$t;
	      while($nt<$trnumber[$s][$g] && $transfrag[$s][$g][$nt]<$trthr) { $nt++;}
	      if($nt>$t) {
		$update=1;
		splice @{$transfrag[$s][$g]}, $t, $nt-$t;
		splice @{$no2tr[$s][$g]}, $t, $nt-$t;
		$trnumber[$s][$g]-=$nt-$t;
	      }
	      else {
		$t++;
	      }
	      
	    }
	    
	    if($update) {
	      %{$tr2no[$s][$g]}=();
	      for(my $t=0;$t<$trnumber[$s][$g];$t++) {
		$tr2no[$s][$g]{$no2tr[$s][$g][$t]}=$t;
	      }
	    }

=cut

	    $trnumber[$s][$g]=update_transfrags($s,$trnumber[$s][$g],$transfrag[$s][$g],$no2tr[$s][$g],$tr2no[$s][$g],$no2gnode[$s][$g],$graphno[$s][$g],$guides);

	    my @trnode; # for each transfrag t, gives its node components trnode[t][0], trnode[t][1], ...
	    my @trset;
	    my @trcode;

	    update_nodes($s,$g,\@trnumber,\@no2tr,\@tr2no,$no2gnode[$s][$g],$graphno[$s][$g],\@compatible,\@trnode,\@trset,\@trcode); # compute transfrag compatibility, node abundances

=begin

	    print STDERR "Transfrags and coverages for nodes:\n";

	    for(my $i=1;$i<$graphno[$s][$g];$i++) {
	      my $nn=scalar(@{$trset[$i]});
	      print STDERR "Node $i nn=$nn :";
	      for(my $t=0;$t<$nn;$t++) {
		print STDERR " ",$trset[$i][$t],"(",$transfrag[$s][$g][$trset[$i][$t]],",",$trcode[$i]{$trset[$i][$t]},")";
	      }
	      print STDERR "\n";
	      
	    }

=cut

	    # parse graph and print transcripts
	    $geneno=find_transcripts($graphno[$s][$g],$no2tr[$s][$g],$tr2no[$s][$g],$transfrag[$s][$g],$no2gnode[$s][$g],\@compatible,\@trnode,\@trset,\@trcode,$refname,$geneno,$label,$s,$readthr,$mintranscriptlen);

	  }
	}
    }


    return($geneno);
}

sub find_guide_pat {
  my ($guides,$no2gnode,$gno,$guidenode)=@_;

  my $guidepat="";
  
  my $inode;
  my $i=1;
  while($i<$gno) {
    $inode=$$no2gnode[$i];
    if(intersect($$inode[0],$$inode[1],$$guides[0][0],$$guides[0][1])) { last;}
    $i++;
  }

  if($i<$gno) { # found start node

    push(@{$guidenode},$i);
    vec($guidepat,$i,1)=0b1;

    my $j=0; # guide coordinate I am looking at
    my $n=scalar(@{$guides});
    while($j<$n) {
      if($$guides[$j][1]<$$inode[1]) {
	if($j==$n-1) { return($guidepat);}
	else { @{$guidenode}=(); return("");}
      }
      elsif($$guides[$j][1]==$$inode[1]) {
	if($j==$n-1) { return($guidepat);}
	$j++;
	my $k=0;
	my $nc=scalar(@{$$inode[3]});
	while($k<$nc && $$inode[3][$k][0]!=$$guides[$j][0]) { $k++;}
	if($k==$nc) { @{$guidenode}=(); return("");}
	vec($guidepat,$$inode[3][$k][2],1)=0b1;
	vec($guidepat,$$inode[3][$k][2]*$gno+$i,1)=0b1;
	$i=$$inode[3][$k][2];
	push(@{$guidenode},$i);
	$inode=$$no2gnode[$i];
      }
      else {
	  if($i<$gno-1 && ($$no2gnode[$i+1][0]==$$inode[1]+1)) { ### !!!  maybe here I need to consider bundledist
	    if(scalar($$inode[3]) && $$inode[3][0][2]==$i+1) { ### !!!  maybe here I need to consider bundledist
	      vec($guidepat,$i+1,1)=0b1;
	      vec($guidepat,($i+1)*$gno+$i,1)=0b1;
	      $i++;
	      push(@{$guidenode},$i);
	      $inode=$$no2gnode[$i];
	    }
	    else { @{$guidenode}=(); return("");}
	  }
	  else { # node doesn't continue
	      if($j==$n-1) { return($guidepat);}
	      else { @{$guidenode}=(); return("");}
	  }
      }
    }
  }
  return($guidepat);
}
    
sub update_transfrags {
  my ($strand,$trnumber,$transfrag,$no2tr,$tr2no,$no2gnode,$gno,$guides)=@_;

  # find guides' patterns
  my @guidepat;
  my @guidenode;
  my $gp=0; # number of guide patterns -> I ignore the ones that do not intersect graph or are single nodes
  my $ng=scalar(@{$guides}); # this is the initial number of guides
  for(my $g=0;$g<$ng;$g++) {

    #print STDERR "guides[$g][0]=",$$guides[$g][0]," guides[$g][1]=",$$guides[$g][1][0][0],"-",$$guides[$g][1][-1][1]," vs ",$$no2gnode[1][0],"-",$$no2gnode[$gno-1][1],"\n";

    if($$guides[$g][0]==$strand && intersect($$guides[$g][1][0][0],$$guides[$g][1][-1][1],$$no2gnode[1][0],$$no2gnode[$gno-1][1])) {
      @{$guidenode[$gp]}=();
      $guidepat[$gp]=find_guide_pat($$guides[$g][1],$no2gnode,$gno,\@{$guidenode[$gp]});
      #print STDERR "Guide pattern $gp: ",unpack("b*",$guidepat[$gp]),"\n";
      #print STDERR "tr2no=", $$tr2no{$guidepat[$gp]}," no nodes=",scalar(@{$guidenode[$gp]}),"\n";exit;
      if((!defined $$tr2no{$guidepat[$gp]}) && scalar(@{$guidenode[$gp]})>1) { # guide isn't among transfrags already and has more than one node, otherwise not worth considering
	$gp++;
      }
    }
  }
	    
  # update transfrags, eliminating the ones supporting guide patterns, and low, unsupported transfrags
  my $update=0;
	     
  my %elim;        # keeps all transfrags that need to be eliminated either because they are bellow the threshold or because they've been used in guides
  my @valguide;    # gives the value of each node in a guide
  my @fillguide;   # for each guide tells me if all its nodes were filled
  my @trguide; # for each node in a guide gives me the transfrags
  my %trnode;

  # sort transfrags so that I can eliminate them in order from guides
  #my @sorttrf= sort {unpack("%32b*",$$no2tr[$b]) <=> unpack("%32b*",$$no2tr[$a])} (keys @{$transfrag});
  my @sorttrf= sort {unpack("%32b*",$$no2tr[$b]) <=> unpack("%32b*",$$no2tr[$a])} (values %{$tr2no});
  for(my $t=0;$t<$trnumber;$t++) {
    my $seen=0;
    # first see if transfrag can contribute to guide
    for(my $g=0;$g<$gp;$g++) {
      if(($$no2tr[$sorttrf[$t]] & $guidepat[$g]) eq $$no2tr[$sorttrf[$t]]) { # the transfrag is completely inside guide
	# transfrag needs to have more than one node inside guide to be considered as contributing to guide coverage
	#print STDERR "Transfrag[",$sorttrf[$t],"]=",unpack("b*",$$no2tr[$sorttrf[$t]])," contributes to guide $g\n";
	my $lastnode=-1;
	my $node=-1;
	my $nn=scalar(@{$guidenode[$g]});
	for(my $n=0;$n<$nn;$n++) {
	  if(vec($$no2tr[$sorttrf[$t]],$guidenode[$g][$n],1)) { # transfrag covers node
	    if($lastnode>-1) { # I have seen a node from transfrag before
	      #print STDERR "  Trf[",$sorttrf[$t],"]=",$$transfrag[$sorttrf[$t]]," fill node[$lastnode]",$guidenode[$g][$lastnode],"\n";
	      if(!defined $valguide[$g][$lastnode]) { $fillguide[$g]++;}
	      $valguide[$g][$lastnode]+=$$transfrag[$sorttrf[$t]];
	      push(@{$trguide[$g][$lastnode]},$sorttrf[$t]);
	      push(@{$trnode{$sorttrf[$t]}[$g]},$lastnode);
	    }
	    $lastnode=$n;
	    $node++;
	  }
	}
	if($node) { # I have seen a node from transfrag before
	  if(!defined $valguide[$g][$lastnode]) { $fillguide[$g]++;}
	  #print STDERR "  Trf[",$sorttrf[$t],"]=",$$transfrag[$sorttrf[$t]]," fill node[$lastnode]",$guidenode[$g][$lastnode],"\n";
	  $valguide[$g][$lastnode]+=$$transfrag[$sorttrf[$t]];
	  push(@{$trguide[$g][$lastnode]},$sorttrf[$t]);
	  #print STDERR "Fill trnode of ",$sorttrf[$t]," for guide $g\n";
	  push(@{$trnode{$sorttrf[$t]}[$g]},$lastnode);
	}
      }
    }
	      
    # add transfrag to pool of being eliminated if less than threshold
    if($$transfrag[$sorttrf[$t]]<$trthr) { $elim{$sorttrf[$t]}=1;} # eliminate transfrag if less than threshold
  }

  #print STDERR "gp=$gp\n";

  # adjust transfrag values according to guides
  my %sharedtr; # for a given transfrag tells me the guides that share it
  my @min;
  for(my $g=0;$g<$gp;$g++) {
    $min[$g]=0;
    if($fillguide[$g]==scalar(@{$guidenode[$g]})) { # all nodes in guide are covered by transfrags
      # first find the minimum coverage value for the guide (min of nodes' coverages)
      my %seen=();
      $min[$g]=$valguide[$g][0];
      my $nn=scalar(@{$guidenode[$g]});
      for(my $n=0;$n<$nn;$n++) {
	if($valguide[$g][$n]<$min[$g]) { $min[$g]=$valguide[$g][$n];}
	my $nt=scalar(@{$trguide[$g][$n]}); # all transfrags in node $n of guide $g
	for(my $t=0;$t<$nt;$t++) {
	  if(!$seen{$trguide[$g][$n][$t]}) {
	    $seen{$trguide[$g][$n][$t]}=1;
	    push(@{$sharedtr{$trguide[$g][$n][$t]}},$g);
	  }
	}
      }
    }
  }

  # update coverages for guides according to transfrag proportions
  my @mult;
  for(my $g=0;$g<$gp;$g++) {

    #print STDERR "fillguide[$g]=",$fillguide[$g]," min[$g]=",$min[$g],"\n";

    if($min[$g]) { # guide was kept (all nodes were covered)
      my $cov=$min[$g];
      # first re-estimate guide coverages according to allowed coverages of transfrags
      my $nn=scalar(@{$guidenode[$g]});
      my $minn=-1;
      for(my $n=0;$n<$nn;$n++) { 
	$valguide[$g][$n]=0;
	my $nt=scalar(@{$trguide[$g][$n]});
	for(my $t=0;$t<$nt;$t++) {
	  if(!defined $mult[$g][$trguide[$g][$n][$t]]) {
	    $mult[$g][$trguide[$g][$n][$t]]=1;
	    my $ns=scalar(@{$sharedtr{$trguide[$g][$n][$t]}});
	    if($ns>1) { # transfrag is shared among multiple guides
	      my $sum=0;
	      for(my $gd=0;$gd<$ns;$gd++) {
		#print STDERR "Transfrag ",$trguide[$g][$n][$t]," shared by guide $gd\n";
		$sum+=$min[$sharedtr{$trguide[$g][$n][$t]}[$gd]];
	      }
	      if($sum) { $mult[$g][$trguide[$g][$n][$t]]=$min[$g]/$sum;}
	    }
	    #print STDERR "Transfrag ",$trguide[$g][$n][$t]," has mult[$g]=",$mult[$g][$trguide[$g][$n][$t]],"\n";
	  }
	  $valguide[$g][$n]+=$$transfrag[$trguide[$g][$n][$t]]*$mult[$g][$trguide[$g][$n][$t]];
	}
	if($valguide[$g][$n]<$cov) { 
	  $cov=$valguide[$g][$n];
	  $minn=$n;
	}
	elsif($minn==-1 && $valguide[$g][$n]==$min[$g]) { 
	  $minn=$n;
	}
      }
      
      # cov is the value of the new transfrag defined by the coverage

      #print STDERR "guide has trnumber=$trnumber, cov=$cov minn=$minn\n";

      $$transfrag[$trnumber]=$cov;
      $$no2tr[$trnumber]=$guidepat[$g];
      $$tr2no{$guidepat[$g]}=$trnumber;
      $trnumber++;

      # update -contributing to guide- transfrag coverages 
      my @nodequota;
      for(my $n=0;$n<$nn;$n++) { $nodequota[$n]=$cov;}
      my $n2fill=$nn;
      while($n2fill) {

	#print STDERR "Fill nodequota[$minn]=",$nodequota[$minn]," with valguide=",$valguide[$g][$minn],"\n";

	# I need to use $cov out of nodequota of $minn node 
	my $nt=scalar(@{$trguide[$g][$minn]});
	my $t=0;
	while($t<$nt && $nodequota[$minn]) {
	  my $trcov=$$transfrag[$trguide[$g][$minn][$t]]*$mult[$g][$trguide[$g][$minn][$t]];

	  #print STDERR " Allowed trcov[",$trguide[$g][$minn][$t],"]=$trcov\n";

	  if($trcov>$nodequota[$minn]) {
	    $trcov=$nodequota[$minn];
	  }
	  $$transfrag[$trguide[$g][$minn][$t]]-=$trcov;
	  if($$transfrag[$trguide[$g][$minn][$t]]<$epsilon) { 
	    $$transfrag[$trguide[$g][$minn][$t]]=0;
	  }
	  if(!$$transfrag[$trguide[$g][$minn][$t]]) { $elim{$trguide[$g][$minn][$t]}=1;}
	  
	  #print STDERR "Get trnode of ",$trguide[$g][$minn][$t]," for guide $g\n";

	  my $ntr=scalar(@{$trnode{$trguide[$g][$minn][$t]}[$g]});
	  for(my $n=0;$n<$ntr;$n++) {
	    $nodequota[$trnode{$trguide[$g][$minn][$t]}[$g][$n]]-=$trcov;
	    $valguide[$g][$trnode{$trguide[$g][$minn][$t]}[$g][$n]]-=$trcov;	    
	    if($nodequota[$trnode{$trguide[$g][$minn][$t]}[$g][$n]]<$epsilon) { 
	      $nodequota[$trnode{$trguide[$g][$minn][$t]}[$g][$n]]=0;
	    }
	    if($valguide[$g][$trnode{$trguide[$g][$minn][$t]}[$g][$n]]<$epsilon) { 
	      $valguide[$g][$trnode{$trguide[$g][$minn][$t]}[$g][$n]]=0;
	    }
	    #print STDERR " nodequota[",$trnode{$trguide[$g][$minn][$t]}[$g][$n],"]=",$nodequota[$trnode{$trguide[$g][$minn][$t]}[$g][$n]],"\n";
	  }
	  $t++;
	}
	$minn=-1;
	$n2fill=0;
	for(my $n=0;$n<$nn;$n++) {
	  if($nodequota[$n]) { 
	    $n2fill++;
	    if($minn==-1) { 
	      $minn=$n;
	      $cov=$valguide[$g][$n];
	    }
	    elsif($valguide[$g][$n]<$cov) {
	      $minn=$n;
	      $cov=$valguide[$g][$n];
	    }
	  }
	}
	#print STDERR "New n2fill=$n2fill minn=$minn\n";
      }
    }
  }

  # eliminate transfrags
  my @tr=sort {$b <=> $a} (keys %elim);
  
  my $update=0;
  my $t=0;
  while($t<=$#tr) {
    my $n=1;
    while($tr[$t]+1==$tr[$t]) { 
      $n++;
    }
    splice @{$transfrag},$tr[$t],$n;
    splice @{$no2tr},$tr[$t],$n;
    $trnumber-=$n;
    $update=1;
    $t+=$n;
  }
  
  if($update) {
    %{$tr2no}=();
    for(my $t=0;$t<$trnumber;$t++) {
      $$tr2no{$$no2tr[$t]}=$t;
    }
  }

  return($trnumber);
}


sub find_transcripts {
    my ($gno,$no2tr,$tr2no,$transfrag,$no2gnode,$compatible,$trnode,$trset,$trcode,$refname,$geneno,$label,$strand,$readthr,$mintranscriptlen)=@_;

    my $ne=0;    # how many edges there are in the graph
    my $maxi=1;
    my @intrcov;
    my @outtrcov;
    my @nodecov;
    my @oneperc;

    for(my $i=0;$i<$gno;$i++) {
      my $inode=$$no2gnode[$i];
      $ne+=scalar(@{$$inode[3]});
      
      if($i) { # for all nodes but the source
	    my $inode=$$no2gnode[$i];
	    my $ilen=$$inode[1]-$$inode[0]+1;

	    #print STDERR "i=$i cov is ",$$inode[7]," len=$ilen:",$$inode[0],"-",$$inode[1],"\n";
	    $nodecov[$i]=$$inode[7]/$ilen;
	    if($nodecov[$i]>$nodecov[$maxi]) {
		$maxi=$i;
	    }

	    #print STDERR "Nodecov[$i]=",$nodecov[$i],"\n";

	    my @trsort= sort { $$trcode[$i]{$b} <=> $$trcode[$i]{$a} } @{$$trset[$i]};
	    my @trin;
	    my @trout;

	    my %cov;

	    my $sum=0;

	    my $j=0;
	    while($j<=$#trsort && $$trcode[$i]{$trsort[$j]} >=4) {

	      if($$trcode[$i]{$trsort[$j]} == 8) {
		$sum+=$$transfrag[$trsort[$j]];
		$intrcov[$i]{$trsort[$j]}=$$transfrag[$trsort[$j]];
		$outtrcov[$i]{$trsort[$j]}=$$transfrag[$trsort[$j]];
	      }
	      elsif(($$trcode[$i]{$trsort[$j]} == 7) || ($$trcode[$i]{$trsort[$j]} == 5)){
		$outtrcov[$i]{$trsort[$j]}=$$transfrag[$trsort[$j]];
		$cov{$trsort[$j]}=$$transfrag[$trsort[$j]];
		push(@trout,$trsort[$j]);
	      }
	      elsif(($$trcode[$i]{$trsort[$j]} == 6) || ($$trcode[$i]{$trsort[$j]} == 4)){
		$intrcov[$i]{$trsort[$j]}=$$transfrag[$trsort[$j]];
		push(@trin,$trsort[$j]);
	      }
	      $j++;
	    }

	    for($j=0;$j<=$#trin;$j++) {
	      $sum+=$intrcov[$i]{$trin[$j]};
	      my $need=$intrcov[$i]{$trin[$j]};
	      my $k=0;
	      while($need && $k<=$#trout) {
		if($$compatible[$trin[$j]][$trout[$k]]) {
		  if($need < $cov{$trout[$k]}) {
		    $cov{$trout[$k]}-=$need;
		    $need=0;
		    if($cov{$trout[$k]}<$epsilon) { $cov{$trout[$k]}=0;}
		    if(!$cov{$trout[$k]}) {
		      splice @trout,$k,1;
		    }
		    last;
		  }
		  else {
		    $need-=$cov{$trout[$k]};
		    if($need<$epsilon) { $need=0;}
		    splice @trout,$k,1;
		  }
		}
		else {
		  $k++;
		}
	      }
	      if($need && scalar(@{$$inode[3]})) { # there are not enough nodes exiting which ar compatible with transcript and this is not last node
		
		# this version just continues a transcript so when it is chosen I have at least compatibilities with it
		$outtrcov[$i]{$trin[$j]}=$need;


=begin 

                # this version adds nodes to nowhere
		if($$trcode[$i]{$trin[$j]}==6) { # this transcript is continued
		  $outtrcov[$i]{$trin[$j]}=$need;
		}
		else { # add node to nowhere because I don't know where I need to go
		  my $t=$$tr2no{$$inode[8]}; # transfrag of current node
		  $outtrcov[$i]{$t}=$need;
		}

=cut

	      }
	    }

	    if(scalar(@trout)) { # I have more out transfrags that are not explained by in ones
	      for(my $k=0;$k<=$#trout;$k++) {
		$sum+=$cov{$trout[$k]};

		# this version just continues a transcript so when it is chosen I have at least compatibilities with it
		if($$inode[4][0][2]) { # if node is not source
		  $intrcov[$i]{$trout[$k]}=$cov{$trout[$k]};


=begin

                # this version adds nodes to nowhere
		if($$trcode[$i]{$trout[$k]}==7) { # this transcript is continued
		  $intrcov[$i]{$trout[$k]}=$cov{$trout[$k]};
		}
		else { # add node to nowhere because I don't know where I need to go
		  my $t=$$tr2no{$$inode[8]}; # transfrag of current node
		  $intrcov[$i]{$t}=$cov{$trout[$k]};
		}

=cut
		}

	      }
	    }
	   
	    if($sum) { $oneperc[$i]=1/$sum;}
	    else { $oneperc[$i]=1;}

=begin

	    print STDERR "Transcripts[$i]=";
	    my $nt=scalar(@{$$trset[$i]});
	    for(my $l=0;$l<$nt;$l++) {
	      print STDERR " ",$$trset[$i][$l],":";
	      if($intrcov[$i]{$$trset[$i][$l]}) { print STDERR $intrcov[$i]{$$trset[$i][$l]};}
	      else { print STDERR "-";}
	      print STDERR ",";
	      if($outtrcov[$i]{$$trset[$i][$l]}) { print STDERR $outtrcov[$i]{$$trset[$i][$l]};}
	      else { print STDERR "-";}
	    }
	    print STDERR "\n";

=cut

	  }
    }

    my @sorttrf= sort {unpack("%32b*",$b) <=> unpack("%32b*",$a)} @{$no2tr};
    my %trforder;
    for(my $t=0;$t<=$#sorttrf;$t++) {
	$trforder{$$tr2no{$sorttrf[$t]}}=$t;
    }

    my %seenedge;

    @sorttrf=();
    @{$trcode}=();

    if($nodecov[$maxi]>=$readthr) { $geneno=parse_trf(0,$ne,$maxi,$no2tr,$tr2no,$no2gnode,$compatible,$trnode,$trset,$trcode,\@oneperc,\@nodecov,$transfrag,\@intrcov,\@outtrcov,\%trforder,\%seenedge,$refname,$geneno,1,$label,$strand,$readthr,$mintranscriptlen,$gno);}

    return($geneno);
}


sub max_compon_size {
  my ($set,$result,$compatible,$cov,$hash,$computed)=@_;

  my $maxsize=0;

  my $n=scalar(@{$set});

  for(my $i=0;$i<$n;$i++) {
    my @agree=();
    for(my $j=$i+1;$j<$n;$j++) {
      if($$compatible[$$set[$i]][$$set[$j]]) {
	push(@agree,$$set[$j]);
      }
    }

    my $size=0;
    if($hash) {
      $size=$$cov[$$set[$i]];
    }
    else {
      $size=$$cov{$$set[$i]};
    }

    my @agres=();
    
    if($#agree>-1) {
      my $id=join(".",@agree);
      my $agreesize=$$computed{$id}[0];
      if(!$agreesize) {
	$agreesize=max_compon_size(\@agree,\@agres,$compatible,$cov,$hash,$computed);
	$$computed{$id}[0]=$agreesize;
	@{$$computed{$id}[1]}=@agres;
      }
      else {
	@agres=@{$$computed{$id}[1]};
      }
      $size+=$agreesize;
    }
    if($size>$maxsize) {
      push(@agres,$$set[$i]);
      $maxsize=$size;
      @{$result}=@agres;
    }
  }

  return($maxsize);
}

sub max_component { # needs to be otimized for speed
  my ($set,$result,$compatible,$cov)=@_;

  my $n=scalar(@{$set});

  my @colors;

  push(@{$colors[0]},$$set[0]);

  my $setno=1;
  my $maxset=0;
  my $maxsetsize=$$cov{$$set[0]};

  my %id2col;  
  $id2col{$$set[0]}=0;
  my @col2id;
  $col2id[0]=$$set[0];
  my @size;
  $size[0]=$maxsetsize;
  
  for(my $i=1;$i<$n;$i++) {
    
    #print STDERR "Consider node ",$$set[$i]," setno=$setno\n";
    my $found=0;
    for(my $j=0;$j<$setno;$j++) {
      my @agree=();
      my $agreesize=get_comp($$set[$i],\@{$colors[$j]},$compatible,$cov,\@agree);
      my $na=scalar(@agree);
      if($na) { # found agreement with some subset of colors[$j]
	$found=1;
	if($na==scalar(@{$colors[$j]})) { # element agrees with whole set
	  my $previd=$col2id[$j];
	  delete $id2col{$previd};
	  $previd.=".".$$set[$i];
	  $id2col{$previd}=$j;
	  $col2id[$j]=$previd;
	  push(@{$colors[$j]},$$set[0]);
	  $size[$j]+=$$cov{$$set[$i]};
	  if($size[$j]>$maxsetsize) {
	    $maxsetsize=$size[$j];
	    $maxset=$j;
	  }
	}
	else { # smaller set
	  push(@agree,$$set[$i]);
	  my $id=join(".",@agree);
	  if(!defined $id2col{$id}) { # I have not seen this set before
	    @{$colors[$setno]}=@agree;
	    $id2col{$id}=$setno;
	    $col2id[$setno]=$id;
	    $size[$setno]=$agreesize+$$cov{$$set[$i]};
	    if($size[$setno]>$maxsetsize) {
	      $maxsetsize=$size[$setno];
	      $maxset=$setno;
	    }
	    $setno++;
	  }
	}
      }
    }
    if(!$found) { # didn't find any agreement -> create set with the element by itself
      push(@{$colors[$setno]},$$set[$i]);
      $id2col{$$set[$i]}=$setno;
      $col2id[$setno]=$$set[$i];
      $size[$setno]=$$cov{$$set[$i]};
      if($size[$setno]>$maxsetsize) {
	$maxsetsize=$size[$setno];
	$maxset=$setno;
      }
      $setno++;
    }
  }
  
  $result=\@{$colors[$maxset]};

  return($maxsetsize);
}
	
sub get_comp {
  my ($elem,$set,$compatible,$cov,$agree)=@_;

  my $size=0;
  
  my $n=scalar(@{$set});
  for(my $i=0;$i<$n;$i++) {
    if($$compatible[$$set[$i]][$elem]) { 
      push(@{$agree},$$set[$i]);
      $size+=$$cov{$$set[$i]};
    }
  }

  return($size);
}

sub parse_trf {
    my ($ve,$ne,$maxi,$no2tr,$tr2no,$no2gnode,$compatible,$trnode,$trset,$trcode,$oneperc,$nodecov,$transfrag,$intrcov,$outtrcov,$trforder,$seenedge,$refname,$geneno,$trno,$label,$strand,$readthr,$mintranscriptlen,$gno)=@_;

    my @path=($maxi);

    my @transcripts=();
    my %computed;
    my @tmp;

    max_compon_size(\@{$$trset[$maxi]},\@tmp,$compatible,$transfrag,1);

    my %pathincov;
    my %pathoutcov;
    my %istranscript;
    my $n=scalar(@tmp);
    for(my $i=0;$i<$n;$i++) {
      $istranscript{$tmp[$i]}=1;
      my $nn=scalar(@{$$trnode[$tmp[$i]]});
      for(my $j=0;$j<$nn;$j++) {
	if($$intrcov[$$trnode[$tmp[$i]][$j]]{$tmp[$i]}) { $pathincov{$maxi}+=$$intrcov[$$trnode[$tmp[$i]][$j]]{$tmp[$i]};}
	if($$outtrcov[$$trnode[$tmp[$i]][$j]]{$tmp[$i]}) { $pathoutcov{$maxi}+=$$outtrcov[$$trnode[$tmp[$i]][$j]]{$tmp[$i]};}
      }
      add_tr2vec($tmp[$i],1,0,$trnode,\@transcripts);
    }
    @tmp=();

    my $pathpat='';
    vec($pathpat,$maxi,1)=0b1;

    $ve=back_to_source_path($maxi,\@path,\%pathincov,\%pathoutcov,\$pathpat,\@transcripts,\%istranscript,$transfrag,$intrcov,$outtrcov,$trnode,$trset,$trcode,$compatible,$ve,$seenedge,$no2gnode,$tr2no,$no2tr,$gno,\%computed);

    if($ve>=0) {

      # need to reorder transcripts
      for(my $t=$#transcripts;$t>=0;$t--) {
	add_tr2vec($transcripts[$t],-1,-1,$trnode,\@tmp);
      }
      @transcripts=@tmp;
      @tmp=();

      $ve=fwd_to_sink_path($maxi,\@path,\%pathincov,\%pathoutcov,\$pathpat,\@transcripts,\%istranscript,$transfrag,$intrcov,$outtrcov,$trnode,$trset,$trcode,$compatible,$ve,$seenedge,$no2gnode,$tr2no,$no2tr,$gno,\%computed);    
    
      if($ve>=0) {
	%pathincov=();
	%pathoutcov=();

	#print STDERR "Path pattern=",unpack("b*",$pathpat),"\n";

	# here add transcripts that might be compatible with path
	for(my $i=0;$i<=$#path;$i++) {
	  if($i) {
	    
	    my @tr=keys %{$$intrcov[$path[$i]]};
	    my $nl=scalar(@tr);
	    for(my $t=0;$t<$nl;$t++) {
	      if(!(defined $istranscript{$tr[$t]}) && (($pathpat & $$no2tr[$tr[$t]]) eq $$no2tr[$tr[$t]])) { # transcript is compatible to path
		#print STDERR "Transcript ",$$in[$path[$i]][$t]," is compatible with path\n";
		$istranscript{$tr[$t]}=1;
		push(@transcripts,$tr[$t]);
	      }
	    }
	  }
	  if($i<$#path) {
	    my @tr=keys %{$$outtrcov[$path[$i]]};
	    my $nr=scalar(@tr);
	    for(my $t=0;$t<$nr;$t++) {
	      if(!(defined $istranscript{$tr[$t]}) && (($pathpat & $$no2tr[$tr[$t]]) eq $$no2tr[$tr[$t]])) { # transcript is compatible to path
		#print STDERR "Transcript ",$$out[$path[$i]][$t]," is compatible with path\n";
		$istranscript{$tr[$t]}=1;
		push(@transcripts,$tr[$t]);
	      }
	    }
	  }
	}

=begin

	print STDERR "Compatible transcripts are:";
	for(my $i=0;$i<=$#transcripts;$i++) {
	  print STDERR " ",$transcripts[$i],"(",$$transfrag[$transcripts[$i]],")";
	}
	print STDERR "\n";

=cut


	my @f; # this keeps the relative flux to first node in path 
	my @nodeflux;
	my @incapacity;
	my @outcapacity;
	my $flux=compute_flux(\@f,\@path,\@transcripts,\%istranscript,$intrcov,$outtrcov,$oneperc,$trnode,\@nodeflux,\@incapacity,\@outcapacity);
        %istranscript=();


=begin

    print STDERR "Print transcript: ";
    for(my $i=0;$i<=$#path;$i++) { print STDERR $path[$i],"-";}
    print STDERR ".\n";

=cut

	if($flux>$epsilon) {
	  ($geneno,$trno)=print_transcript($flux,\@path,\@f,\@nodeflux,$nodecov,$no2gnode,$refname,$geneno,$trno,$label,$strand,$readthr,$mintranscriptlen);
    
	  #print STDERR "ve=$ve ne=$ne\n";

	  if($ve<$ne) { # if I didn't parse all edges at least once
	     
	    #print STDERR "Node coverages:\n";
	    for(my $i=1;$i<$gno;$i++) {
	      #print STDERR "Nodecov[$i]=",$$nodecov[$i],"\n";
	      if($$nodecov[$i]>$$nodecov[$maxi]) {
		$maxi=$i;
	      }
	    }
	    
	    #print STDERR "nodecov[$maxi]=",$$nodecov[$maxi]," flux=$flux readthr=$readthr\n";
	     
	    if($$nodecov[$maxi]>=$readthr) { # if I still have nodes that are above coverage threshold
	
	      update_coverages(\@path,\@nodeflux,\@incapacity,\@outcapacity,$trnode,$trset,$intrcov,$outtrcov,$transfrag,\@transcripts,$trforder,$oneperc);

=begin

	      print STDERR "After update new cov of nodes are:\n";
	      for(my $x=1;$x<$gno;$x++) {
		print STDERR "Node[$x]:";
		my $nn=scalar(@{$$trset[$x]});
		for(my $t=0;$t<$nn;$t++) {
		  print STDERR " ",$$trset[$x][$t],"(",$$transfrag[$$trset[$x][$t]],",";
		  if($$intrcov[$x]{$$trset[$x][$t]}) {
		    print STDERR $$intrcov[$x]{$$trset[$x][$t]};
		  }
		  else {
		    print STDERR "-";
		  }
		  print STDERR ",";
		  if($$outtrcov[$x]{$$trset[$x][$t]}) {
		    print STDERR $$outtrcov[$x]{$$trset[$x][$t]};
		  }
		  else {
		    print STDERR "-";
		  }
		  print STDERR ")";
		}
		print STDERR "\n";
	      }

=cut      

	      @path=();
	      @nodeflux=();
	      @transcripts=();
	      @incapacity=();
	      @outcapacity=();
	      $geneno=parse_trf($ve,$ne,$maxi,$no2tr,$tr2no,$no2gnode,$compatible,$trnode,$trset,$trcode,$oneperc,$nodecov,$transfrag,$intrcov,$outtrcov,$trforder,$seenedge,$refname,$geneno,$trno,$label,$strand,$readthr,$mintranscriptlen,$gno);
	    }
	  }
	}
      }
    }
    
    return($geneno);
}

sub update_coverages {
  my ($path,$nodeflux,$incapacity,$outcapacity,$trnode,$trset,$intrcov,$outtrcov,$transfrag,$transcripts,$trforder,$onepercent)=@_;

  my @inquota=@{$nodeflux};
  my @outquota=@{$nodeflux};
  my $nfilled=0;

  $inquota[$$path[0]]=0;
  $outquota[$$path[-1]]=0;

  my $n=scalar(@{$path}); # how many nodes I need to update their quotas

  # these are all the transcripts in the path that are compatible, now I have to see which ones to pick
  my @sortedtrf= sort { $$trforder{$a} <=> $$trforder{$b} } @{$transcripts};

=begin

  print STDERR "Sorted transcripts and orders:";
  for(my $i=0;$i<=$#sortedtrf;$i++) {
    print STDERR " ",$sortedtrf[$i],"(",$$trforder{$sortedtrf[$i]},")";
  }
  print STDERR "\n";

=cut

  my $t=0;
  my $nt=scalar(@sortedtrf);
  while($t<$nt && $nfilled<$n) { # as long as I didn't fill up the node quota
    my $nn=scalar(@{$$trnode[$sortedtrf[$t]]}); # nn is the number of nodes in the transcript

    my $trleft=0; # remembers how much of the transcript I still have from transfrag

    #print STDERR "Update transfrag ",$sortedtrf[$t]," with covg=",$$transfrag[$sortedtrf[$t]],":\n";

    my $spliced=0;
    my $frac=1;

    for(my $i=0;$i<$nn;$i++) {

=begin

      print STDERR " node=",$$trnode[$sortedtrf[$t]][$i]," oneperc=",$$onepercent[$$trnode[$sortedtrf[$t]][$i]]," intrcov=";
      if($$intrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]}) {
	print STDERR $$intrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]};
      }
      else { print STDERR "-";}
      print STDERR " inquota=",$inquota[$$trnode[$sortedtrf[$t]][$i]]," outtrcov=";
      if($$outtrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]}) {
	print STDERR $$outtrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]};
      }
      else { print STDERR "-";}
      print STDERR " outquota=",$outquota[$$trnode[$sortedtrf[$t]][$i]],"\n";

=cut

      if(!$spliced && !$$intrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]} && !$$outtrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]}) { # this is an inner transfrag
	my $in=1;
	my $out=1;
	if($i) { # transfrag enters node
	  if($$incapacity[$$trnode[$sortedtrf[$t]][$i]]) {
	    $in=$$nodeflux[$$trnode[$sortedtrf[$t]][$i]]/$$incapacity[$$trnode[$sortedtrf[$t]][$i]];
	  }
	}
	if($i<$nn-1) { # transfrag exits node
	  if($$outcapacity[$$trnode[$sortedtrf[$t]][$i]]) {
	    $out=$$nodeflux[$$trnode[$sortedtrf[$t]][$i]]/$$outcapacity[$$trnode[$sortedtrf[$t]][$i]];
	  }
	}
	if($in<$out) { if($in<$frac) { $frac=$in;} }
	elsif($out<$frac) { $frac=$out;}
      }

      if($$intrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]}) { # transfrag enters this node; otherwise I can ignore node
	$spliced=1;
	
	if($inquota[$$trnode[$sortedtrf[$t]][$i]]) { # I still have room to fill this node's inquota
	  
	  if($$intrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]}*$$onepercent[$$trnode[$sortedtrf[$t]][$i]]<$inquota[$$trnode[$sortedtrf[$t]][$i]]) {
	    $inquota[$$trnode[$sortedtrf[$t]][$i]]-=$$intrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]}*$$onepercent[$$trnode[$sortedtrf[$t]][$i]];
	    if($inquota[$$trnode[$sortedtrf[$t]][$i]]<$epsilon) { $inquota[$$trnode[$sortedtrf[$t]][$i]]=0;}
	    $$intrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]}=0;

	    # since I deleted all incoming transcript I don't have anything for trleft

	  }
	  else {
	    $$intrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]}-=$inquota[$$trnode[$sortedtrf[$t]][$i]]/$$onepercent[$$trnode[$sortedtrf[$t]][$i]];
	    if($$intrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]}<$epsilon) { 
	      $$intrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]}=0;
	    }

	    $inquota[$$trnode[$sortedtrf[$t]][$i]]=0;
	  }
	  
	  if(!$inquota[$$trnode[$sortedtrf[$t]][$i]] && !$outquota[$$trnode[$sortedtrf[$t]][$i]]) { # first time both quotas were filled
	    $nfilled++;
	  }
	}

	if($$intrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]}>$trleft) {
	  $trleft=$$intrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]};
	}
	elsif(!$$intrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]}) {
	  delete $$intrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]};
	}
      }
      
      if($$outtrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]}) { # transfrag exits this node; otherwise I can ignore node
	$spliced=1;
	if($outquota[$$trnode[$sortedtrf[$t]][$i]]) { # I still have room to fill this node
	  if($$outtrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]}*$$onepercent[$$trnode[$sortedtrf[$t]][$i]]<$outquota[$$trnode[$sortedtrf[$t]][$i]]) {
	    $outquota[$$trnode[$sortedtrf[$t]][$i]]-=$$outtrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]}*$$onepercent[$$trnode[$sortedtrf[$t]][$i]];
	    if($outquota[$$trnode[$sortedtrf[$t]][$i]]<$epsilon) { $outquota[$$trnode[$sortedtrf[$t]][$i]]=0;}
	    $$outtrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]}=0;
	  }
	  else {
	    $$outtrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]}-=$outquota[$$trnode[$sortedtrf[$t]][$i]]/$$onepercent[$$trnode[$sortedtrf[$t]][$i]];
	    if($$outtrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]}<$epsilon) {
	      $$outtrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]}=0;
	    }
	    $outquota[$$trnode[$sortedtrf[$t]][$i]]=0;
	  }
	  if(!$inquota[$$trnode[$sortedtrf[$t]][$i]] && !$outquota[$$trnode[$sortedtrf[$t]][$i]]) { # first time both quotas were filled
	    $nfilled++;
	  }
	}
	if($$outtrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]}>$trleft) {
	  $trleft=$$outtrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]};
	}
	elsif(!$$outtrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]}) {
	  delete $$outtrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]};
	}
      }
    }
    
    if($trleft) {
      $$transfrag[$sortedtrf[$t]]=$trleft;
    }
    else {
      if(!$spliced) { 
	$trleft=$$transfrag[$sortedtrf[$t]]*(1-$frac);
	if($trleft<$epsilon) { $trleft=0;}
	$$transfrag[$sortedtrf[$t]]=$trleft;
      }
      
      if(!$trleft) {
      
	delete $$transfrag[$sortedtrf[$t]];

	# also delete from trset (there are nodes that don't go into trcov so they need to be cleaned up
	for(my $i=0;$i<$nn;$i++) {
	  my $k=0;
	  my $nk=scalar(@{$$trset[$$trnode[$sortedtrf[$t]][$i]]});
	  
	  #print STDERR "Have to delete transfrag ",$sortedtrf[$t]," from trset[",$$trnode[$sortedtrf[$t]][$i],"] ";
	    
	  while($k<$nk) {
	    if($$trset[$$trnode[$sortedtrf[$t]][$i]][$k] == $sortedtrf[$t]) {
	      splice @{$$trset[$$trnode[$sortedtrf[$t]][$i]]},$k,1;
	      last;
	    }
	    $k++;
	  }
	}
      }
    }
    #print STDERR " trleft=$trleft nfilled=$nfilled t=$t nt=$nt\n";

    $t++;
  }

  # check for unspliced fragments
  while($t<$nt) { 
    my $nn=scalar(@{$$trnode[$sortedtrf[$t]]}); # nn is the number of nodes in the transcript
    my $frac=1;
    my $spliced=0;
    for(my $i=0;$i<$nn;$i++) {
      if(!$$intrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]} && !$$outtrcov[$$trnode[$sortedtrf[$t]][$i]]{$sortedtrf[$t]}) { # this is an inner transfrag
	my $in=1;
	my $out=1;
	if($i) { # transfrag enters node
	  if($$incapacity[$$trnode[$sortedtrf[$t]][$i]]) {
	    $in=$$nodeflux[$$trnode[$sortedtrf[$t]][$i]]/$$incapacity[$$trnode[$sortedtrf[$t]][$i]];
	  }
	}
	if($i<$nn-1) { # transfrag exits node
	  if($$outcapacity[$$trnode[$sortedtrf[$t]][$i]]) {
	    $out=$$nodeflux[$$trnode[$sortedtrf[$t]][$i]]/$$outcapacity[$$trnode[$sortedtrf[$t]][$i]];
	  }
	}
	if($in<$out) { if($in<$frac) { $frac=$in;} }
	elsif($out<$frac) { $frac=$out;}
      }
      else {
	$spliced=1;
	last;
      }
    }
  
    if(!$spliced) { 

      my $trleft=$$transfrag[$sortedtrf[$t]]*(1-$frac);
      if($trleft<$epsilon) { $trleft=0;}
      $$transfrag[$sortedtrf[$t]]=$trleft;
    
      if(!$trleft) {
      
	delete $$transfrag[$sortedtrf[$t]];

	# also delete from trset (there are nodes that don't go into trcov so they need to be cleaned up
	for(my $i=0;$i<$nn;$i++) {
	  my $k=0;
	  my $nk=scalar(@{$$trset[$$trnode[$sortedtrf[$t]][$i]]});
	  
	  #print STDERR "Have to delete transfrag ",$sortedtrf[$t]," from trset[",$$trnode[$sortedtrf[$t]][$i],"] ";
	  
	  while($k<$nk) {
	    if($$trset[$$trnode[$sortedtrf[$t]][$i]][$k] == $sortedtrf[$t]) {
	      splice @{$$trset[$$trnode[$sortedtrf[$t]][$i]]},$k,1;
	      last;
	    }
	    $k++;
	  }
	}
      }
    }
    
    $t++;

  }

}

sub print_transcript {
  my ($flux,$path,$f,$nodeflux,$nodecov,$no2gnode,$refname,$geneno,$trno,$label,$strand,$readthr,$mintranscriptlen)=@_;

  my $cov=0;
  my @exon;
  my $len=0;
  my $prevnode=$$no2gnode[0];
  my $nex=0;

  my $n=scalar(@{$path});
  for(my $i=0;$i<$n;$i++) { # this is the transcript; I am going through the path here maybe I can compute node fluxes here too
    
    #print STDERR "flux=$flux f[",$$path[$i],"]=",$$f[$i]," nex=$nex\n";
    
    # $$nodeflux[$$path[$i]]=$flux/$$f[$i];

    my $node=$$no2gnode[$$path[$i]];

    if($$node[0]>$$prevnode[1]+1) { # this is a new exon
      $exon[$nex][0]=$$node[0];
      $exon[$nex][1]=$$node[1];
      $nex++;
    }
    else {
      $exon[$nex-1][1]=$$node[1];
    }
    $len+=$$node[1]-$$node[0];

    $prevnode=$node;

    my $usedcov=$$node[7]*$$nodeflux[$$path[$i]];
    $$nodecov[$$path[$i]]-=$usedcov/($$node[1]-$$node[0]+1);
    $cov+=$usedcov;
  }

  $cov/=$len;

  #print STDERR "cov=$cov len=$len\n";

  if($cov>=$readthr && $len>=$mintranscriptlen) { # print transcript here (I might also choose to keep all transcripts and only print them at the end (depending on the coverages, abundances relative to the other transcripts, or do some elimination technique as I said before
    my $sign="-";
    if($strand) { $sign="+";}
    my $fcov=sprintf ("%.6f",$cov);
    if($trno==1) { $geneno++;}

    print "$refname\tRlink\ttranscript\t",$exon[0][0],"\t",$exon[$nex-1][1],"\t1000\t$sign\t.\tgene_id \"",$label,".$geneno\"; transcript_id \"",$label,".$geneno.$trno\"; cov \"$fcov\";\n";
	
    for(my $n=0;$n<$nex;$n++) {

      print "$refname\tRlink\texon\t",$exon[$n][0],"\t",$exon[$n][1],"\t1000\t$sign\t.\tgene_id \"",$label,".$geneno\"; transcript_id \"",$label,".$geneno.",$trno,"\"; exon_number \"",$n+1,"\";cov \"$fcov\";\n"; # I could put exon coverage here
    }

    $trno++;
  }

  return($geneno,$trno);
}


sub compute_flux {
  my ($f,$path,$transcripts,$istranscript,$intrcov,$outtrcov,$oneperc,$trnode,$nodeflux,$incapacity,$outcapacity)=@_;

  $$incapacity[$$path[0]]=0;
  $$outcapacity[$$path[-1]]=0;

  # compute nodes capacities for all compatible transcripts
  
  my $n=scalar(@{$transcripts});
  for(my $t=0;$t<$n;$t++) {
    
    my $nn=scalar(@{$$trnode[$$transcripts[$t]]});
    for(my $i=0;$i<$nn;$i++) {
      if($$intrcov[$$trnode[$$transcripts[$t]][$i]]{$$transcripts[$t]}) { 
	$$incapacity[$$trnode[$$transcripts[$t]][$i]]+=$$intrcov[$$trnode[$$transcripts[$t]][$i]]{$$transcripts[$t]};
	$$incapacity[$$trnode[$$transcripts[$t]][$i]]=0 unless $$incapacity[$$trnode[$$transcripts[$t]][$i]]>$epsilon;
      }
      if($$outtrcov[$$trnode[$$transcripts[$t]][$i]]{$$transcripts[$t]}) {
	$$outcapacity[$$trnode[$$transcripts[$t]][$i]]+=$$outtrcov[$$trnode[$$transcripts[$t]][$i]]{$$transcripts[$t]};
	$$outcapacity[$$trnode[$$transcripts[$t]][$i]]=0 unless $$outcapacity[$$trnode[$$transcripts[$t]][$i]]>$epsilon;
      }
    }
  }


  $n=scalar(@{$path});

  my $update=1;

  while($update==1) {

    $update=0;

    # check if there are nodes that are not covered
    for(my $i=1;$i<$n-1;$i++) {
      if(!$$incapacity[$$path[$i]] && $$outcapacity[$$path[$i]]) { # transcripts exit node when none enters it
	my @tr=keys %{$$outtrcov[$$path[$i]]};
	my $nr=scalar(@tr);
	for(my $t=0;$t<$nr;$t++) {
	  if($$istranscript{$tr[$t]}) { # need to remove transcript
	    my $nn=scalar(@{$$trnode[$tr[$t]]});
	    for(my $j=0;$j<$nn;$j++) {
	      if($$intrcov[$$trnode[$tr[$t]][$j]]{$tr[$t]}) {
		$$incapacity[$$trnode[$tr[$t]][$j]]-=$$intrcov[$$trnode[$tr[$t]][$j]]{$tr[$t]};
		$$incapacity[$$trnode[$tr[$t]][$j]]=0 unless $$incapacity[$$trnode[$tr[$t]][$j]]>$epsilon;
		if(!$$incapacity[$$trnode[$tr[$t]][$j]] && $$trnode[$tr[$t]][$j]<$$path[$i] && $$trnode[$tr[$t]][$j]>$$path[0]) {
		  $update++;
		}
	      }
	      if($$outtrcov[$$trnode[$tr[$t]][$j]]{$tr[$t]}) {
		$$outcapacity[$$trnode[$tr[$t]][$j]]-=$$outtrcov[$$trnode[$tr[$t]][$j]]{$tr[$t]};
		$$outcapacity[$$trnode[$tr[$t]][$j]]=0 unless $$outcapacity[$$trnode[$tr[$t]][$j]]>$epsilon;
		if(!$$outcapacity[$$trnode[$tr[$t]][$j]] && $$trnode[$tr[$t]][$j]<$$path[$i] && $$trnode[$tr[$t]][$j]>$$path[0]) {
		  $update++;
		}
	      }
	    }
	    $$istranscript{$tr[$t]}=0;
	  }
	}
      }
      
      if(!$$outcapacity[$$path[$i]] && $$incapacity[$$path[$i]]) { # transcripts enters node when none exits it
	my @tr=keys %{$$intrcov[$$path[$i]]};
	my $nl=scalar(@tr);
	for(my $t=0;$t<$nl;$t++) {
	  if($$istranscript{$tr[$t]}) { # need to remove transcript
	    my $nn=scalar(@{$$trnode[$tr[$t]]});
	    for(my $j=0;$j<$nn;$j++) {
	      if($$intrcov[$$trnode[$tr[$t]][$j]]{$tr[$t]}) { 
		$$incapacity[$$trnode[$tr[$t]][$j]]-=$$intrcov[$$trnode[$tr[$t]][$j]]{$tr[$t]};
		$$incapacity[$$trnode[$tr[$t]][$j]]=0 unless $$incapacity[$$trnode[$tr[$t]][$j]]>$epsilon;
		if(!$$incapacity[$$trnode[$tr[$t]][$j]] && $$trnode[$tr[$t]][$j]<$$path[$i] && $$trnode[$tr[$t]][$j]>$$path[0]) {
		  $update++;
		}
	      }
	      if($$outtrcov[$$trnode[$tr[$t]][$j]]{$tr[$t]}) {
		$$outcapacity[$$trnode[$tr[$t]][$j]]-=$$outtrcov[$$trnode[$tr[$t]][$j]]{$tr[$t]};
		$$outcapacity[$$trnode[$tr[$t]][$j]]=0 unless $$outcapacity[$$trnode[$tr[$t]][$j]]>$epsilon;
		if(!$$outcapacity[$$trnode[$tr[$t]][$j]] && $$trnode[$tr[$t]][$j]<$$path[$i] && $$trnode[$tr[$t]][$j]>$$path[0]) {
		  $update++;
		}
	      }
	    }
	    $$istranscript{$tr[$t]}=0;
	  }
	}
      }
    }
  }

  #print STDERR "n=$n\n";

  $$f[0]=1;
  my $flux=1;
  if($n>1) {

    my $i=0;
    while($i<$n && !$$outcapacity[$$path[$i]]) { $$f[$i]=1;$i++;}

    if($i<$n) {
      $$f[$i]=1;
      #print STDERR " before outcapacity[",$$path[$i],"]=",$$outcapacity[$$path[$i]];
      $$outcapacity[$$path[$i]]*=$$oneperc[$$path[$i]]; # the capacity leaving first node

      #print STDERR " outcapacity=",$$outcapacity[$$path[$i]],"\n";
      $flux=$$outcapacity[$$path[$i]];
    }
    else { $flux=0;}

    my $previ=$i;
    $i++;

    #print STDERR "i=$i\n";

    while($i<$n) {

      while($i<$n && !$$incapacity[$$path[$i]]) { $$f[$i]=1; $i++;}
      
      if($i<$n) {
	#print STDERR " before incapacity[",$$path[$i],"]=",$$incapacity[$$path[$i]];
  
	$$incapacity[$$path[$i]]*=$$oneperc[$$path[$i]];
	$$f[$i]=$$outcapacity[$$path[$previ]]*$$f[$previ]/$$incapacity[$$path[$i]];
	#print STDERR "f[$i]=",$$f[$i],"\n";

	if($i<$n-1) {

	  #print STDERR " before outcapacity[",$$path[$i],"]=",$$outcapacity[$$path[$i]];
      
	  $$outcapacity[$$path[$i]]*=$$oneperc[$$path[$i]];
	  my $tmpnodeflux=$$outcapacity[$$path[$i]]*$$f[$i];
	  if($tmpnodeflux<$flux) { $flux=$tmpnodeflux;}
	}
	$previ=$i;
	#print STDERR " incapacity[",$$path[$i],"]=",$$incapacity[$$path[$i]]," outcapacity=",$$outcapacity[$$path[$i]],"\n";
      }

      $i++;
    }
  }

  #print STDERR "flux=$flux\n";

  if($flux) { # compute node coverages
    for(my $i=0;$i<$n;$i++) {

      #print STDERR "compute flux: flux[",$$path[$i],"]=f[$i]=",$$f[$i]," incap=",$$incapacity[$$path[$i]]," outcap=",$$outcapacity[$$path[$i]],"\n";

      if($n==1 || $$incapacity[$$path[$i]] || $$outcapacity[$$path[$i]]) {
	$$nodeflux[$$path[$i]]=$flux/$$f[$i];
      }
      else {
	$$nodeflux[$$path[$i]]=0;
      }
    }
  }

  return($flux);
}


sub fwd_to_sink_path {
  my ($i,$path,$pathincov,$pathoutcov,$pathpat,$transcripts,$istranscript,$transfrag,$intrcov,$outtrcov,$trnode,$trset,$trcode,$compatible,$ve,$seenedge,$no2gnode,$tr2no,$no2tr,$gno,$computed)=@_;

  # find all children -> if no children then go back
  my $inode=$$no2gnode[$i];
    
  #my $nc=scalar(@{$$no2gnode[$i][3]});
  my $nc = scalar(@{$$inode[3]}); # no of parents of i

  if($nc) { # child is present

    my $maxc=0; # max child
    my $maxtrfno=0; # maximum number of transcripts
    my @maxtrf_in; # keeps transcripts to add for maxchild
    my @maxtrf_out; # keeps transcripts to add for maxchild
    my @maxrm; # transcripts that will need to be removed from transcripts

    my $n=scalar(@{$transcripts});

    #print STDERR "Node $i has $nc children\n";

    # choose child that keeps maximum number of transcripts on a *direct* path to child!!!
    for(my $child=0;$child<$nc;$child++) { 
      my %seen; # remembers transcripts that were added before to the $tanscripts set	
      my $c=$$inode[3][$child][2]; 
      my $trfno=0;
      my @trf_in;
      my @trf_out;
      my @rm;

      my $edge=$$inode[3][$child][8];
      vec($edge,$i,1)=0b1;
      vec($edge,$c*$gno+$i,1)=0b1;
      my $te=$$tr2no{$edge};
      
      my %torm;

      my $skip=0;
      # check if transcripts are compatible with child
      my $j=0;
      while($j<$n && $$trnode[$$transcripts[$j]][-1]>=$i) { # transcript could be incompatible to child
	if(!$$compatible[$te][$$transcripts[$j]]) { # edge to child and transcripts are not compatible -> don't consider this child
	  if(can_be_removed($$transcripts[$j],$trnode,$pathincov,$pathoutcov,$intrcov,$outtrcov,$$path[0],$i)) {
	    $trfno-=$$transfrag[$$transcripts[$j]];
	    push(@rm,$j);
	    $torm{$$transcripts[$j]}=1;
	  }
	  else {
	    #print STDERR "Transcript ",$$transcripts[$j]," and edge $te are incompatible\n";
	    $skip=1;
	    last;
	  }
	}
	$seen{$$transcripts[$j]}=1;
	$j++;
      }
      
      if($skip) {
	next; # go to next child
      }
      

      # compute how many of the trset[c] transcripts are compatible with the previous kept transcripts and with the edge to child
      my $nr=scalar(@{$$trset[$c]});

      for(my $t=0;$t<$nr;$t++) {
	if(!$seen{$$trset[$c][$t]}) { # transcript was not added before
	  if($$compatible[$$trset[$c][$t]][$te] && onpath($$trset[$c][$t],$$pathpat,$$path[0],$i,$trnode,$no2gnode,$no2tr,$gno)) {
	    my $skip=0;
	    my $j=0;
	    while($j<$n && $$trnode[$$trset[$c][$t]][0]<=$$trnode[$$transcripts[$j]][-1]) { # transcript could be incompatible to in node
	      if(!(defined $torm{$$transcripts[$j]})) { # transcript wasn't marked to be removed
		if(!$$compatible[$$trset[$c][$t]][$$transcripts[$j]]) { # in transcript and transcripts are not compatible
		  $skip=1;
		  last;
		}
	      }
	      $j++;
	    }
	    if(!$skip) { 
	      if($$intrcov[$c]{$$trset[$c][$t]}) {
		push(@trf_in,$$trset[$c][$t]);
	      }
	      elsif($$outtrcov[$c]{$$trset[$c][$t]}) {
		push(@trf_out,$$trset[$c][$t]);
	      }
	    }
	  }
	}
      }
	  
      my @keeptrf=();
      $trfno+=max_compon_size(\@trf_in,\@keeptrf,$compatible,$$intrcov[$c],0);


=begin

      print STDERR "Result=";
      for(my $r=0;$r<=$#keeptrf;$r++) {
	print STDERR " ",$keeptrf[$r],"(",$$transfrag[$keeptrf[$r]],")";
      }
      print STDERR "\n";

      print STDERR "Child $c has trfno=$trfno\n";

=cut
	  
      if(!$maxc || $trfno>$maxtrfno) { # if I haven't seen any parent so far
	$maxc=$c;
	$maxtrfno=$trfno;
	@maxtrf_in=@keeptrf; # these could be sped up by remembering addresses instead of copying all the array
	@maxtrf_out=@trf_out;
	@maxrm=@rm;
      }
    }

    #print STDERR "Max child is $maxc\n";

    if(!$maxc) { return(-1);}

    push(@{$path},$maxc);
    my $edge=$$path[-2]."-".$$path[-1];
    if(!$$seenedge{$edge}) {
      $ve++;
      $$seenedge{$edge}=1;
    }
    
    vec($$pathpat,$$path[-1],1)=0b1;
    vec($$pathpat,$$path[-1]*$gno+1,1)=0b1;

    # remove @maxrm from transcripts shouldn't I also update coverages??
    for(my $t=$#maxrm;$t>=0;$t--) {
      #print STDERR "Remove transcript ",$$transcripts[$maxrm[$t]],"\n";
      delete $$istranscript{$$transcripts[$maxrm[$t]]};
      splice @{$transcripts},$maxrm[$t],1;

      my $nn=scalar(@{$$trnode[$maxrm[$t]]});
      for(my $i=0;$i<$nn;$i++) {
	if($$intrcov[$$trnode[$maxrm[$t]][$i]]{$maxrm[$t]}) { $$pathincov{$$trnode[$maxrm[$t]][$i]}-=$$intrcov[$$trnode[$maxrm[$t]][$i]]{$maxrm[$t]};}
	if($$outtrcov[$$trnode[$maxrm[$t]][$i]]{$maxrm[$t]}) { $$pathoutcov{$$trnode[$maxrm[$t]][$i]}-=$$outtrcov[$$trnode[$maxrm[$t]][$i]]{$maxrm[$t]};}
      }
    }
    
    # add transcripts from trf to the set of kept transcripts
    for(my $t=0;$t<=$#maxtrf_in;$t++) {
      #print STDERR "Add transcript ",$maxtrf_in[$t]," to transcripts\n";
      add_tr2vec($maxtrf_in[$t],-1,-1,$trnode,$transcripts);
      $$istranscript{$maxtrf_in[$t]}=1;
      my $nn=scalar(@{$$trnode[$maxtrf_in[$t]]});
      for(my $i=0;$i<$nn;$i++) {
	if($$intrcov[$$trnode[$maxtrf_in[$t]][$i]]{$maxtrf_in[$t]}) { $$pathincov{$$trnode[$maxtrf_in[$t]][$i]}+=$$intrcov[$$trnode[$maxtrf_in[$t]][$i]]{$maxtrf_in[$t]};}
	if($$outtrcov[$$trnode[$maxtrf_in[$t]][$i]]{$maxtrf_in[$t]}) { $$pathoutcov{$$trnode[$maxtrf_in[$t]][$i]}+=$$outtrcov[$$trnode[$maxtrf_in[$t]][$i]]{$maxtrf_in[$t]};}
      }
    }


    # add transcript that continue the node
    my $t=0;
    while($t<=$#maxtrf_out) {
      my $u=0;
      while($u<=$#maxtrf_in && $$compatible[$maxtrf_in[$u]][$maxtrf_out[$t]]) { $u++;}
      if($u<=$#maxtrf_in) { # transcript t is incompatible with some transcripts in maxtrf_out
	splice @maxtrf_out,$t,1;
      }
      else {$t++;}
    }
    if(scalar(@maxtrf_out)) {
      my @keeptrf=();
      max_compon_size(\@maxtrf_out,\@keeptrf,$compatible,$$outtrcov[$maxc],0);
      for($t=0; $t<=$#keeptrf;$t++) {
	add_tr2vec($keeptrf[$t],-1,-1,$trnode,$transcripts);
      }
    }


    $ve=fwd_to_sink_path($maxc,$path,$pathincov,$pathoutcov,$pathpat,$transcripts,$istranscript,$transfrag,$intrcov,$outtrcov,$trnode,$trset,$trcode,$compatible,$ve,$seenedge,$no2gnode,$tr2no,$no2tr,$gno,$computed);
  }

  return($ve);
}


sub back_to_source_path {
  my ($i,$path,$pathincov,$pathoutcov,$pathpat,$transcripts,$istranscript,$transfrag,$intrcov,$outtrcov,$trnode,$trset,$trcode,$compatible,$ve,$seenedge,$no2gnode,$tr2no,$no2tr,$gno,$computed)=@_;

  # find all parents -> if parent is source then go back
  my $inode=$$no2gnode[$i];
    
  if($$inode[4][0][2]) { # parent is not source

    my $maxp=0; # max parent
    my $maxtrfno=0; # maximum number of transcripts
    my @maxtrf_out; # keeps transcripts to add for maxparent
    my @maxtrf_in; # keeps transcripts to add for maxparent
    my @maxrm; # transcripts that will need to be removed from transcripts

#    my $maxpcomp=0; # max parent
#    my $maxtrfnocomp=0; # maximum number of transcripts
#    my @maxtrfcomp; # keeps transcripts to add for maxparent

    my $n=scalar(@{$transcripts});
    my $np = scalar(@{$$inode[4]}); # no of parents of i

    #print STDERR "Node $i has $np parents\n";

    # choose parent that keeps maximum number of transcripts on a *direct* path to parent!!!
    for(my $parent=0;$parent<$np;$parent++) { 
      my %seen; # remembers transcripts that were added before to the $tanscripts set	
      my $p=$$inode[4][$parent][2]; # CHECK here that it indeed de-references correctly
      my $trfno=0;
      my @trf_out;
      my @trf_in;
      my @rm;

      my $edge=$$inode[4][$parent][8];
      vec($edge,$i,1)=0b1;
      vec($edge,$i*$gno+$p,1)=0b1;
      my $te=$$tr2no{$edge};
      
      my %torm;

      my $skip=0;
      # check if transcripts are compatible with parent
      my $j=0;
      while($j<$n && $$trnode[$$transcripts[$j]][0]<=$i) { # transcript could be incompatible to parent
	if(!$$compatible[$te][$$transcripts[$j]]) { # edge to parent and transcripts are not compatible -> don't consider this parent
	  if(can_be_removed($$transcripts[$j],$trnode,$pathincov,$pathoutcov,$intrcov,$outtrcov,$i,$$path[-1])) {
	    $trfno-=$$transfrag[$$transcripts[$j]];
	    #print STDERR "add pos $j to rm\n";
	    push(@rm,$j);
	    $torm{$$transcripts[$j]}=1;
	  }
	  else {
	    #print STDERR "Edge $p-$i: $te is incompatible to transcript ",$$transcripts[$j],"\n";
	    $skip=1;
	    last;
	  }
	}
	$seen{$$transcripts[$j]}=1;
	$j++;
      }

      if($skip) {
	next;
      }
      
      # compute how many of the trset[p] transcripts are compatible with the previous kept transcripts and with edge to the parent
      my $nl=scalar(@{$$trset[$p]});

      for(my $t=0;$t<$nl;$t++) {
	if(!$seen{$$trset[$p][$t]}) { # transcript was not added before

	  if($$compatible[$$trset[$p][$t]][$te] && onpath($$trset[$p][$t],$$pathpat,$i,$$path[-1],$trnode,$no2gnode,$no2tr,$gno)) { # transcript needs to be compatible with edge and with path so far
	    my $skip=0;
	    my $j=0;
	    while($j<$n && $$trnode[$$trset[$p][$t]][-1]>=$$trnode[$$transcripts[$j]][0]) { # transcript could be incompatible to in node
	      if(!(defined $torm{$$transcripts[$j]})) { # transcript wasn't marked to be removed
		if(!$$compatible[$$trset[$p][$t]][$$transcripts[$j]]) { # in transcript and transcripts are not compatible
		  $skip=1;
		  last;
		}
	      }
	      $j++;
	    }
	    if(!$skip) { 
	      if($$outtrcov[$p]{$$trset[$p][$t]}) {
		push(@trf_out,$$trset[$p][$t]);
	      }
	      elsif($$intrcov[$p]{$$trset[$p][$t]}) {
		push(@trf_in,$$trset[$p][$t]);
	      }
	    }
	  }
	}
      }
	  
      my @keeptrf=();
      #print STDERR "trf(parent $p)=",join(" ",@trf),"\n";
      $trfno+=max_compon_size(\@trf_out,\@keeptrf,$compatible,$$outtrcov[$p],0);


=begin

      print STDERR "Result=";
      for(my $r=0;$r<=$#keeptrf;$r++) {
	print STDERR " ",$keeptrf[$r],"(",$$transfrag[$keeptrf[$r]],")";
      }
      print STDERR "\n";

      print STDERR "Parent $p has trfno=$trfno\n";

=cut
	  
      if(!$maxp || $trfno>$maxtrfno) { # if I haven't seen any parent so far
	$maxp=$p;
	$maxtrfno=$trfno;
	@maxtrf_out=@keeptrf; # these could be speed up by remembering addresses instead of copying all the array
	@maxtrf_in=@trf_in;
	@maxrm=@rm;
      }
    }

    #print STDERR "Max parent is $maxp\n";

    if(!$maxp) { return(-1);}

    unshift(@{$path},$maxp);
    my $edge=$$path[0]."-".$$path[1];
    if(!$$seenedge{$edge}) {
      $ve++;
      $$seenedge{$edge}=1;
    }
    vec($$pathpat,$$path[0],1)=0b1;
    vec($$pathpat,$$path[1]*$gno+$$path[0],1)=0b1;


    # remove @maxrm from transcripts 
    for(my $t=$#maxrm;$t>=0;$t--) {
      #print STDERR "Remove transcript[",$maxrm[$t],"] ",$$transcripts[$maxrm[$t]]," from transcripts\n";
      delete $$istranscript{$$transcripts[$maxrm[$t]]};
      splice @{$transcripts},$maxrm[$t],1;
      # update path coverages
      my $nn=scalar(@{$$trnode[$maxrm[$t]]});
      for(my $i=0;$i<$nn;$i++) {
	if($$intrcov[$$trnode[$maxrm[$t]][$i]]{$maxrm[$t]}) { $$pathincov{$$trnode[$maxrm[$t]][$i]}-=$$intrcov[$$trnode[$maxrm[$t]][$i]]{$maxrm[$t]};}
	if($$outtrcov[$$trnode[$maxrm[$t]][$i]]{$maxrm[$t]}) { $$pathoutcov{$$trnode[$maxrm[$t]][$i]}-=$$outtrcov[$$trnode[$maxrm[$t]][$i]]{$maxrm[$t]};}
      }
    }


    # add transcripts from trf_out to the set of kept transcripts
    for(my $t=0;$t<=$#maxtrf_out;$t++) {
      #print STDERR "Add transcript ",$maxtrf_out[$t]," to transcripts\n";
      add_tr2vec($maxtrf_out[$t],1,0,$trnode,$transcripts);
      $$istranscript{$maxtrf_out[$t]}=1;
      # update path coverage
      my $nn=scalar(@{$$trnode[$maxtrf_out[$t]]});
      for(my $i=0;$i<$nn;$i++) {
	if($$intrcov[$$trnode[$maxtrf_out[$t]][$i]]{$maxtrf_out[$t]}) { $$pathincov{$$trnode[$maxtrf_out[$t]][$i]}+=$$intrcov[$$trnode[$maxtrf_out[$t]][$i]]{$maxtrf_out[$t]};}
	if($$outtrcov[$$trnode[$maxtrf_out[$t]][$i]]{$maxtrf_out[$t]}) { $$pathoutcov{$$trnode[$maxtrf_out[$t]][$i]}+=$$outtrcov[$$trnode[$maxtrf_out[$t]][$i]]{$maxtrf_out[$t]};}
      }
    }


    # add transcript that continue the node
    my $t=0;
    while($t<=$#maxtrf_in) {
      my $u=0;
      while($u<=$#maxtrf_out && $$compatible[$maxtrf_in[$t]][$maxtrf_out[$u]]) { $u++;}
      if($u<=$#maxtrf_out) { # transcript t is incompatible with some transcripts in maxtrf_out
	splice @maxtrf_in,$t,1;
      }
      else {$t++;}
    }
    if(scalar(@maxtrf_in)) {
      my @keeptrf=();
      max_compon_size(\@maxtrf_in,\@keeptrf,$compatible,$$intrcov[$maxp],0);
      for($t=0; $t<=$#keeptrf;$t++) {
	add_tr2vec($keeptrf[$t],1,0,$trnode,$transcripts);
      }
    }


    $ve=back_to_source_path($maxp,$path,$pathincov,$pathoutcov,$pathpat,$transcripts,$istranscript,$transfrag,$intrcov,$outtrcov,$trnode,$trset,$trcode,$compatible,$ve,$seenedge,$no2gnode,$tr2no,$no2tr,$gno,$computed);
  }

  return($ve);
}

sub onpath {
  my ($t,$pathpat,$mini,$maxi,$trnode,$no2gnode,$no2tr,$gno)=@_;
  
  if($$trnode[$t][0]<$mini) { # mini can be reached through transcript
    if(!vec($$no2gnode[$mini][5],$$trnode[$t][0],1)) { 
      return(0);
    }  
  }

  if($$trnode[$t][-1]>$maxi) { # from maxi I can reach end of transcript
    if(!vec($$no2gnode[$maxi][6],$$trnode[$t][-1],1)) { 
      return(0);
    }  
  }

  my $first=1;

  my $n=scalar(@{$$trnode[$t]});
  for(my $i=0;$i<$n;$i++) {
    if($$trnode[$t][$i]>=$mini && $$trnode[$t][$i]<=$maxi) {
      if(!vec($pathpat,$$trnode[$t][$i],1)) { 
	return(0);
      }
      if($first) {
	$first=0;
	if($i && $$trnode[$t][$i]>$mini && vec($$no2tr[$t],$$trnode[$t][$i]*$gno+$$trnode[$t][$i-1],1)) {
	  return(0);
	}
      }
    }
    if($$trnode[$t][$i]>$maxi) {
      if($i && $$trnode[$t][$i-1]<$maxi && vec($$no2tr[$t],$$trnode[$t][$i]*$gno+$$trnode[$t][$i-1],1)) {
	return(0);
      }
      last;
    }
  }

  return(1);
}


sub can_be_removed {
  my ($t,$trnode,$pathincov,$pathoutcov,$incov,$outcov,$mini,$maxi)=@_;

  my $n=scalar(@{$$trnode[$t]});

  #print STDERR "Check if $t can be removed where mini=$mini maxi=$maxi n=$n\n";

  if($n==1) { return(0);}

  for(my $i=0;$i<$n;$i++) {
    if($$pathincov{$$trnode[$t][$i]}) { 
      if($$trnode[$t][$i]>$mini && $$trnode[$t][$i]<=$maxi) {
	if((defined $$incov[$$trnode[$t][$i]]{$t}) && ($$pathincov{$$trnode[$t][$i]}==$$incov[$$trnode[$t][$i]]{$t})) { 
	  return(0);
	}
      }
    }
    if($$pathoutcov{$$trnode[$t][$i]}) {
      if($$trnode[$t][$i]>=$mini && $$trnode[$t][$i]<$maxi) {
	if((defined $$outcov[$$trnode[$t][$i]]{$t}) && ($$pathoutcov{$$trnode[$t][$i]}==$$outcov[$$trnode[$t][$i]]{$t})) { 
	  return(0);
	}
      }
    }
  }

  return(1);
}

sub find_transfragnodes {
    my ($t,$gno,$trpat,$trnode,$trset,$trcode)=@_;  

=begin

 Codes are assigned as follows:
       node       code
       n          0 
    ..>n..>       1
    ..>n          2
       n..>       3
    -->n          4
       n-->       5
    -->n..>       6
    ..>n-->       7
    -->n-->       8
 where .. means the edge is missing but there is a connection, and -- represents the existence of the edge

=cut
           
    #print STDERR "Find nodes for transfrag $t\n";

    # create trnode; 
    for(my $i=1;$i<$gno;$i++) {
	if(vec($trpat,$i,1)) { # node appears in transfrag

	  #print STDERR "found in node $i\n";
	  push(@{$$trnode[$t]},$i);

	  # create $trset
	  @{$$trset[$i]}=() unless defined $$trset[$i];

	}
    }

    my $n=scalar(@{$$trnode[$t]});
    if($n>1) { # I should eliminate the if if I want all transfrags (including the one node ones) in trset
      for(my $i=0;$i<$n;$i++) {

	# if I need to insert transfrag into node in the order given by the first node :
	# add_tr2vec($t,1,0,$trnode,$$trset[$$trnode[$t][$i]]);

	push(@{$$trset[$$trnode[$t][$i]]},$t); # here only transfrags that have more than one node get inserted

	my $ein=0;
	if($i && vec($trpat,$gno*$$trnode[$t][$i]+$$trnode[$t][$i-1],1)) {
	  $ein=1;
	}

	my $eout=0;
	if($i<$n-1 && vec($trpat,$gno*$$trnode[$t][$i+1]+$$trnode[$t][$i],1)) {
	  $eout=1;
	}
      
	if($i) {
	  if($i<$n-1) { # I can have both in and out edges
	    if($ein && $eout) {
	      $$trcode[$$trnode[$t][$i]]{$t}=8;
	    }
	    elsif($ein) {
	      $$trcode[$$trnode[$t][$i]]{$t}=6;
	    }	    
	    elsif($eout) {
	      $$trcode[$$trnode[$t][$i]]{$t}=7;
	    }	    
	    else {
	      $$trcode[$$trnode[$t][$i]]{$t}=1;
	    }
	  }
	  else { # I am at last node
	    if($ein) { 
	      $$trcode[$$trnode[$t][$i]]{$t}=4;
	    }
	    else {
	      $$trcode[$$trnode[$t][$i]]{$t}=2;
	    }
	  }
	}
	else { # I am at first node
	  if($eout) {
	    $$trcode[$$trnode[$t][$i]]{$t}=5;
	  }
	  else {
	    $$trcode[$$trnode[$t][$i]]{$t}=3;
	  }
	}
      }	 
    }
    else {
      $$trcode[$$trnode[$t][0]]{$t}=0;
    }
  
}

sub add_tr2vec {  # NEED TO CHECK THIS FOR REFERENCING
    my ($t,$dir,$pos,$trnode,$where)=@_;

    my $n=scalar(@{$where});
    while($n>0 && $dir*$$trnode[$t][$pos]<$dir*$$trnode[$$where[$n-1]][$pos]) {
	$$where[$n]=$$where[$n-1];
	$n--;
    }
    $$where[$n]=$t;
}


sub exit_node2tr2 {
    my ($i,$t,$node,$lastnode,$trnode,$transfrag)=@_;

    my $n=scalar(@{$$node[9]});
    while($n>0) {
	my $l=scalar(@{$$trnode[$$node[10]]})-1;
	if($$trnode[$$node[10]][$l]<$lastnode) {
	    $$node[10][$n]=$$node[10][$n-1];
	    $n--;
	}
	else { last;}
    }
    $$node[10][$n]=$t;
    $$node[12]+=$$transfrag[$t];
}
 
sub enter_tr2node {
    my ($i,$t,$node,$firstnode,$trnode,$transfrag)=@_;
   
    my $n=scalar(@{$$node[9]});
    while($n>0 && $$trnode[$$node[9][$n-1]][0]>$firstnode) {
	$$node[9][$n]=$$node[9][$n-1];
	$n--;
    }
    $$node[9][$n]=$t;
    $$node[11]+=$$transfrag[$t];
}

    
sub update_nodes {
    my ($s,$g,$trnumber,$no2tr,$tr2no,$no2gnode,$gno,$compatible,$trnode,$trset,$trcode)=@_; 

    for(my $i=0;$i<$gno;$i++) {

      $$no2gnode[$i][8]='' unless $$no2gnode[$i][8];
      vec($$no2gnode[$i][8],$i,1)=0b1; 
      if(!defined $$tr2no[$s][$g]{$$no2gnode[$i][8]}) { # I have to add the node pattern among the transfrag patterns in order to compute compatibilities
	$$tr2no[$s][$g]{$$no2gnode[$i][8]}=$$trnumber[$s][$g];
	$$trnode[$$trnumber[$s][$g]][0]=$i;
	$$trcode[$i]{$$trnumber[$s][$g]}=0;
	@{$$trset[$i]}=() unless defined $$trset[$i];
	$$no2tr[$s][$g][$$trnumber[$s][$g]]=$$no2gnode[$i][8];
	$$trnumber[$s][$g]++;
      }

      # add the edges to transfrag patterns
      my $nc=scalar(@{$$no2gnode[$i][3]});
      for(my $j=0;$j<$nc;$j++) {
	my $childedge=$$no2gnode[$i][8];
	vec($childedge,$$no2gnode[$i][3][$j][2],1)=0b1;
	vec($childedge,$$no2gnode[$i][3][$j][2]*$gno+$i,1)=0b1;
	if(!defined $$tr2no[$s][$g]{$childedge}) { # I have to add the node pattern among the transfrag patterns in order to compute compatibilities
	  $$tr2no[$s][$g]{$childedge}=$$trnumber[$s][$g];
	  $$trnode[$$trnumber[$s][$g]][0]=$i;
	  $$trnode[$$trnumber[$s][$g]][1]=$$no2gnode[$i][3][$j][2];
	  $$trcode[$i]{$$trnumber[$s][$g]}=5;
	  @{$$trset[$i]}=() unless defined $$trset[$i];
	  $$trcode[$$no2gnode[$i][3][$j][2]]{$$trnumber[$s][$g]}=4;
	  @{$$trset[$$no2gnode[$i][3][$j][2]]}=() unless defined $$trset[$i];
	  $$no2tr[$s][$g][$$trnumber[$s][$g]]=$childedge;
	  $$trnumber[$s][$g]++;
	}
      }

      $$no2gnode[$i][7]=0 unless $$no2gnode[$i][7];
    }

    for(my $t1=0;$t1<$$trnumber[$s][$g];$t1++) {
	# first find nodes associated with transfrag
	if(!defined $$trnode[$t1]) { 
	    find_transfragnodes($t1,$gno,$$no2tr[$s][$g][$t1],$trnode,$trset,$trcode);
	}

	my $n1=scalar(@{$$trnode[$t1]});
	$$compatible[$t1][$t1]=1;

	for(my $t2=$t1+1;$t2<$$trnumber[$s][$g];$t2++) {
	    if(!defined $$trnode[$t2]) { 
		find_transfragnodes($t2,$gno,$$no2tr[$s][$g][$t2],$trnode,$trset,$trcode);
	    }   
	
	    my $n2=scalar(@{$$trnode[$t2]});

	    # here check compatibility between t1 and t2;
	    my $i1=0;
	    my $i2=0;

	    $$compatible[$t1][$t2]=1;
	    $$compatible[$t2][$t1]=1;

	    while($i1<$n1 && $i2<$n2) {

	      if($$trnode[$t1][$i1]==$$trnode[$t2][$i2]) { 
		if($i1==$n1-1 || $i2==$n2-1) { # one transcript finishes -> no need to check anymore
		  $i1=$n1;$i2=$n2;
		}
		else { # advance the smallest one
		  if($$trnode[$t1][$i1+1]<$$trnode[$t2][$i2+1]) {
		    $i1++;
		  }
		  else {
		    $i2++;
		  }
		}
	      }
	      elsif($$trnode[$t1][$i1]<$$trnode[$t2][$i2]) {
		$i1++;
		if(conflict(\$i1,$$trnode[$t2][$i2],\@{$$trnode[$t1]},$n1,$no2gnode,$$no2tr[$s][$g][$t1],$gno,$t1,$t2)) {
		  $$compatible[$t1][$t2]=0;
		  $$compatible[$t2][$t1]=0;
		  last;
		}
	      }
	      else {
		$i2++;
		if(conflict(\$i2,$$trnode[$t1][$i1],\@{$$trnode[$t2]},$n2,$no2gnode,$$no2tr[$s][$g][$t2],$gno,$t1,$t2)) {
		  $$compatible[$t1][$t2]=0;
		  $$compatible[$t2][$t1]=0;
		  last;
		}
	      }
	    }
	  }
      }


}
	

	    
sub conflict {
  my ($i, $node, $trnode,$n,$no2gnode,$trpat,$gno,$t1,$t2)=@_;

  while($$i<$n && $node>$$trnode[$$i]) { $$i++;}

  if($$i<$n && $node==$$trnode[$$i]) { return(0);}

  my $node1=$$trnode[$$i-1];

  if($$i==$n) { # transcript is at the end; $node > $node1
    if(vec($$no2gnode[$node1][6],$node,1)) { # node is among children of node1 -> should be equivalent to testing that node1 is among parents of node
      return(0);
    }
    return(1);
  }
  
  my $node2=$$trnode[$$i]; # node1 < node < node2
  
  if(vec($$no2gnode[$node1][6],$node,1) && vec($$no2gnode[$node][6],$node2,1)) {
    if(vec($trpat,$node2*$gno+$node1,1)) {
      return(1);
    }
    else { return(0);}
  }

  return(1);
}
  
  

sub parse {
    my ($ve,$ne,$e,$i,$j,$edge,$no2gnode,$refname,$geneno,$trno,$label,$strand,$readthr,$mintranscriptlen)=@_; # $ve = number of visited edges; when ve == ne then I can exit the parse graph
    
    if($ve==$ne) { return($geneno);} # done all parsing -> explained all edges; some transcripts would be too low in coverage to be printed -> have to deal with this later

    # choose largest weight edge, and start building transcript from here
    #my ($i,$j)=@{$$e[0]}; # edge i->j has largest coverage 

    my $fmax=$$edge[$i][$j][4]; # $$edge[$i][$j][5] keeps flux from i to j in j; 

    # if coverage of edge is below threshold then there is no point in continuing parsing:
    if($$edge[$i][$j][0]<$readthr) { return($geneno);}

    my @f;
    $f[$j]=1; # relation to f[$j] of the chosen node
    $f[0]=1;  # the source doesn't matter
    my @path=($i,$j);  


    #print STDERR "bef bts fmax=$fmax\n";

    $fmax=back_to_source($i,$j,\@f,\@path,$fmax,$edge,$no2gnode); # -> work with unshift here to pre-append nodes to path

    $fmax=fwd_to_sink($j,\@f,\@path,$fmax,$edge,$no2gnode); # -> work with push here to append nodes to path

    #print STDERR "aft fts fmax=$fmax\n";

    # fmax gives the maximum flux in rapport to $f[$j] that is carried from source to sink

    my $cov=0;
    my @exon;
    my $len=0;
    my $prevnode=$$no2gnode[0];
    my $nex=0;
    for(my $n=1;$n<=$#path;$n++) { # this is the transcript

	if($$edge[$path[$n-1]][$path[$n]][3]) { # edge was not visited before
	    $ve++;
	    $$edge[$path[$n-1]][$path[$n]][3]=0;
	}

	my $node=$$no2gnode[$path[$n]];

	update_edge($path[$n-1],$path[$n],$prevnode,$node,\@f,$fmax,$edge,$e,$ne); # here I need to sort the edges after updating the coverage

	if($$node[0]>$$prevnode[1]+1) { # this is a new exon
	    $exon[$nex][0]=$$node[0];
	    $exon[$nex][1]=$$node[1];
	    $nex++;
	}
	else {
	    $exon[$nex-1][1]=$$node[1];
	}
	$len+=$$node[1]-$$node[0];

	$prevnode=$node;

	$cov+=$$node[7]*$fmax/$f[$path[$n]];
    }

    $cov/=$len;

    #print STDERR "cov=$cov len=$len\n";

    if($cov>=$readthr && $len>=$mintranscriptlen) { # print transcript here (I might also choose to keep all transcripts and only print them at the end (depending on the coverages, abundances relative to the other transcripts, or do some elimination technique as I said before
	my $sign="-";
	if($strand) { $sign="+";}
	my $fcov=sprintf ("%.6f",$cov);
	if($trno==1) { $geneno++;}

	print "$refname\tRlink\ttranscript\t",$exon[0][0],"\t",$exon[$nex-1][1],"\t1000\t$sign\t.\tgene_id \"",$label,".$geneno\"; transcript_id \"",$label,".$geneno.$trno\"; cov \"$fcov\";\n";
	
	for(my $n=0;$n<$nex;$n++) {

	    print "$refname\tRlink\texon\t",$exon[$n][0],"\t",$exon[$n][1],"\t1000\t$sign\t.\tgene_id \"",$label,".$geneno\"; transcript_id \"",$label,".$geneno.",$trno,"\"; exon_number \"",$n+1,"\";cov \"$fcov\";\n"; # I could put exon coverage here
	}

	$trno++;
    }

    #my @se = sort {$$edge[$$b[0]][$$b[1]][0] <=> $$edge[$$a[0]][$$a[1]][0]} @{$e};
    #$e=\@se;

    #print STDERR "After one transcript the updated edges are as follows:\n";

    my $maxi=0;
    my $maxj=0;
    for(my $l=0;$l<$ne;$l++) {
	if($$edge[$$e[$l][0]][$$e[$l][1]][0]>$$edge[$maxi][$maxj][0]) {
	    $maxi=$$e[$l][0];
	    $maxj=$$e[$l][1];
	}

	#print STDERR "edge[",$$e[$l][0],"][",$$e[$l][1],"]=",$$edge[$$e[$l][0]][$$e[$l][1]][0]," ",$$edge[$$e[$l][0]][$$e[$l][1]][4]," ",$$edge[$$e[$l][0]][$$e[$l][1]][5],"\n";
	
    }


    #print STDERR "maxi=$maxi maxj=$maxj\n";


    $geneno=parse($ve,$ne,$e,$maxi,$maxj,$edge,$no2gnode,$refname,$geneno,$trno,$label,$strand,$readthr,$mintranscriptlen);

    return($geneno);
}


sub update_edge {
    my ($i,$j,$nodei,$nodej,$f,$fmax,$edge,$e,$ne)=@_; # here I need to make sure the edges are sorted after updating the coverage

    #print STDERR "Update edge: $i-$j fmax=$fmax f[$i]=",$$f[$i]," f[$j]=",$$f[$j]," before update [4]=",$$edge[$i][$j][4]," [5]=",$$edge[$i][$j][5]," after update [4]=";

    $$edge[$i][$j][4]-=$fmax/$$f[$j]; if($$edge[$i][$j][4]<0) { $$edge[$i][$j][4]=0;} # hopefully this will happen only because of rounding errors
    $$edge[$i][$j][5]-=$fmax/$$f[$i]; if($$edge[$i][$j][5]<0) { $$edge[$i][$j][5]=0;} # hopefully this will happen only because of rounding errors
    $$edge[$i][$j][0]=($$edge[$i][$j][4]*$$nodej[7]+$$edge[$i][$j][5]*$$nodei[7])/($$nodei[1]-$$nodei[0]+$$nodej[1]-$$nodej[0]+2); # nodei[7] should be undef if the i is the source -> should be equivalent to 0


    #print STDERR $$edge[$i][$j][4]," [5]=",$$edge[$i][$j][5],"\n";
    

}   


sub back_to_source {
    my ($i,$j,$f,$path,$fmax,$edge,$no2gnode)=@_;

    if($i) { # $i is not the source

	$$f[$i]=$$edge[$i][$j][1]*$$f[$j]/$$edge[$i][$j][2]; # flux is conserved on the same path -> this should be improved in the future so it is not computed at every step

	my $maxfi=0;  # max flux entering i
	my $maxpi=-1; # parent of i achieving maximum flux
	my $nodei=$$no2gnode[$i];
	my $np = scalar(@{$$nodei[4]}); # no of parents of i
	for(my $p=0;$p<$np;$p++) { 
	    my $parent=$$nodei[4][$p];
	    # choosing based on edge coverage:
	    #if($$edge[$$parent[2]][$i][0]>$maxfi) { # I need to choose the edge with the maximum flux transported
	     # $maxfi=$$edge[$$parent[2]][$i][0];
	    # choosing based on flux proportions:
	    if($$edge[$$parent[2]][$i][4]>$maxfi) { # choose the edge most likely that a transcript would follow
	      $maxfi=$$edge[$$parent[2]][$i][4];
	      $maxpi=$$parent[2];
	    }
	}
	if($maxfi) { # found a path to source
	    if($$edge[$maxpi][$i][4]*$$f[$i]<$fmax) { # this edge can not transport as much flux as the maximum so far
		$fmax=$$edge[$maxpi][$i][4]*$$f[$i];
	    }
	    unshift(@{$path},$maxpi); # add parent to path
	    $fmax=back_to_source($maxpi,$i,$f,$path,$fmax,$edge,$no2gnode);
	}
	else {
	    print STDERR "Couldn't find path to source from node $i!\n";exit;
	}
    }

    return($fmax);
}

sub fwd_to_sink {
    my ($j,$f,$path,$fmax,$edge,$no2gnode)=@_; 

    my $nodej=$$no2gnode[$j];
    my $nc = scalar(@{$$nodej[3]}); # no of children of j
    if($nc) { # $j has children (otherwise is the end of transcript
	my $maxfj=0;  # max flux exiting j
	my $maxcj=-1; # child of j achieving maximum flux
	my $nodej=$$no2gnode[$j];
	for(my $c=0;$c<$nc;$c++) {
	    my $child=$$nodej[3][$c];
	    # choose based on edge coverage:
	    #if($$edge[$j][$$child[2]][0]>$maxfj) { # I need to choose the edge with the maximum flux transported
	     # $maxfj=$$edge[$j][$$child[2]][0];
	    # choose based on the most likely transcript:
	      if($$edge[$j][$$child[2]][5]>$maxfj) { # I need to choose the edge that the transcript will most likely follow
		$maxfj=$$edge[$j][$$child[2]][5];
		$maxcj=$$child[2];
	    }
	}
	if($maxfj) { # found a path to sink
	    if($$edge[$j][$maxcj][5]*$$f[$j]<$fmax) { # this edge can not transport as much flux as the maximum so far
		$fmax=$$edge[$j][$maxcj][5]*$$f[$j];
	    }
	    push(@{$path},$maxcj); # add child to path
	    $$f[$maxcj]=$$edge[$j][$maxcj][2]*$$f[$j]/$$edge[$j][$maxcj][1]; # flux is conserved on the same path
	    $fmax=fwd_to_sink($maxcj,$f,$path,$fmax,$edge,$no2gnode);
	}
	else {
	    print STDERR "Couldn't find path to sink from node $j!\n";exit;
	}
    }
    
    return($fmax);
}

sub compute_probs {
    my ($gno,$no2gnode,$edge,$e)=@_; # $gno = no of nodes in graph; $no2gnode = gets the actual graph node for a certain graph node number

    my $maxi=0;
    my $maxj=0;
    my $ne=0;

    for(my $j=1;$j<$gno;$j++) { # for each node but the source compute probabilities to transition to children and from parents 

	my $jnode=$$no2gnode[$j];
	my $jlen=$$jnode[1]-$$jnode[0]+1;

	#print STDERR "Solve graphnode $j:",$$jnode[0],"-",$$jnode[1]," with coverage ",$$jnode[7]," len=",$$jnode[1]-$$jnode[0]+1,"and parent sum(n1)=",$$jnode[10]," and children sum(n1)=",$$jnode[11],"\n";

	# first compute probabilities to transition from parents (all nodes have parents; but not all nodes have children)
	my $sum=0;
	my $np=scalar(@{$$jnode[4]}); # no of parents for node j
	for(my $i=0;$i<$np;$i++) { 

	  #print STDERR "Check $i th parent (",$$jnode[4][$i][2],") with n1=",$$jnode[8][$i][0]," n2=",$$jnode[8][$i][1]," n3=",$$jnode[8][$i][2],"\n";
	    if($$jnode[8][$i][1]<$$jnode[10]-$$jnode[8][$i][0]) { # this shouldn't be possible because all the n1's of the other nodes are specific
		$$jnode[8][$i][1]=$$jnode[10]-$$jnode[8][$i][0];
	    }
	    my $frac=$$jnode[8][$i][0]/($$jnode[8][$i][0]+$$jnode[8][$i][1]); # fraction of transcripts entering j from i

	    $$jnode[8][$i][0]+=$frac*$$jnode[8][$i][2];
	    $sum+=$$jnode[8][$i][0];

	    #print STDERR "Adjusted n1=",$$jnode[8][$i][0]," n2=",$$jnode[8][$i][1],"\n";

	}
	for(my $i=0;$i<$np;$i++) { 

	    my $inode=$$jnode[4][$i]; # i is the parent
	    my $ilen=$$inode[1]-$$inode[0]+1; # ilen is 0 in case of source

	    @{$$e[$ne]}=($$inode[2],$j);
	    
	    $$edge[$$inode[2]][$j][1]=$$jnode[8][$i][0]/$sum; # estimated fraction of transcripts entering j through edge i->j 
	    $$edge[$$inode[2]][$j][0]+=$$edge[$$inode[2]][$j][1]*$$jnode[7]/($ilen+$jlen); # estimated number of transcripts explained by edge i->j (computed as read/bp coverage for exons i & j 
	    #$$edge[$$inode[2]][$j][0]+=$$edge[$$inode[2]][$j][1]*$$jnode[7]/$jlen; # estimated number of transcripts explained by edge i->j (computed as read/bp coverage for exons i & j 

	    #print STDERR "edge[",$$inode[2],"][$j][0]=",$$edge[$$inode[2]][$j][0],"\n";


	    $$edge[$$inode[2]][$j][3]=1;  # edge not visited before
	    $$edge[$$inode[2]][$j][4]=$$edge[$$inode[2]][$j][1]; # this should be improved in the future!!; probably this keeps how much is retained from the initial fraction??

	    if($$edge[$$inode[2]][$j][0]>$$edge[$maxi][$maxj][0]) {
		$maxi=$$inode[2];
		$maxj=$j;
	    }

	    $ne++;
	}

	# now compute probabilities to transition to children
	$sum=0;
	my $nc=scalar(@{$$jnode[3]}); # no of children for node j
	for(my $k=0;$k<$nc;$k++) {

	  #print STDERR "Check $k th child (",$$jnode[3][$k][2],") with n1=",$$jnode[9][$k][0]," n2=",$$jnode[9][$k][1]," n3=",$$jnode[9][$k][2],"\n";

	    if($$jnode[9][$k][1]<$$jnode[11]-$$jnode[9][$k][0]) { # this shouldn't be possible because all the n1's of the other nodes are specific
		$$jnode[9][$k][1]=$$jnode[11]-$$jnode[9][$k][0];
	    }
	    my $frac=$$jnode[9][$k][0]/($$jnode[9][$k][0]+$$jnode[9][$k][1]);

	    $$jnode[9][$k][0]+=$frac*$$jnode[9][$k][2];
	    $sum+=$$jnode[9][$k][0];
	    #print STDERR "Adjusted n1=",$$jnode[9][$k][0]," n2=",$$jnode[9][$k][1],"\n";

	}
	
	for(my $k=0;$k<$nc;$k++) {
	    my $knode=$$jnode[3][$k];  # k is the child
	    my $klen=$$knode[1]-$$knode[0]+1;

	    $$edge[$j][$$knode[2]][2]=$$jnode[9][$k][0]/$sum; # estimated fraction of transcripts exiting j through edge j->k
	    $$edge[$j][$$knode[2]][5]=$$edge[$j][$$knode[2]][2]; # keeps how much I've to use from the flux in this node
	    $$edge[$j][$$knode[2]][0]+=$$edge[$j][$$knode[2]][2]*$$jnode[7]/($klen+$jlen);
	    #$$edge[$j][$$knode[2]][0]+=$$edge[$j][$$knode[2]][2]*$$jnode[7]/$jlen; 

	    #print STDERR "edge[$j][",$$knode[2],"][0]=",$$edge[$j][$$knode[2]][0],"\n";


	    if($$edge[$j][$$knode[2]][0]>$$edge[$maxi][$maxj][0]) {
		$maxj=$$knode[2];
		$maxi=$j;
	    }

	}

    }

    #my @se= sort {$$edge[$$b[0]][$$b[1]][0] <=> $$edge[$$a[0]][$$a[1]][0]} @e;  # sort edges according to their weight from biggest to smallest

    return($ne,$maxi,$maxj);
}

sub update_graph_freq { # this counts number of reads to compute probabilities; another possibility would be to compute read coverages to compute probabilities; 
    my ($pattern,$nh,$jnodes,$gno,$no2gnode)=@_; # $pattern = fragment pattern; $nh = number of fragment hits; $gno = number of nodes in graph; $no2gnode = gets the actual graph node for a certain graph node number

    my %seen; # some nodes in pattern could be duplicated because read and pair could span the same node in graph

    # for each node in graph that the read spans : consider this i -> j (j is the node covered by the read, and i is the parent)
    #    for(my $j=1;$j<$gno;$j++) { # I don't neet to start at 0 which is the source, because a read is never covering it -> I don't need to go through all nodes but only through the ones that the fragment spans -> this should save me time
    
    my $nj=scalar(@{$jnodes});
    for(my $j=0;$j<$nj;$j++) {

      #print STDERR "Consider node ",$$jnodes[$j],"\n";

	if(!$seen{$$jnodes[$j]} ) {

	    # if(vec($pattern,$j,1)) { # read covers child $j : I don't need to test this anymore since I know it to be true from how the jnodes were created
	    my $jnode=$$no2gnode[$$jnodes[$j]]; # graph node of child j, or parent j
	    my $jpat=""; # jnode pattern
	    vec($jpat,$$jnodes[$j],1)=0b1;

	    # first consider the case of edge i->j where i is the parent of j
	    my $np = scalar(@{$$jnode[4]}); # number of parents of j

	    #print STDERR $$jnodes[$j]," has $np parents\n";

	    if($np==1) { # there is only one parent of j
		$$jnode[8][0][0]+=1/$nh; # n1++ : read is i->j specific; i has index 0 among j's parents
		$$jnode[10]+=1/$nh; # sum all n1++ 

		#print STDERR "One parent: add ",1/$nh," to n1 and sum(n1)\n";

	    }
	    elsif(unpack("%32b*",$pattern & $$jnode[5])) { # more than one parent and fragment spans parents of jnode; unpack("%32b*",$bitvec) efficiently counts the number of bits set in $bitvec

		for(my $i=0;$i<$np;$i++) { # for each parent of jnode
		    my $ij=""; # edge from i to j
		    my $inode=$$jnode[4][$i];
		    my $min=$$inode[2];
		    my $max=$$jnodes[$j];

		    if($min>$max) { $min=$$jnodes[$j];$max=$$inode[2];}

		    vec($ij,$max*$gno+$min,1)=0b1;

		    my $ipat=""; # inode pattern
		    vec($ipat,$$inode[2],1)=0b1;

		    #if(unpack("b*",$edgevec) =~ m/1/) { print STDERR "bit set in string\n";}

		    # establish read counts
		    if(unpack("%32b*",$pattern & $ij)) { # f & ij : read spans edge i->j
			$$jnode[8][$i][0]+=1/$nh; # n1++ : read is i->j specific
			$$jnode[10]+=1/$nh; # sum all n1++ 
			
			#print STDERR "Parent add ",1/$nh," to n1(",$$inode[2],")\n";
		    }
		    elsif(unpack("%32b*", ($pattern & ($pattern ^ ( $ipat | $jpat | $$jnode[6] | $$inode[5]))))) {   # f & non-parents of i = f & non-p(i) = f & non-(i|j|children(j)|parents(i)) = f & (f ^ (i|children(i)|parents(i))) : because x & non-y = x & (x ^ y)
			# here I am in the case that the reads are not i specific
			$$jnode[8][$i][1]+=1/$nh; # n2 ++ : read doesn't go through i->j

			#print STDERR "Parent Add ",1/$nh," to n2(",$$inode[2],")\n";

		    }
		    elsif(unpack("%32b*",$pattern & $ipat)) { # f spans i : f & i
			if(unpack("%32b*",($$jnode[5] & $$inode[6])^ $ij)) { # there is not-ij= path from i to j but not through i->j
			    $$jnode[8][$i][2]+=1/$nh; # n3++ : can't decide if read is i->j specific or not
			    #print STDERR "Parent Add ",1/$nh," to n3(",$$inode[2],")\n";

			}
			else {
			    $$jnode[8][$i][0]+=1/$nh; # n1++ : read is i->j specific
			    $$jnode[10]+=1/$nh; # sum all n1++ 
			    
			    #print STDERR "Parent Add ",1/$nh," to n1(",$$inode[2],")\n";
			}
		    }
		    else { # here I could check if f & parents(i) but then to decide I need to see if there is a path from parents(i) to j not going through i->j which I don't know how to do now
			$$jnode[8][$i][2]+=1/$nh; # n3++ : can't decide if read is i->j specific or not
			#print STDERR "Parent add ",1/$nh," to n3(",$$inode[2],")\n";
		    }
		}

	    } #	else fragment only spans children of j -> i 


	    # now treat the case where jnode is the parent so I consider the edge j->k with k being the child
	    my $nc = scalar(@{$$jnode[3]}); # number of children of j
	    if($nc==1) { # there is only one child of j => the transcripts leaving from j can only go here
		$$jnode[9][0][0]+=1/$nh; # n1++ : read is j->k specific ; k has index 0 among j's children
		$$jnode[11]+=1/$nh; # sum all n1++ 
		
		#print STDERR "One Child Add ",1/$nh," to n1\n";

	    }
	    elsif(unpack("%32b*",$pattern & $$jnode[6])) { # fragment covers some children of j

		for(my $k=0;$k<$nc;$k++) { # for each child of jnode
		    my $jk=""; # edge from j to k
		    my $knode=$$jnode[3][$k];
		    my $max=$$knode[2];
		    my $min=$$jnodes[$j];
		    if($min>$max) { $max=$$jnodes[$j];$min=$$knode[2];}
		    vec($jk,$max*$gno+$min,1)=0b1;

		    my $kpat=""; # knode pattern
		    vec($kpat,$$knode[2],1)=0b1;
		    # establish read counts
		    if(unpack("%32b*",$pattern & $jk)) { # f & jk : read spans edge j->k
			$$jnode[9][$k][0]+=1/$nh; # n1++ : read is j->k specific
			$$jnode[11]+=1/$nh; # sum all n1++ 
			
			#print STDERR "Child Add ",1/$nh," to n1(",$$knode[2],")\n";

		    }
		    elsif(unpack("%32b*", ($pattern & ($pattern ^ ( $kpat | $jpat | $$jnode[5] | $$knode[6]))))) {   # f & non-children of k = f & non-c(k) = f & non-(j|k|children(k)|parents(j)) = f & (f ^ (k|children(k)|parents(k))) : because x & non-y = x & (x ^ y)
			# here I am in the case that the reads are not k specific
			$$jnode[9][$k][1]+=1/$nh; # n2 ++ : read doesn't go through j->k
			
			#print STDERR "Child Add ",1/$nh," to n2(",$$knode[2],")\n";

		    }
		    elsif (unpack("%32b*", $pattern & $kpat)) { # f spans k : f & k
			if(unpack("%32b*",($$knode[5] & $$jnode[6])^ $jk)) { # there is not-jk= path from j to k but not through j->k
			    $$jnode[9][$k][2]+=1/$nh; # n3++ : can't decide if read is i->j specific or not
			    #print STDERR "Child Add ",1/$nh," to n3(",$$knode[2],")\n";
			}
			else {
			    $$jnode[9][$k][0]+=1/$nh; # n1++ : read is j->k specific
			    $$jnode[11]+=1/$nh; # sum all n1++ 
			    
			    #print STDERR "Child add ",1/$nh," to n1(",$$knode[2],")\n";
			}
		    }
		    else { # here I could check if f & children(k) but then to decide I need to see if there is a path from j to children(k) not going through j->k which I don't know how to do now
			$$jnode[9][$k][2]+=1/$nh; # n3++ : can't decide if read is i->j specific or not
			#print STDERR "Child add ",1/$nh," to n3(",$$knode[2],")\n";
		    }
		}

	    } # else fragment only spans parents of j -> k 
	    
	    $seen{$$jnodes[$j]}++;
	}
    }
}
		
		
		

sub get_fragment_pattern {
    my ($read,$n,$np,$readgroup,$merge,$group2bundle,$bundle2graph,$graphno,$no2gnode)=@_;

    my @f; 

    my @firstnode;
    my @r;
    # get read pattern
    if($$read[$n][1]) {
	@r=get_read_pattern($read,$n,$readgroup,$merge,$group2bundle,$bundle2graph,$graphno,$no2gnode,\@firstnode);
    }

    my @p;
    # get pair pattern if pair exists and it hasn't been deleted
    if($np>-1 && $$read[$np][1] ) {
      #print STDERR "Consider pair at $np\n";
	for(my $sno=0;$sno<3;$sno+=2) {
	    if($$read[$n][1] && defined $r[$sno][1]) { # read was valid

	      #print STDERR "Read is valid\n";
		# get pair pattern
		@p=get_read_pattern($read,$np,$readgroup,$merge,$group2bundle,$bundle2graph,$graphno,$no2gnode,\@firstnode);
		
		# check if there are conflicts between patterns and graph -> could this happen if there are no reads covering  a gap in bundle? --> yes
		#die "Read and pair have different graph assingments: $r[$sno][1] vs. $p[$sno][1]!\n" unless $r[$sno][1]==$p[$sno][1];
		my $differentgraphs=0;
		if((!defined $p[$sno][1]) || $r[$sno][1]!=$p[$sno][1]) { $differentgraphs=1;}
		my $gnode;
		my $conflictpattn;
		if(!$differentgraphs) {
		    $gnode=$$no2gnode[$sno][$r[$sno][1]][$firstnode[$sno]];
		    $conflictpattn=$$gnode[5]; # this is the parents pattern of firstnode in pair

		    if($firstnode[$sno]<0) {
		      print STDERR "This isn't supposed to happen!!! n=$n np=$np sno=$sno firstnode[$sno]=",$firstnode[$sno],"\n";
		      print STDERR "r[$sno][1]=",$r[$sno][1]," p[$sno][1]=",$p[$sno][1],"\n";
		      print STDERR "read[$n] at ",$$read[$n][0][0],"-",$$read[$n][0][1];
		      print STDERR " pair[$np] at ",$$read[$np][0][0],"-",$$read[$np][0][1],"\n";
		      print STDERR " gnode=$gnode conflictpattn=",unpack("b*",$conflictpattn),"\n";
		    }

		    vec($conflictpattn,$firstnode[$sno],1)=0b1; ### here I get negative value in offset sometimes; is it possible that firstnode[$sno] is -1 and this is why?
		}
		if($differentgraphs || (($conflictpattn & $r[$sno][0]) ne $r[$sno][0])) { # there is a conflict -> pair parents should contain read pattern
		  #print STDERR "Read $r[$sno][0] is not among parents of $firstnode[$sno]: $conflictpattn!\n";
		    # keep patterns separately
		    $f[0][$sno][0]=$r[$sno][0];
		    $f[0][$sno][1]=$r[$sno][1];
		    #@{$f[0][$sno][2]}=@{$r[$sno][2]};
		    if(defined $p[$sno][1]) {
		      $f[1][$sno][0]=$p[$sno][0];
		      $f[1][$sno][1]=$p[$sno][1];
		      #@{$f[1][$sno][2]}=@{$p[$sno][2]};
		    }
		}
		else { # noconflicts -> patterns can be merged
		    $f[0][$sno][0]=$r[$sno][0] | $p[$sno][0];
		    $f[0][$sno][1]=$r[$sno][1];
		    #@{$f[0][$sno][2]}=@{$r[$sno][2]};
		    #push(@{$f[0][$sno][2]},@{$p[$sno][2]});
		}
	    }
	}
    }
    else { # return only read pattern
	if($$read[$n][1]) {
	    for(my $sno=0;$sno<3;$sno+=2) {
		if(defined $r[$sno][1]) { # read is valid for this strand
		    $f[0][$sno][0]=$r[$sno][0];
		    $f[0][$sno][1]=$r[$sno][1];	
		    #@{$f[0][$sno][2]}=@{$r[$sno][2]};
		}
	    }
	}
    }

    return(@f);
}

sub get_read_pattern {
    my ($read,$n,$readgroup,$merge,$group2bundle,$bundle2graph,$graphno,$no2gnode,$firstnode)=@_;

    my @f; # f[0] is pattern for - strand, f[2] is for + strand; f[s][0] gives the pattern, f[s][1] gives the graphno, $f[s][2] gives the gnodes associated with this read on strand s

    #my @firstnode=(-1,-1,-1);
    @{$firstnode}=(-1,-1,-1);
    my @lastgnode; # $lastgnode[0] is for - strand; [2] is for + strand -> I need these in order to add the edges to the read pattern; check this: if it's not correct than storage was wrong!


    my $ncoord=scalar(@{$$read[$n][3]});

    #print STDERR "Get read pattern for read $n starting at ",$$read[$n][3][0][0],"-",$$read[$n][3][0][1]," with $ncoord coordinates\n";


    my $k=0; # need to keep track of coordinates already added to coverages of graphnodes

    my @valid=(1,1,1);
    my $ngno=scalar(@{$$readgroup[$n]});

    #print STDERR "ngno=$ngno\n";

    for(my $i=0;$i<$ngno;$i++) {
	if($valid[0]+$valid[2]) { # there are still stranded bundles associated with the read
	    my $gr=$$readgroup[$n][$i];

	    #print STDERR "gr=$gr\n";

	    die "Can't identify group number for read $n, coordinates $i!\n" unless defined $gr;
	    while(defined $$merge{$gr}) { $gr=$$merge{$gr};}
	    for(my $sno=0;$sno<3;$sno+=2) {
		if($valid[$sno]) {
		    my $id=$gr.":".$sno;
		    my $bnode=$$group2bundle{$id}; 

		    #print STDERR "id=$id bnode=$bnode:",$$bnode[0],"-",$$bnode[1],"\n";

		    die "Group $id has no bundle node associated!\n" unless !$i || defined($bnode);

		    if(defined($bnode) && defined($$bundle2graph{$bnode})) { # group $id has a bundle node associated with it and bundle was processed
			my $nbnode=scalar(@{$$bundle2graph{$bnode}});
			my $j=0;
			while($j<$nbnode && $k<$ncoord) {
			    my $ngraph=$$bundle2graph{$bnode}[$j][0];
			    my $gnode=$$bundle2graph{$bnode}[$j][1];
			    
			    #print STDERR "j=$j: ngraph=$ngraph gnode=$gnode\n";


			    my $node = $$no2gnode[$sno][$ngraph][$gnode];

			    # compute graphnode coverage by read here (assuming they come in order on the genomic line)
			    #print STDERR "Node=$node:",$$node[0],"-",$$node[1]," k=$k n=$n\n";
			    
			    my $intersect=0;

			    while($k<$ncoord) {
				my $bp = intersect($$node[0],$$node[1],$$read[$n][3][$k][0],$$read[$n][3][$k][1]);
				if($bp) {
				  $intersect=1;
				  #print STDERR "Node and read:",$$read[$n][3][$k][0],"-",$$read[$n][3][$k][1]," intersect with $bp bp\n";
				  $$node[7]+=$bp/$$read[$n][1];
				  
				  if($$read[$n][3][$k][1]<=$$node[1]) { $k++;}
				  else {last;}
				}
				else {last;}
			    }

			    if($intersect) { # read intersects gnode
			      if( defined $lastgnode[$sno]) { # this is not the first time I see a gnode for the read
				#print STDERR "lastgnode=",$lastgnode[$sno],"\n";
				
				die "Conflict with graph number[$sno]: $f[$sno][1] vs. $ngraph!\n" unless $f[$sno][1]==$ngraph;
				# need to update the edge pattern
				my $min=$lastgnode[$sno];
				my $max=$gnode;
				if($min>$gnode) { 
				  $max=$min;
				  $min=$gnode;
				}
				vec($f[$sno][0],$max*$$graphno[$sno][$ngraph]+$min,1)=0b1; # added edge from previous gnode to current one for the read
				#print STDERR "Read pattern is:",unpack("b*",$f[$sno][0]),"\n";
			      }
			      else { # first time considering read
				$f[$sno][1]=$ngraph;	
				$$firstnode[$sno]=$gnode;

=begin			      

				# if gnode is linked to source -> add edge to source here in pattern: this shouldn't matter
				my $parent=$$node[4][0];

				#print STDERR "First time read parent is ",$$parent[2],"\n";

				if($$parent[2]==0) { # parent is source
				  vec($f[$sno][0],$gnode*$$graphno[$sno][$ngraph],1)=0b1;
				  #print STDERR "Read pattern is:",unpack("b*",$f[$sno][0]),"\n";
				}

=cut
			      }
			      $lastgnode[$sno]=$gnode;
			      $f[$sno][0]="" unless $f[$sno][0];
			      vec($f[$sno][0],$gnode,1)=0b1; # here I could remember kids as well to speed things up
			      

			      #print STDERR "Read pattern is:",unpack("b*",$f[$sno][0]),"\n";

			      # push(@{$f[$sno][2]},$gnode); -> not needed in the new implementation
			    }
			    $j++;
			  }
		    }
		    else { # I don't need to process this strand for the read anymore because bundle was not processed or read doesn't belong to that stranded bundle
			$valid[$sno]=0;
			next;
		    }
		}
	    }
	}
    }

    return(@f);
}
    
sub intersect {
    my ($n1,$n2,$r1,$r2)=@_;

    #print STDERR "Check intersect between $n1-$n2 and $r1-$r2\n";

    if($r1>$n2 || $r2<$n1) { return(0);}

    my $beg = $n1 < $r1 ? $r1 : $n1;
    my $end = $n2 < $r2 ? $n2 : $r2;
	
    return($end-$beg+1);
}


sub create_graphnode {
    my ($start,$end,$nodeno,$bundlenode,$bundle2graph,$ngraph)=@_;

    #print STDERR "Create graphnode($nodeno):$start-$end\n";

    if($nodeno) {
      push(@{$$bundle2graph{$bundlenode}},[$ngraph,$nodeno]); # for each bundlenode we remember several graph nodes' nos that we split the bundlenode into 
    }

    my @graphnode=($start,$end,$nodeno); # this is a graph node -> 0: start; 1: end; 2: graph node no; 3: children addresses; 4: parents addresses; 5: parents pattern (all ways to get there); 6: children patterns

    {my @empty=(); $graphnode[3]=\@empty;}
    {my @empty=(); $graphnode[4]=\@empty;}

    return(\@graphnode);
}



sub create_graph {
    my ($sno,$ngraph,$bundle,$sjunction,$bundle2graph,$no2gnode)=@_;
    
    my $source = create_graphnode(0,-1,0); # start of the graph; not real content associated with it; just the children are important to leaving from here
    $$no2gnode[$sno][$ngraph][0]=$source;

    my $njunctions=scalar(@{$$sjunction[0]});
    my $njs=0; # index of sorted junction starts
    my $nje=0; # index of sorted junction ends

    my $graphno=1; # number of nodes in graph

    my %ends; # remembers graph nodes' addresses ending at current position

    my $bundlenode=$$bundle[2]; # bundlenode initialized at start node of the bundle 

    while($bundlenode) {


	my $currentstart=$$bundlenode[0]; # current start is bundlenode's start
	my $endbundle=$$bundlenode[1]; # initialize end with bundlenode's end for now

	#print STDERR "New bundle:$currentstart-$endbundle\n";
	

	my $graphnode=create_graphnode($currentstart,$endbundle,$graphno,$bundlenode,$bundle2graph,$ngraph); # creates a $graphno graphnode  with start at bundle start, and end at bundle end
	$$no2gnode[$sno][$ngraph][$graphno]=$graphnode;
	$graphno++;
	my $end=0;
	while($nje<$njunctions && $$sjunction[1][$nje][2]<=$currentstart) { # read all junction ends at or before the current start -> assuming there are any (at this point, smaller junction ends should not be relevant to this bundle/currentstart

	  #print STDERR "Junction end here: ",$$sjunction[1][$nje][2],"\n";

	    if($$sjunction[1][$nje][2]==$currentstart && $$sjunction[1][$nje][0]+1 == $sno) { # junction ends at current start and is on the same strand and not deleted
		$end=1;
	    }
	    $nje++;
	}

	#print STDERR "end=$end nje=$nje njunctions=$njunctions\n";

	if($end) { # I have nodes finishing here
	    my $node=pop(@{$ends{$currentstart}});
	    #print STDERR "Popped node:",$$node[0],"-",$$node[1],"\n";
	    while(defined $node) { 
		push(@{$$node[3]},$graphnode); # this node is the child of previous node
		push(@{$$graphnode[4]},$node);  # this node has as parent the previous node

		#print STDERR "Parent=",$$node[2]," Child=",$$graphnode[2],"\n";

		$node=pop(@{$ends{$currentstart}});
	    }
	    #print STDERR "Done $end\n";
	}
	else { # this node comes from source directly
	    
	    push(@{$$source[3]},$graphnode); # this node is the child of source
	    push(@{$$graphnode[4]},$source);  # this node has source as parent 

	    #print STDERR "Parent=",$$source[2]," Child=",$$graphnode[2],"\n";

	}

	my $completed=0;
	do {

	  #print STDERR "In do: njunctions=$njunctions nje=$nje njs=$njs\n";

	    while($nje<$njunctions && ($$sjunction[1][$nje][0]+1)!=$sno) { $nje++;}
	    while($njs<$njunctions && ((($$sjunction[0][$njs][0]+1)!=$sno) || ($$sjunction[0][$njs][1]<$currentstart))) { $njs++;} # junctions that start before the current graphnode and I haven't seen them before are part of a different bundle

	    my $minjunction = -1; # process next junction -> either a start or an ending whichever has the first position on the genome
	    #print STDERR "Setting minjunction\n";

	    #print STDERR "current njs junction[$njs](",$$sjunction[0][$njs][0]+1,")=",$$sjunction[0][$njs][1],"-",$$sjunction[0][$njs][2],"\n";
	    #print STDERR "current nje junction[$nje](",$$sjunction[1][$nje][0]+1,")=",$$sjunction[1][$nje][1],"-",$$sjunction[1][$nje][2],"\n";


	    if(($nje<$njunctions && ($$sjunction[1][$nje][2]<$endbundle)) || ($njs<$njunctions && ($$sjunction[0][$njs][1]<=$endbundle))) {	    
		if($nje<$njunctions) { # there are still junctions endings
		    if($njs<$njunctions) { # there are still junctions starting
			$minjunction = $$sjunction[0][$njs][1] > $$sjunction[1][$nje][2] ? 1 : 0;
		    }
		    else {
			$minjunction = 1;
		    }
		}
		else {
		    $minjunction = 0;
		}
	    }

	    #print STDERR "minjunction=$minjunction njs=$njs nje=$nje\n";

	    if($minjunction == 0 ) { # found a start junction here
		$$graphnode[1]=$$sjunction[0][$njs][1]; # set the end of current graphnode to here

		my $pos=$$sjunction[0][$njs][1];

		#print STDERR "minj=0 pos=$pos Adjust end of graphnode(",$$graphnode[2],") to ",$$graphnode[1],"\n";

		while($njs<$njunctions && $$sjunction[0][$njs][1]==$pos ) { # remember ends here

		  #print STDERR "njs=$njs Junction on ",$$sjunction[0][$njs][0]+1,":",$$sjunction[0][$njs][1],"-",$$sjunction[0][$njs][2],"\n";
		    
		    if(($$sjunction[0][$njs][0]+1) == $sno) { 
			push(@{$ends{$$sjunction[0][$njs][2]}},$graphnode);
		    }
		    $njs++;
		}


		if($pos<$endbundle) { # there is still place for another node in this bundle (I might put a limit of length here for the graphnode -> because otherwise one can assume this is just a pre-mRNA fragment)
		  #print STDERR "Create nextnode:",$pos+1,"-",$endbundle,"\n";
		    my $nextnode = create_graphnode($pos+1,$endbundle,$graphno,$bundlenode,$bundle2graph,$ngraph);
		    $$no2gnode[$sno][$ngraph][$graphno]=$nextnode;
		    push(@{$$graphnode[3]},$nextnode); # make nextnode a child of current graphnode
		    push(@{$$nextnode[4]},$graphnode); # make graphnode a parent of nextnode

		    #print STDERR "Parent=",$$graphnode[2]," Child=",$$nextnode[2],"\n";

		    $graphno++;
		    $graphnode=$nextnode;
		}
		else {
		    $completed=1;
		}
	    }
	    elsif($minjunction == 1) { # found a junction end here

	      my $pos=$$sjunction[1][$nje][2];
	      #print STDERR "minj=1\n";
	      while($nje<$njunctions && $$sjunction[1][$nje][2]==$pos) { # read all junction ends at the current start
		$nje++;
	      }
	      
	      if($$graphnode[0]<$pos) { # last created node starts before the position of the new node I want to create

		$$graphnode[1]=$pos-1; # set end of current graphnode here
		
		#print STDERR "pos=$pos Adjust end of graphnode(",$$graphnode[2],") to ",$$graphnode[1],"\n";

		my $nextnode = create_graphnode($pos,$endbundle,$graphno,$bundlenode,$bundle2graph,$ngraph);
		$$no2gnode[$sno][$ngraph][$graphno]=$nextnode;
		push(@{$$graphnode[3]},$nextnode); # make nextnode a child of current graphnode
		push(@{$$nextnode[4]},$graphnode); # make graphnode a parent of nextnode

		#print STDERR "Parent=",$$graphnode[2]," Child=",$$nextnode[2],"\n";

		$graphno++;
		$graphnode=$nextnode;
	      }

	      #print STDERR "Check how many nodes end here at $pos\n";
	      my $node=pop(@{$ends{$pos}});
	      while(defined $node) { 
		push(@{$$node[3]},$graphnode); # this node is the child of previous node
		push(@{$$graphnode[4]},$node);  # this node has as parent the previous node
		
		#print STDERR "Parent=",$$node[2]," Child=",$$graphnode[2],"\n";

		$node=pop(@{$ends{$pos}});
	      }
	      
	    }

	    #print STDERR "njunctions=$njunctions nje=$nje njs=$njs endbundle=$endbundle\n";

	} while(($nje<$njunctions && ($$sjunction[1][$nje][2]<$endbundle)) || ($njs<$njunctions && ($$sjunction[0][$njs][1]<=$endbundle)));
	    

	if(!$completed) {
	    $$graphnode[1]=$endbundle;

	    #print STDERR "Make end of graphnode(",$$graphnode[2],") to ",$$graphnode[1],"\n";

	}

		
	$bundlenode=$$bundlenode[2]; # advance to next bundle
    }

    # finished reading bundle -> now create the parents' and children's patterns
    my %visit;    
    traverse_dfs($source,"",$graphno,\%visit,\@{$$no2gnode[$sno][$ngraph]});

    #print STDERR "source[3]=",$$source[3]," source[5]=",$$source[5],"\n";
    #traverse($source,"",$graphno);

=begin 

    for(my $i=0;$i<$graphno;$i++) {
      my $inode=$$no2gnode[$sno][$ngraph][$i];
      print STDERR "Node $i=",$$inode[2],": parents=",unpack("b*",$$inode[5])," children=",unpack("b*",$$inode[6]),"\n";
    }

=cut

    return($graphno);
}



# shorter? version of traverse
sub traverse_dfs {
    my ($node,$parents,$gno,$visited,$no2gnode)=@_;


    #print STDERR "Traverse node(",$$node[2],"):",$$node[0],"-",$$node[1],"\n";

    $$node[5]='' unless $$node[5];
    $$node[6]='' unless $$node[6];

    $$node[5] = $$node[5] | $parents;

    if($$visited{$$node[2]}) {

	for(my $n=0;$n<$gno;$n++) {

	    if(vec($parents,$n,1)) { # add $$node[6] (children) to all parents of node

		my $pnode=$$no2gnode[$n];
		$$pnode[6]=$$pnode[6] | $$node[6];
	    }
	    elsif(vec($$node[6],$n,1)) { # add $parents to all children of node; a node can not be both parent and child
		my $cnode=$$no2gnode[$n];
		$$cnode[5] = $$cnode[5] | $parents;
	    }
	}
    }
    else {

	$$visited{$$node[2]}=1;

	vec($parents,$$node[2],1)=0b1; # add the node to the parents
	
	my $n = scalar(@{$$node[3]}); # these are all the kids of the node;

	#print STDERR "node has $n kids\n";

	for(my $i=0; $i< $n; $i++) { 
	    my $childnode=$$node[3][$i];
	    
	    my $min;
	    my $max;
	    if($$childnode[2]<$$node[2]) {
		$min=$$childnode[2];
		$max=$$node[2];
	    }
	    else {
		$min=$$node[2];
		$max=$$childnode[2];
	    }

	    my $childparents = $parents;
	    vec($childparents,$max*$gno+$min,1)=0b1; # add edge from node to child to the set of parents from child

	    vec($$node[6],$max*$gno+$min,1)=0b1; # add edge from node to child to the set of node children 

	    $$node[6] = $$node[6] | traverse_dfs($childnode,$childparents,$gno,$visited,$no2gnode);
	
	}
    }

    my $children = $$node[6];
    vec($children,$$node[2],1)=0b1;
    return($children);
}



# longer version of traverse
sub traverse {
    my ($node,$parents,$gno)=@_;

    # CHECK THIS: if node was visited before I still need to add the new parents to all its kids!!

    $$node[5]='' unless $$node[5];
    $$node[6]='' unless $$node[6];

    $$node[5] = $$node[5] | $parents;

    vec($parents,$$node[2],1)=0b1; # add the node to the parents

    my $n = scalar(@{$$node[3]}); # these are all the kids of the node;
    for(my $i=0; $i< $n; $i++) { 
	my $childnode=$$node[3][$i];

	my $min;
	my $max;
	if($$childnode[2]<$$node[2]) {
	    $min=$$childnode[2];
	    $max=$$node[2];
	}
	else {
	    $min=$$node[2];
	    $max=$$childnode[2];
	}

	my $childparents = $parents;
	vec($childparents,$max*$gno+$min,1)=0b1; # add edge from node to child to the set of parents from child

	vec($$node[6],$max*$gno+$min,1)=0b1; # add edge from node to child to the set of node children 

	$$node[6] = $$node[6] | traverse($childnode,$childparents,$gno);
	
    }

    my $children = $$node[6];
    vec($children,$$node[2],1)=0b1;
    return($children);
}


sub create_bundle {
    my ($bundle,$nbundle,$sno,$group)=@_;

    my $bno=$$nbundle[$sno];


    $$nbundle[$sno]++;

    my @bnode; # bundlenodes (bnode) are defined like this: 0:start; 1: end; 2: next bnode in list
    $bnode[0]=$$group[0];
    $bnode[1]=$$group[1];
    $bnode[2]=0;
    $bnode[3]=$$group[4];
    $$bundle[$sno][$bno][0]+=$$group[1]-$$group[0]+1;
    $$bundle[$sno][$bno][1]+=$$group[4];
    $$bundle[$sno][$bno][2]=\@bnode;

    return($bno);
}

sub add_group_to_bundle {
    my ($group,$bundle,$currbnode)=@_;

    #print STDERR "Compare group:",$$group[0],"-",$$group[1]," to currbnode:",$$currbnode[0],"-",$$currbnode[1],"\n";

    if($$group[0]>$$currbnode[1]) { # group after last bnode 
	my @bnode;
	$bnode[0]=$$group[0];
	$bnode[1]=$$group[1];
	$bnode[2]=0;
	$bnode[3]=$$group[4];
	$$currbnode[2]=\@bnode;
	$currbnode=\@bnode;
	$$bundle[0]+=$$group[1]-$$group[0]+1;
	$$bundle[1]+=$$group[4];
    }
    else { # group overlaps bnode
	if($$currbnode[1] < $$group[1]) {
	    $$bundle[0]+= $$group[1] - $$currbnode[1];
	    $$currbnode[1]= $$group[1];
	}
	$$bundle[1]+=$$group[4];
	$$currbnode[3]+=$$group[4];
    }
    return($currbnode);
}
	

sub set_strandcol {
    my ($prevgroup,$group,$grcol,$eqcol,$equalcolor)=@_;

    #print STDERR "Set strand col $grcol for ",$$prevgroup[0],"-",$$prevgroup[1],"with color ",$$prevgroup[2]," and group ",$$group[0],"-",$$group[1],"with color ",$$group[2],"\n";

    my $zerocol=$$eqcol{$$prevgroup[2]};
    if(defined $zerocol) {
      #print STDERR "zerocol=$zerocol\n";
	while(defined $$eqcol{$zerocol}) {
	    $zerocol=$$eqcol{$zerocol};
	    #print STDERR "new zerocol=$zerocol\n";
	}
	$$eqcol{$$prevgroup[2]}=$zerocol;
	
	if($zerocol<$grcol) {
	    $$equalcolor{$grcol}=$zerocol;
	    $$group[2]=$zerocol;
	}
	elsif($grcol<$zerocol) {
	    $$equalcolor{$zerocol}=$grcol;
	    $$eqcol{$$prevgroup[2]}=$grcol;
	}
    }
    else {
	$$eqcol{$$prevgroup[2]}=$grcol;
    }
}


sub get_min_start {
    my ($currgroup) = @_;

    my $nextgr=0;

    if($$currgroup[0]) {
	if($$currgroup[1]) {
	    if($$currgroup[2]) {
		my $twogr = $$currgroup[0][0] < $$currgroup[1][0] ? 0 : 1;
		$nextgr = $$currgroup[$twogr] < $$currgroup[2][0] ? $twogr : 2;
		return($nextgr);
	    }
	    else {
		$nextgr = $$currgroup[0][0] < $$currgroup[1][0] ? 0 : 1;
		return($nextgr);
	    }
	}
	else {
	    if($$currgroup[2]) {
		$nextgr = $$currgroup[0][0] < $$currgroup[2][0] ? 0 : 2;
		return($nextgr);
	    }
	    else {
		return(0);
	    }
	}
    }
    else {
	if($$currgroup[1]) {
	    if($$currgroup[2]) {
		$nextgr = $$currgroup[1][0] < $$currgroup[2][0] ? 1 : 2;
		return($nextgr);
	    }
	    else {
		return(1);
	    }
	}
	else {
	    return(2);
	}
    }
    
    return($nextgr);
}
	


sub add_read_to_group {
    my ($n,$read,$col,$group,$ngroup,$allcurrgroup,$startgroup,$readgroup,$eqcol,$merge)=@_;

    my $sno=$$read[$n][0]+1; # 0: negative strand; 1: zerostrand; 2: positive strand
    my $readcol=$col;

    # check if I've seen read's pair and if yes get its readcol; at the least get read's pair strand if available
    my $np=$$read[$n][2]; # pair read number

    if($np>-1 && $$read[$np][1]) { # read pair exists and it wasn't deleted
	if($np<$n) { # there is a pair and it came before the current read in sorted order of position
	    # first group of pair read is: $$readgroup[$np][0]
	    my $grouppair=$$readgroup[$np][0];
	    while(defined $$merge{$grouppair}) {
		$grouppair=$$merge{$grouppair};
	    }
	    $$readgroup[$np][0]=$grouppair;
	    
	    $readcol=$$group[$$readgroup[$np][0]][2];
	    while(defined $$eqcol{$readcol}) { # get smallest color
		$readcol=$$eqcol{$readcol};
	    }
	    #print STDERR "Adjust color of group ",$$readgroup[$np][0]," to $readcol\n";
	    $$group[$$readgroup[$np][0]][2]=$readcol;
	}
	else { # it's the first time I see the read in the fragment

	    # see if I have the correct read strand
	    my $snop=$$read[$np][0]+1;
	    if($sno!=$snop) { # different strands
		if($sno==1 && $snop!=1) { 
		    $$read[$n][0]=$$read[$np][0];
		    $sno=$snop;
		}
		elsif($snop==1 && $sno!=1) { 
		    $$read[$np][0]=$$read[$n][0];
		}
		else { # conflicting strands -> un-pair reads in the hope that one is right
		    $$read[$n][2]=-1;
		    $$read[$np][2]=-1;
		}
	    }

	    $col++;
	}
    }
    else {
	$col++;
    }

    my $currgroup=$$allcurrgroup[$sno];


    #if($n==10450) { print STDERR "strand=$sno Add read $n starting at ",$$read[$n][3][0][0]," where currgroup=",$currgroup," and col=$col\n";}


    if($currgroup) { # this type of group - negative, unknown, or positive - was created before



	# set currgroup first

      #if($n==10450) { print STDERR "Current group is at coords: ",$$currgroup[0],"-",$$currgroup[1]," with color = ",$$currgroup[2],"\n";}

	my $lastgroup=0; 
	while($currgroup && $$read[$n][3][0][0]>$$currgroup[1]) {
	    $lastgroup=$currgroup;
	    $currgroup=$$currgroup[5];
	}
	if(!$currgroup || $$read[$n][3][0][1]<$$currgroup[0]) { $currgroup=$lastgroup;} # currgroup is null only if we reached end of currgroup list

      #if($n==10450) { print STDERR "After search current group is at coords: ",$$currgroup[0],"-",$$currgroup[1],"\n";}


	# now process each group of coordinates individually
	my $thisgroup=$currgroup;
	my $ncoord=scalar(@{$$read[$n][3]});
	my $lastpushedgroup=-1;
	for(my $i=0;$i<$ncoord;$i++) {

	  #if($n==10450) { print STDERR "Read coord $i: ",$$read[$n][3][$i][0],",",$$read[$n][3][$i][1]," while thisgroup=",$thisgroup," is at coords: ",$$thisgroup[0],"-",$$thisgroup[1],"\n";}

	    # skip groups that are left behind
	    while($thisgroup && $$read[$n][3][$i][0]>$$thisgroup[1]) {
		$lastgroup=$thisgroup;
		$thisgroup=$$thisgroup[5];
	    }

	  #if($n==10450) { print STDERR "thisgroup=",$thisgroup,"\n"; }

	    if($thisgroup && $$read[$n][3][$i][1]>=$$thisgroup[0] ) { # read overlaps group

	      #if($n==10450) { print STDERR "Read overlaps group thisgroup=",$thisgroup," at coords: ",$$thisgroup[0],"-",$$thisgroup[1],"with color ",$$thisgroup[2]," check against coords: ",$$read[$n][3][$i][0],",",$$read[$n][3][$i][1],"\n";}
		
		if($$read[$n][3][$i][0]<$$thisgroup[0]) {
		    $$thisgroup[0]=$$read[$n][3][$i][0];
		}
		
		# find end of new group
		my $nextgroup=$$thisgroup[5];
		while($nextgroup && $$read[$n][3][$i][1]>=$$nextgroup[0]) {
		    merge_fwd_groups($thisgroup,$nextgroup,$merge,$eqcol,1);
		    $nextgroup=$$thisgroup[5];
		}
		if($$read[$n][3][$i][1]>$$thisgroup[1]) {
		    $$thisgroup[1]=$$read[$n][3][$i][1];
		}
		

		# get smallest color of group
		while(defined $$eqcol{$$thisgroup[2]}) {
		    $$thisgroup[2]=$$eqcol{$$thisgroup[2]};
		}
		
	      #print STDERR "New color of this group:",$$thisgroup[0],"-",$$thisgroup[1]," is ",$$thisgroup[2]," where readcol=$readcol\n";
    
		if($readcol!=$$thisgroup[2]) { # read color is different from group color
		    if($readcol<$$thisgroup[2]) { # set group color to current read color
			$$eqcol{$$thisgroup[2]}=$readcol;

			#if($n==10450) { print STDERR "eqcol[",$$thisgroup[2],"]=$readcol\n";}

			$$thisgroup[2]=$readcol;
		    }
		    else { # read color is bigger than group
			$$eqcol{$readcol}=$$thisgroup[2];

			#if($n==10450) { print STDERR "eqcol[$readcol]=",$$thisgroup[2],"\n";}

			$readcol=$$thisgroup[2];
		    }
		}
		
		if($$thisgroup[3] != $lastpushedgroup) {
		    push(@{$$readgroup[$n]},$$thisgroup[3]);
		    $lastpushedgroup=$$thisgroup[3];
		}
		$$thisgroup[4]+=($$read[$n][3][$i][1]-$$read[$n][3][$i][0]+1)/$$read[$n][1];
	    }
	    else { # read is at the end of groups, or read is not overlapping other groups -> $lastgroup should be not null here

	      #if($n==10450) { print STDERR "Read at the end of groups, or not overlapping read: create group at ",$$read[$n][3][$i][0],"-",$$read[$n][3][$i][1]," with color $readcol where lastgroup is $lastgroup\n";}

		$$group[$ngroup][0]=$$read[$n][3][$i][0];
		$$group[$ngroup][1]=$$read[$n][3][$i][1];
		$$group[$ngroup][2]=$readcol;
		$$group[$ngroup][3]=$ngroup;
		$$group[$ngroup][4]=($$read[$n][3][$i][1]-$$read[$n][3][$i][0]+1)/$$read[$n][1]; # defines read coverage here
		$$lastgroup[5]=\@{$$group[$ngroup]};
		$$group[$ngroup][5]=$thisgroup;
		push(@{$$readgroup[$n]},$ngroup);
		$lastpushedgroup=$ngroup;
		$ngroup++;		
		$thisgroup=$lastgroup;
	    }
	}
    }
    else { # create new group of this type

      #if($n==10450) { print STDERR "Create new group of this type with color $readcol";}

	my $ncoord=scalar(@{$$read[$n][3]});
	my $lastgroup=0;
	for(my $i=0;$i<$ncoord;$i++) {

	  #if($n==10450) { print STDERR $$read[$n][3][$i][0],"-",$$read[$n][3][$i][1]," ";}

	    $$group[$ngroup][0]=$$read[$n][3][$i][0];
	    $$group[$ngroup][1]=$$read[$n][3][$i][1];
	    $$group[$ngroup][2]=$readcol;
	    $$group[$ngroup][3]=$ngroup;
	    $$group[$ngroup][4]+=($$read[$n][3][$i][1]-$$read[$n][3][$i][0]+1)/$$read[$n][1]; # defines read coverage here
	    $$group[$ngroup][5]=0;
	    if($lastgroup) { 
		$$lastgroup[5]=\@{$$group[$ngroup]};
	    }
	    else {
		$currgroup=\@{$$group[$ngroup]};
	    }
	    $lastgroup=\@{$$group[$ngroup]};

	    push(@{$$readgroup[$n]},$ngroup);
	    $ngroup++;	
	}
	
      #if($n==10450) { print STDERR "\n";}
    }

    #if($n==10450) { print STDERR "Updated current group is at coords: ",$$currgroup[0],"-",$$currgroup[1]," with color ",$$currgroup[2],"\n";}

    $$allcurrgroup[$sno]=$currgroup;
    
    if(!$$startgroup[$sno]) { $$startgroup[$sno]=$currgroup;}
    
    return($col,$ngroup);
}

sub merge_fwd_groups {
    my ($group1,$group2,$merge,$eqcol,$del)=@_;


    #print STDERR "Merge groups: ",$$group1[0],"-",$$group1[1]," col=",$$group1[2]," and",$$group2[0],"-",$$group2[1]," col=",$$group2[2],"\n";

    # get end of group (group1 is assumed to come before group2)
    $$group1[1]=$$group2[1];
    
    # get smallest color of group
    while(defined $$eqcol{$$group1[2]}) {
	$$group1[2]=$$eqcol{$$group1[2]};
    }
    while(defined $$eqcol{$$group2[2]}) {
	$$group2[2]=$$eqcol{$$group2[2]};
    }

    if($$group1[2]<$$group2[2]) { 
	$$eqcol{$$group2[2]}=$$group1[2];
	#print STDERR "in merge set eqcol[",$$group2[2],"]=",$$eqcol{$$group1[2]},"\n";
    }
    elsif($$group1[2]>$$group2[2]) {
	$$eqcol{$$group1[2]}=$$group2[2];
	#print STDERR "in merge set eqcol[",$$group1[2],"]=",$$group2[2],"\n";
	$$group1[2]=$$group2[2];
    }

    $$group1[4]+=$$group2[4];
    $$group1[5]=$$group2[5];

    $$merge{$$group2[3]}=$$group1[3];

    if($del) { @{$group2}=();}
}


sub clean_junctions {
    my ($junction,$junctionthr)=@_;

    my $njunctions=scalar(@{$junction});

    my $n=$njunctions;
    #print STDERR "njunctions=$njunctions\n";
    
    for(my $i=0;$i<$njunctions;$i++) {

      #print STDERR "Consider junction: ",$$junction[$i][1],"-",$$junction[$i][2]," with support: ",$$junction[$i][4],"\n";

	if($$junction[$i][4]<$junctionthr) {
	    $$junction[$i][0]=0;
	    $n--;

	    #print STDERR "Deleted junction $i\n";

	}
    }
    #print STDERR "junctions after cleaning: $n\n";
}

sub process_read {
    my ($currentend,$readlist,$hashread,$junction,$hashjunction,$junctionsupport,$readname,$flag,$readstart,$cigar,$pairstart,$strand,$nh,$hi)=@_;

    my $nreads=scalar(@{$readlist});
    my $n=$nreads;

    my $id=$readname.":".$readstart.":".$hi;
    $$hashread{$id}=$n;

    $$readlist[$nreads++][0]=$strand;
    $$readlist[$n][1]=$nh;
    $$readlist[$n][2]=-1;
    
    # check if I've seen it's pair 
    my $pairid=$readname.":".$pairstart.":".$hi;    
    if(defined $$hashread{$pairid}) {

	my $np=$$hashread{$pairid};

	$$readlist[$n][2]=$np;
	# add read to pair
	$$readlist[$np][2]=$n;

	# check if the reads have nonconflicting strands --> better keep this until working with junction because some of the spliced reads might be wrong!

=begin

	if($$readlist[$n][0]!=$$readlist[$np][0]) { # different strands
	    if($$readlist[$n][0]==0) { # read is unstranded
		$$readlist[$n][0]=$$readlist[$np][0];
	    }
	    elsif($$readlist[$np][0]==0) { # pair is unstranded
		$$readlist[$np][0]=$$readlist[$n][0];
	    }
	    else { # conflicting reads -> delete both of them? or unpair them in the chance that one is right? (chose this for now)
		$$readlist[$n][2]=-1;
		$$readlist[$np][2]=-1;
	    }
	}

=cut

    }

	

    # process cigar string

    my $start=$readstart;
    my $end=$start-1;
    my @m=($cigar =~ /(\d+)([A-Z])/g);
    my $leftsupport=0;
    my $support=0;
    my $rightstart=0;
    my $leftend=0;
    my $maxleftsupport=0;
    my @newjunction;
    my @leftsup;
    my @rightsup;
    my $njunc=0;

    for(my $j=0;$j<=$#m;$j+=2) {
	if($m[$j+1] eq 'M') {
	    $end+=$m[$j];
	    $support+=$m[$j];
	}
	elsif($m[$j+1] eq 'D') {
	  $end+=$m[$j];


=begin

	  if($start<=$end) { push(@{$$readlist[$n][3]},[$start,$end]);}
	  $start=$end+$m[$j]+1;
	  $end=$start-1;

=cut

	}
	elsif($m[$j+1] eq 'N') { # junction start here

	    if($start<=$end) { push(@{$$readlist[$n][3]},[$start,$end]);}

	    $start=$end+$m[$j]+1;

	    # deal with junction specific stuff
	    if($leftsupport) { # new complete junction
	      if($leftsupport >$maxleftsupport) { # the read is split into more than one exon
		$maxleftsupport=$leftsupport;
	      }
	      my $nj=add_junction($leftend,$rightstart,$maxleftsupport,$support,$junction,$hashjunction,$junctionsupport,$strand,$nh);
	      $newjunction[$njunc]=$nj;
	      $rightsup[$njunc]=$support;
	      $leftsup[$njunc]=$maxleftsupport;
	      $njunc++;

	      push(@{$$readlist[$n][4]},$nj);
		
	    }
	    $rightstart=$start;
	    $leftend=$end;
	    $leftsupport=$support;
	    $support=0;

	    $end=$start-1;
	}
    }
	    
    if($start<=$end) { push(@{$$readlist[$n][3]},[$start,$end]);}

    if($leftsupport) { # update junction
      if($leftsupport >$maxleftsupport) { # the read is split into more than one exon -> if the exon on the left is small than the junction shouldn't be penalized in support
	$maxleftsupport=$leftsupport;
      }
      my $nj=add_junction($leftend,$rightstart,$maxleftsupport,$support,$junction,$hashjunction,$junctionsupport,$strand,$nh);
      $newjunction[$njunc]=$nj;
      $rightsup[$njunc]=$support;
      $leftsup[$njunc]=$maxleftsupport;
      $njunc++;
      push(@{$$readlist[$n][4]},$nj);

    }

    my $maxrightsupport=$rightsup[$njunc-1];
    for(my $j=$njunc-2;$j>=0;$j--) {
      if($rightsup[$j]>$maxrightsupport) { $maxrightsupport=$rightsup[$j];}
      else { # adjust support on the right so that I don't exclude very short internal exons
	if($rightsup[$j]<$junctionsupport && $leftsup[$j]>= $junctionsupport && $maxrightsupport>=$junctionsupport) {
	  $$junction[$newjunction[$j]][4]+=1/$nh;
	}
      }
    }

    if($end<$currentend) { return($currentend);}

    return($end);
}

sub add_junction {
    my ($start,$end,$leftsupport,$rightsupport,$junction,$hashjunction,$junctionsupport,$strand,$nh)=@_;

    my $njunctions=scalar(@{$junction});
    
    my $n=$njunctions;
    
    my $id=$start.":".$end.":".$strand;

    if(defined $$hashjunction{$id}) {
	$n=$$hashjunction{$id};
    }
    else { # this is a new junction
	$$junction[$njunctions++][0]=$strand;
	$$junction[$n][1]=$start;
	$$junction[$n][2]=$end;
	$$junction[$n][4]=0;
	$$hashjunction{$id}=$n;
    }

    #print STDERR "Add junction $n from $start to $end with nh=$nh id=$id and leftsupport=$leftsupport, rightsupport=$rightsupport\n";


    $$junction[$n][3]+=1/$nh;
    if(($leftsupport >= $junctionsupport) && ($rightsupport >= $junctionsupport)) {
	$$junction[$n][4]+=1/$nh;

    }

    return($n);
}

