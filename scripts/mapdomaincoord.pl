#!/usr/bin/perl

($proteinsearches,$protdom)=@ARGV;

open(F,$proteinsearches);
while(<F>) {
    chomp;
    @a=split;
    $name=substr($a[0],1);
    $beg{$name}=$a[1];
    $end{$name}=$a[2];
    <F>;
}
close(F);

open(F,$protdom);
while(<F>) {
    chomp;
    @a=split;
    @b=split(/_/,$a[0]);
    $newname=$b[0];
    for($i=1;$i<$#b;$i++) {
	$newname.="_".$b[$i];
    }
    
    if($beg{$a[0]}<$end{$a[0]}) {
	$newbeg=$beg{$a[0]}+($a[1]-1)*3;
	$newend=$beg{$a[0]}+($a[2]-1)*3+2;
    }
    else {
	$newbeg=$beg{$a[0]}-($a[1]-1)*3;
	$newend=$beg{$a[0]}-($a[2]-1)*3-2;
    }
    printf("%s %20d %10d %20s %5d %5d %5s %10s %s\n",$newname,$newbeg,$newend,$a[3],$a[4],$a[5],$a[6],$a[7],$a[8],$a[9]);
}
close(F);

