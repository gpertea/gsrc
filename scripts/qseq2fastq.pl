#!/usr/bin/perl
use strict;
my ($minq, $maxq)=(255,0);
while (<>) {
    chomp;
    my @t = split(/\t/);
    next unless $t[10]==1;
    print '@'."$t[0]:$t[2]:$t[3]:$t[4]:$t[5]#$t[6]/$t[7]\n";
    die ("Error: uneven length !\n") unless length($t[9])==length($t[8]);
    $t[8]=~tr/./N/;
    $t[8]=~s/^N+//;
    my $d=length($t[9])-length($t[8]);
    $t[9]=substr($t[9],$d) if $d;
    $t[8]=~s/N+$//;
    $d=length($t[9])-length($t[8]);
    $t[9]=substr($t[9],0,-$d) if $d;
    die ("Error: uneven length after trimming!\n") 
         unless length($t[9])==length($t[8]);
    my $q=qv($t[9]);
    print "$t[8]\n";
    print "+\n";
    print "$q\n";
}

print STDERR "minq=$minq('".chr($minq)."'), maxq=$maxq('".chr($maxq)."')\n";
$minq-=31;
$maxq-=31;
print STDERR " After conversion:\n";
print STDERR "minq=$minq('".chr($minq)."'), maxq=$maxq('".chr($maxq)."')\n";

sub qv {
 my @s=unpack('C*', $_[0]);
 # assuming phred64 data was given, convert it to phred33
 my $qual='';
 foreach my $c (@s) {
   #my $newc=$c-26;
   my $newc=$c-31;
   $qual .= chr($newc);
   $minq=$c if $c<$minq;
   $maxq=$c if $c>$maxq;
   }
 return $qual; 
}
