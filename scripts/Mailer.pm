#!/usr/bin/perl
use strict 'vars';
#use Mail::Mailer;
use POSIX 'uname';
use IO::File;

#------- retrieve the local host name

my $thisHost= (&POSIX::uname)[1]; # can't trust HOST envvar under condor, because for getenv=TRUE
                          # it will either inherit HOST and HOSTNAME from submitter env
                          # or not set anything for getenv=FALSE
                          # --> so we use POSIX uname() call here
chomp($thisHost);
$thisHost=lc($thisHost);
my ($thisDomain)=($thisHost=~/^[\w\-]+\.(.+)/);

#my $def_domain = 'jimmy.harvard.edu';
my $def_domain = $thisDomain;
$def_domain=~s/^dfci\./jimmy\./;
my $def_host = $thisHost;

#----------------------
sub send_mail {
 my $hash=shift;
 $hash->{'from'}=$ENV{'USER'}.'@'.$def_host
    unless defined($hash->{'from'});
 my $to=$hash->{'to'};
 unless ($to) {
    $hash->{'to'}=$ENV{'USER'}.'@'.$def_domain;
    }
   else {
    $hash->{'to'}=$to.'@'.$def_domain unless $to =~ m/@/;
    }
    
 my $file;
 local *ADDFILE;
 if (defined($hash->{file})) {
   #warning: assumes it's a decent, short text file!
   local $/=undef; #read whole file
   open(ADDFILE, '<'.$hash->{file}) || return "Cannot open file ".$hash->{file}."\n";
   $file=<ADDFILE>;
   close ADDFILE;
   }
 my $body = $hash->{'body'};
 $body.=$file;
 my $fh;
 $fh = IO::File->new('| /usr/lib/sendmail -t -oi') ||
    die("Mailer.pm error: Cannot open the sendmail pipe\n");
 $fh->print("To: $hash->{to}\n");
 $fh->print("From: $hash->{from}\n");
 $fh->print("Subject: $hash->{subj}\n\n");
 $fh->print($body);
 $fh->close();

 #-- create the Mail::Mailer object and send the message:
 # my $mailer = Mail::Mailer->new();
 # $mailer->open({ 'From'    => $hash->{'from'},
 #                'To'      => $hash->{'to'},
 #                'Cc'      => $hash->{'cc'},
 #                'Subject' => $hash->{'subj'}
 #              })
 #    or die "Mailer.pm error: Can't open. $!\n";
 #print $mailer $body;
 #$mailer->close();
}
