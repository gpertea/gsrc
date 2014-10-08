#package minicgi;
use strict;

our ($VERSION, @ISA, @EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw( cgi_param cgi_start cgi_end cgi_redir 
        cgi_htpl ptag ptag_start ptag_end);

=head1 NAME

minicgi - a simple replacement of the the bloated CGI module

=cut

our ( $contype, $doctype, %qparams );
open(STDERR, ">&STDOUT"); #just so we can see the error messages!
#my $ptag_level=0;
my @ptag_stack;
$contype='html';
#$doctype='HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd"';
$doctype='HTML PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd"';
my $reqmethod=$ENV{REQUEST_METHOD};
my $qrystr;
if ($reqmethod) {
 if ($reqmethod eq 'GET') {
    $qrystr = $ENV{QUERY_STRING};
 } elsif ($reqmethod eq 'POST') {
    read(STDIN, $qrystr, $ENV{CONTENT_LENGTH});
 }
}

if ($qrystr) {
 my @nvpairs = split(/&/, $qrystr);
 foreach my $nv (@nvpairs) {
    my ($n, $v) = split(/=/, $nv);
    $v =~ tr/+/ /; 
    $v =~ s/%([\dA-Fa-f][\dA-Fa-f])/pack("C", hex($1))/eg;
    ## unpack non-printables (represented as '%00' - '%FF'  hex values
    ## in URLs) back to original values
    $qparams{$n} = $v;
    }
}    

return 1;


=cut

=head2 cgi_param([<param_name>..])
 
 Returns names or values of CGI GET or POST parameters.
 - if <param_name> is given, its value is returned if there
 - if no parameters are given, the list of all 
   CGI parameter names is returned
 - if more than one <param_name> is given, a list of values
   is returned in list context (or only the first value 
   in scalar context)

=cut

sub cgi_param {
my @k=@_;
if (@k==0 && wantarray()) {
 my @r=keys(%qparams);
 return @r;
 }
my @r = map { $qparams{$_} } @k;
return wantarray() ? @r : $r[0]; 
}


=head2 cgi_start([<option> => 'strvalue', ..])
 
 Starts a document -- also sends the content-type header.
 The <head> and the opening <html> and <body> tags are also sent
 unless the 'nohead' option is given.
 
 Options accepted: 
     ctype =>'html'/'text'/'gif'/'jpg'/'png' 
               (sets content-type; default: html)
 
 The following options are only considered if ctype is html:

   doctype =>"custom-doctype-string" (default: XHTML 1.0 Strict)
     nohead=>1 this will disable sending of any html data after the doctype string
     title =>"pagetitle" sets the content of the <title> tag inside <head>
    onload =>"javascript_code" code to use for <body onLoad="javascript_code">)
     body =>"body_attrs.." other attributes for the BODY tag <body body_attrs..>
    head =>"additional_head_content" (you can place all your 
             javascript and CSS related stuff in here, literally)
                
               
=cut

sub cgi_start {
 my (%hp)=@_;
 my $ctype=$contype;
 my $dtype=$doctype;
 my ($nohead, $title, $head, $onload, $battrs);
 my @k=keys(%hp);
 if (@k) {  
  $ctype=$hp{'ctype'} || $contype;
  $dtype=$hp{'doctype'} || $doctype;
  $title=$hp{'title'};
  $head=$hp{'head'};
  $battrs=$hp{'body'};
  $nohead=1 unless $head;
  $onload=$hp{'onload'};
  }
  
 #---
 if ($ctype eq 'html') {
    print STDOUT "Content-type: text/html\n\n";
    }
  elsif ($ctype eq 'text') {
    print STDOUT "Content-type: text/plain\n\n";
    return;
    }
  else {
    print STDOUT "Content-type: image/$ctype\n\n";
    return;
    }
 #only for html ctype:   
 print STDOUT "<!DOCTYPE $dtype>\n";
 #add optional head content and body onLoad="<script>" here :
 return if ($nohead);
 print STDOUT "<html>\n";
 print STDOUT "<head>\n";
 #--print meta content type info
 print STDOUT q{ <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
 <meta http-equiv="Content-Style-Type" content="text/css" />
 <meta http-equiv="Content-Script-Type" content="text/javascript" />}."\n";
 if ($title) {
   $title=~tr/<>//d;
   print STDOUT "<title>$title</title>";
   }
 print STDOUT $head if ($head);
 print STDOUT "</head>\n";
 $battrs=' '.$battrs if $battrs;
 print STDOUT ($onload ? qq{<body onload="$onload"$battrs>\n} : "<body$battrs>\n");
}


=head2 cgi_end()

 Simply sends the html closing tags for the page:
 </body> and </html>

=cut 

sub cgi_end {
 print STDOUT "\n</body>\n</html>";
}


=head2 cgi_redir(<url>)

 Redirects the response to a the given <url>.

=cut 

sub cgi_redir {
 #print STDOUT "HTTP/1.1 301 Moved Permanently\n\n";
 #print STDOUT "Status: 301 Moved Permanently\n\n";
 my $url=$_[0];
 if (index($url,'/')==0) {
  $url='http://'.$ENV{'HTTP_HOST'}.$url;
  }
 print STDOUT "Location: $url\n\n";
}

=head2 cgi_htpl(filename, [varname=>replvalue, ...])

 Loads the provided html template file and replaces all 
 template variables accordingly.
 
 Template variables in the file must follow this format:

   ${varname}

 If a variable is found in the template file but its 
 replacement value is not provided among cgi_htpl parameters, 
 it will be simply replaced with empty space.

=cut 

sub cgi_htpl {
 my $fname=shift(@_);
 if ($fname!~m/^[\/\\]/) {
  # if it doesn't start with an absolute path marker, take it as relative to DOCUMENT_ROOT
  $fname=$ENV{'DOCUMENT_ROOT'}.'/'.$fname;
  }
 local *HTPL;
 open(HTPL, $fname) || die("Error opening html template file $fname!\n");
 local $/=undef;
 my $ftxt=<HTPL>;
 close(HTPL);
 $ftxt=~s/<!\s*--.*?--+\s*>//sg; #remove commented html
 if (@_) {
   my (%hv)=@_;
   my @tplvars=($ftxt=~m/\$\{(\w+)\}/sg);
   my %h; @h{@tplvars}=(); @tplvars=keys(%h);
   foreach my $v (@tplvars) {
      my $r=$hv{$v} || '';
      $ftxt=~s/\$\{$v\}/$r/sg;
      }
   }
  else { #no values given, just remove all vars
   $ftxt=~s/\$\{\w+\}//sg;
   }
 print STDOUT $ftxt;
}


=head2 ptag('tag', 'text' [, HTMLattrs..])
 | ptag('tag/', [, HTMLattrs..])

 Print HTML tag block with text content and optional HTML attributes.
 By default text is printed between start and closing tags like this:
 
 <tag HTMLattrs.. >
  text
 </tag>
 
 A self-closing no content tag will be printed if the tag ends with 
 the '/' character, in which case all other parameters are expected 
 to be HTML attributes for the tag, e.g.
 
   ptag('img/', 'class="topleft"', 'src="/img/jhmilogo_wbg.png"');
 
=cut

sub ptag {
 my $tag=shift(@_);
 my $noct;
 { local $/="/"; $noct=chomp($tag); }
 my $ct=shift(@_) unless $noct;
 print STDOUT '<'.$tag;
 print STDOUT ' '.join(' ', @_).' ' if @_; #attrs
 print $noct ? '/>' :
       '>'.$ct.'</'.$tag.">\n";
}

=head2 ptag_start('tag' [, HTMLattrs..])
  ... ptag_end()

 ptag_start('tag') opens a HTML tag, with optional HTMLattrs
 Tags are placed in a stack such that ptag_end() will 
 automatically pop the last open tag and close it.
 Example:
  
  ptag_start('div', 'id="addrdiv"');
  print STDOUT "Address is $addr";
  ptag_end();
 
=cut

sub ptag_start {
 my $tag=shift(@_);
 my $i=scalar(@ptag_stack);
 push(@ptag_stack, $tag);
 print STDOUT ' 'x$i if $i;
 print STDOUT '<'.$tag;
 print STDOUT ' '.join(' ', @_).' ' if @_; #attrs
 print ">\n";
}

sub ptag_end {
 return unless @ptag_stack;
 my $tag=pop(@ptag_stack);
 my $i=scalar(@ptag_stack);
 print STDOUT ' 'x$i if $i;
 print STDOUT '<'.$tag."/>\n";
}
