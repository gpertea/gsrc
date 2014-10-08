ls ../original/  |  perl -ane 'next unless /(\d+)_\w+_L001_R(\d)/; print "ln -s ../original/$F[0] $1_$2.fastq\n";' | sh
