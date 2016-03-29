#!/usr/bin/perl

use ValueExtractor;

# directory where runs are kept
my $OUTPUT_DIR='spllt_starpu';

# my $facto_time_string = '[>] [spllt_stf_factorize] time:';
my $facto_time_string = '\[>\] \[factorize\] time:';

# print "$facto_time_string\n";

# Loop over problems input on stdin
while(<>) {
   chomp;
   if(m/^#/) { next; }
   
   my $prob = $_;
   my $fprob = $prob;
   $fprob =~ s/\//_/g;

   print "$_\n";
   print "$fprob\n";
   
   $t = get_value("${OUTPUT_DIR}", $fprob, "_NCPU-4_NB-256", $facto_time_string, "mean");
   print "$t\n";
}
