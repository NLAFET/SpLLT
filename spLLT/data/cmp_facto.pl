#!/usr/bin/perl

use ValueExtractor;

# directory where runs are kept
my $spllt_output_dir='spllt_starpu';
my $ma87_output_dir='ma87';

# my $facto_time_string = '[>] [spllt_stf_factorize] time:';
my $spllt_facto_time_string = '\[>\] \[factorize\] time:';
my $ma87_facto_time_string  = 'Factor took';

# print "$facto_time_string\n";

# Loop over problems input on stdin
while(<>) {
   chomp;
   if(m/^#/) { next; }
   
   my $prob = $_;
   my $fprob = $prob;
   $fprob =~ s/\//_/g;

   # print "$_\n";
   # print "$fprob\n";
   
   $t_spllt = get_value("${spllt_output_dir}", $fprob, "_NCPU-4_NB-256", $spllt_facto_time_string, "mean");
   $t_ma87  = get_value("${ma87_output_dir}" , $fprob, "_NCPU-4_NB-256", $ma87_facto_time_string,  "mean");
   printf("%10.3e", $t_spllt);
   printf("%10.3e", $t_ma87 );

   # print "$t_spllt\n";
   # print "$t_ma87\n";
   printf("\n");
}
