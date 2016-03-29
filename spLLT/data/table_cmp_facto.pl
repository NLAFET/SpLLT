#!/usr/bin/perl

use ValueExtractor;
use LatexPrint;

# directory where runs are kept
my $spllt_output_dir='cn255/spllt_starpu';
my $ma87_output_dir='cn255/ma87';

# my $facto_time_string = '[>] [spllt_stf_factorize] time:';
my $spllt_facto_time_string = '\[>\] \[factorize\] time:';
my $ma87_facto_time_string  = 'Factor took';

# print "$facto_time_string\n";

my $suffix = "_NCPU-27_NB-256";

printf("%-35s%12s%12s", " ", " & spLLT", "& MA87");
printf(" \\\\\n");

# Loop over problems input on stdin
while(<>) {
   chomp;
   if(m/^#/) { next; }
   
   my $prob = $_;
   my $fprob = $prob;
   $fprob =~ s/\//_/g;

   # print "$_\n";
   # print "$fprob\n";
   
   $t_spllt = get_value("${spllt_output_dir}", $fprob, "_NCPU-27_NB-256", $spllt_facto_time_string, "mean");
   $t_ma87  = get_value("${ma87_output_dir}" , $fprob, "_NCPU-27_NB-256", $ma87_facto_time_string , "mean");

   printf("%-35s", latex_escape($prob));
   # printf("%30s", $prob);
   # printf("%12.3e", $t_spllt);
   printf(" & ");
   printf("%12.3e", $t_spllt);
   # latex_float($t_spllt, 3);
   printf(" & ");
   printf("%12.3e", $t_ma87 );

   # print "$t_spllt\n";
   # print "$t_ma87\n";
   # printf("\n");
   printf(" \\\\\n");
}
