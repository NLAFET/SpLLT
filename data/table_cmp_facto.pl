#!/usr/bin/perl

use ValueExtractor;
use LatexPrint;
use List::Util qw(min);

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
   my @t_spllt;
   my @t_ma87;
   
   # push @t_spllt, {'64' => get_value("${spllt_output_dir}", $fprob, "_NCPU-27_NB-64" , $spllt_facto_time_string, "mean")};
   # push @t_spllt, {'128' => get_value("${spllt_output_dir}", $fprob, "_NCPU-27_NB-128", $spllt_facto_time_string, "mean")};
   # push @t_spllt, {'256' => get_value("${spllt_output_dir}", $fprob, "_NCPU-27_NB-256", $spllt_facto_time_string, "mean")};
   # push @t_spllt, {'512' => get_value("${spllt_output_dir}", $fprob, "_NCPU-27_NB-512", $spllt_facto_time_string, "mean")};

   # push @t_ma87 , {'64' => get_value("${ma87_output_dir}" , $fprob, "_NCPU-28_NB-64" , $ma87_facto_time_string , "mean")};
   # push @t_ma87 , {'128' => get_value("${ma87_output_dir}" , $fprob, "_NCPU-28_NB-128", $ma87_facto_time_string , "mean")};
   # push @t_ma87 , {'256' => get_value("${ma87_output_dir}" , $fprob, "_NCPU-28_NB-256", $ma87_facto_time_string , "mean")};
   # push @t_ma87 , {'512' => get_value("${ma87_output_dir}" , $fprob, "_NCPU-28_NB-512", $ma87_facto_time_string , "mean")};

   # my $best_t_spllt = $t_spllt->{ ( sort {$a <=> $b} keys %$t_spllt )[0] };
   # my $best_t_ma87  = $t_ma87->{ ( sort {$a <=> $b} keys %$t_ma87 )[0] };

   push @t_spllt, get_value("${spllt_output_dir}", $fprob, "_NCPU-27_NB-64" , $spllt_facto_time_string, "mean");
   push @t_spllt, get_value("${spllt_output_dir}", $fprob, "_NCPU-27_NB-128", $spllt_facto_time_string, "mean");
   push @t_spllt, get_value("${spllt_output_dir}", $fprob, "_NCPU-27_NB-256", $spllt_facto_time_string, "mean");
   push @t_spllt, get_value("${spllt_output_dir}", $fprob, "_NCPU-27_NB-512", $spllt_facto_time_string, "mean");

   push @t_ma87 , get_value("${ma87_output_dir}" , $fprob, "_NCPU-28_NB-64" , $ma87_facto_time_string , "mean");
   push @t_ma87 , get_value("${ma87_output_dir}" , $fprob, "_NCPU-28_NB-128", $ma87_facto_time_string , "mean");
   push @t_ma87 , get_value("${ma87_output_dir}" , $fprob, "_NCPU-28_NB-256", $ma87_facto_time_string , "mean");
   push @t_ma87 , get_value("${ma87_output_dir}" , $fprob, "_NCPU-28_NB-512", $ma87_facto_time_string , "mean");

   my $best_t_spllt = min @t_spllt;
   my $best_t_ma87  = min @t_ma87;

   printf("%-35s", latex_escape($prob));
   # printf("%30s", $prob);
   # printf("%12.3e", $t_spllt);
   printf(" & ");
   # printf("%12.3e", $t_spllt);
   latex_float($best_t_spllt, 3, ($best_t_spllt<$best_t_ma87));

   printf(" & ");
   # printf("%12.3e", $t_ma87 );
   latex_float($best_t_ma87, 3, ($best_t_ma87<$best_t_spllt));

   # print "$t_spllt\n";
   # print "$t_ma87\n";
   # printf("\n");
   printf(" \\\\\n");
}
