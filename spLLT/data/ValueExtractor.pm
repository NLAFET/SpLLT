package ValueExtractor;
use strict;
use warnings;

use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
$VERSION = 1.00;
@ISA = qw(Exporter);
@EXPORT = qw(get_value);
@EXPORT_OK = qw();

use List::Util qw( sum min max );
use Switch;

# Usage:
#     get_value(dir, prob, nproc, pattern[, type])
#
# Will look in the file ${dir}${run}/${prob}${post} for the last line matching
# ${pattern} and pull the last number off that line, averaging it over files.
# ${run} here varies over '', '_1', '_2', ... '_5'.
# 
# Type can be "first", "min", "max", "mean", or "median". default is "mean".
sub get_value
{
   my $nparam = @_;
   my $dir=$_[0];
   my $fprob=$_[1];
   my $post=$_[2];
   my $pattern=$_[3];
   my $which = "mean";
   if($nparam >= 5) { $which=$_[4]; }

   my @times;
   for(my $run=-1; $run<=9; $run++) {
      my $file;
      if($run==-1) { 
         $file = "$dir/${fprob}$post";
      } else {
         $file = "$dir\_$run/$fprob$post";
      }
      if ( -e "$file" ) {
         my $t = `grep -P "$pattern" $file`;
         # print "$t\n";
         if(${^CHILD_ERROR_NATIVE} == 0) {
            chomp ($t);
            $t = `echo '$t' | tail -1`;
            $t =~ s/${pattern}[^0-9]*([0-9E.+-]*).*$/$1/g;
            chomp($t);
            if($t eq "") { next; }
            if($which eq "first") { return $t; }
            push @times, $t;
         }
      }
   }

   if(@times==0) {
      return "Inf";
   }
   switch ($which) {
      case "mean"    { return sum(@times)/@times; }
      case "median"  {
         my @times_sorted = sort @times;
         return $times_sorted[@times/2];
      }
      case "min"     { return min(@times); }
      case "max"     { return max(@times); }
      case "first"   { return $times[0]; }
   }

   print "Unkown operation '$which'";
   return "Inf";
}

1;
