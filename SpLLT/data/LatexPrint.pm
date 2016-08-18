package LatexPrint;
use strict;
use warnings;

use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
$VERSION = 1.00;
@ISA = qw(Exporter);
@EXPORT = qw(latex_escape latex_decimal latex_exp latex_float);
@EXPORT_OK = qw();

use Math::SigFigs;

sub latex_bold_italic {
   my $bold    = $_[0];
   my $italic  = $_[1];

   if($bold) {
      printf("\\bf");
   } else {
      if($italic) {
         printf("\\it");
      }
   }
}

sub latex_escape
{
   $_ = $_[0];
   s/_/\\_/g;
   return $_;
}

sub latex_decimal
{
   my $val=$_[0];
   my $bold=$_[1];
   my $italic=$_[2];

   latex_bold_italic($bold, $italic);
   printf("%d", $val);
}

sub latex_float
{
   my $val     = $_[0];
   my $sf      = $_[1];
   my $bold    = $_[2];
   my $italic  = $_[3];

   latex_bold_italic($bold, $italic);
   if($val < 1) {
      printf("%.2f", $val);
   } else {
      my $str = FormatSigFigs($val, $sf);
      $str =~ s/\.$//;
      print $str;
   }
}

sub latex_exp
{
   my $val=$_[0];
   my $bold=$_[1];
   my $italic=$_[2];

   my $str = sprintf("%10.2e", $val);
   my @parts = split /e/, $str;

   if($bold) {
      printf("\\bf%.2f\$\\mathbf{\\times10^{%d}}\$", $parts[0], $parts[1]);
   } else {
      if($italic) {
         printf("\\it%.2f\$\\mathit{\\times10^{%d}}\$", $parts[0], $parts[1]);
      } else {
         printf("%.2f\$\\times10^{%d}\$", $parts[0], $parts[1]);
      }
   }
}

1;
