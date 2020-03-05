#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

my %opts = ();
getopts('h:f:', \%opts);

if ($opts{'h'} || !exists $opts{'f'}){
   usage();
}

my $fasta = $opts{'f'};
my %fasta = store_fasta($fasta);

print join("\t", "Accession", "Length", "A", "C", "G", "T", "Unknown"), "\n";

foreach my $acc (keys %fasta){
   my $seq = $fasta{$acc};
   $seq = uc($seq);
   my $len = length($seq);
   my $a = my $c = my $g = my $t = $seq;
   $a =~ s/[^A]//g;
   $c =~ s/[^C]//g;
   $g =~ s/[^G]//g;
   $t =~ s/[^T]//g;
   my $ulen = $len - length($a) - length($c) - length($g) - length($t);
   print join("\t", $acc, $len, length($a), length($c), length($g), length($t), $ulen), "\n";
}

sub store_fasta {
   my ($infile) = @_;
   my %store = ();
   open(IN, '<', $infile) || die "Could not open $infile: $!\n";
   my $current_id = '';
   while(<IN>){
      chomp;
      if (/^>(.*?)\s/){
         $current_id = $1;
      } else {
         $store{$current_id} .= $_;
      }
   }
   close(IN);
   return(%store);
}

sub usage {
print STDERR <<EOF;
Usage: $0 -i file -f file

Where:   -f         FASTA file
         -h         this helpful usage message

EOF
exit();
}

__END__
