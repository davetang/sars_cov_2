#!/usr/bin/env perl

use strict;
use warnings;
use BeginPerlBioinfo;
use Getopt::Std;

my %opts = ();
getopts('h:f:', \%opts);

if ($opts{'h'} || !exists $opts{'f'}){
   usage();
}

my $fasta = $opts{'f'};
my %lookup = store_fasta($fasta);

foreach my $acc (keys %lookup){
   my $dna = $lookup{$acc};
	my $protein = '';
	my $codon;
	for(my $i=0; $i < (length($dna) - 2) ; $i += 3) {
		$codon = substr($dna,$i,3);
		$protein .= codon2aa($codon);
	}
	print ">$acc\n$protein\n";
}

sub store_fasta {
   my ($infile) = @_;
   my %store = ();
   open(IN, '<', $infile) || die "Could not open $infile: $!\n";
   my $current_id = '';
   while(<IN>){
      chomp;
      if (/^>(.*)/){
         $current_id = $1;
         my @s = split(/\s/, $current_id);
         $current_id = $s[0];
      } else {
         $store{$current_id} .= $_;
      }
   }
   close(IN);
   return(%store);
}

sub usage {
print STDERR <<EOF;
Usage: $0 -f file

Where:   -f         FASTA file
         -h         this helpful usage message

EOF
exit();
}


exit(0);

