#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

my %opts = ();
getopts('i:h:p:l:f:', \%opts);

if ($opts{'h'} ||
    !exists $opts{'i'} ||
    !exists $opts{'p'} ||
    !exists $opts{'l'} ||
    !exists $opts{'f'}
){
   usage();
}

my $infile = $opts{'i'};
my $fasta = $opts{'f'};
my $perc_thres = $opts{'p'};
my $len_thres = $opts{'l'};

my %lookup = store_defline($fasta);

open(IN, '<', $infile) || die "Could not open $infile: $!\n";
while(<IN>){
   chomp;
   next if /^#/;
	my ($query, $subject, $perc_id, $align_len, $mismatches, $gap_opens, $qstart, $qend, $sstart, $send, $evalue, $bit_score) = split(/\t/);
   if ($lookup{$subject}){
      if ($perc_id < $perc_thres || $align_len < $len_thres){
         next;
      } else {
         print join("\t", $subject, $align_len, $perc_id, $mismatches, $lookup{$subject}), "\n";
      }
   } else {
      die "Could not find $subject in $fasta\n";
   }

}
close(IN);

sub store_defline {
   my ($infile) = @_;
   my %store = ();
   open(IN, '<', $infile) || die "Could not open $infile: $!\n";
   while(<IN>){
      chomp;
      next unless /^>/;
      if (/>(.*?)\s(.*)/){
         my $query = $1;
         my $def = $2;
         $store{$query} = $def;
      }
   }
   close(IN);
   return(%store);
}

sub usage {
print STDERR <<EOF;
Usage: $0 -i blastout.tsv -f seq.fa -p 90 -l 100

Where:   -i                BLAST output outfmt 7
         -f                FASTA file used for BLAST db
         -p                Percentage identity threshold
         -l                Alignment length threshold
         -h                this helpful usage message

EOF
exit();
}

__END__
