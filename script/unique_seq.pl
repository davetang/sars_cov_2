#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

my %opts = ();
getopts('i:h:f:', \%opts);

if ($opts{'h'} ||
    !exists $opts{'f'}
){
   usage();
}

my $infile = $opts{'i'};
my $fasta_file = $opts{'f'};

my %fasta = store_fasta($fasta_file);
my %seq_hash = ();

foreach my $id (keys %fasta){
   if (exists $seq_hash{$fasta{$id}}){
      $seq_hash{$fasta{$id}}++;
   } else {
      $seq_hash{$fasta{$id}} = 1;
   }
}

foreach my $seq (sort {$seq_hash{$b} <=> $seq_hash{$a}} keys %seq_hash){
   print "$seq\t$seq_hash{$seq}\n";
}

sub store_fasta {
   my ($infile) = @_;
   my %store = ();
   my $fh;
   if ($infile =~ /\.gz$/){
      open($fh, '-|', "gunzip -c $infile") || die "Could not open $infile: $!\n";
   } else {
      open($fh, '<', $infile) || die "Could not open $infile: $!\n";
   }
   my $current_id = '';
   while(<$fh>){
      chomp;
      if (/^>(.*?)\s/){
         $current_id = $1;
      } else {
         $store{$current_id} .= $_;
      }
   }
   close($fh);
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

__END__
