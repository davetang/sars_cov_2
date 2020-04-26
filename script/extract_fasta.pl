#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

my %opts = ();
getopts('i:h:f:', \%opts);

if ($opts{'h'} ||
    !exists $opts{'i'} ||
    !exists $opts{'f'}
){
   usage();
}

my $infile = $opts{'i'};
my $fasta = $opts{'f'};

my %lookup = store_fasta($fasta);

open(IN, '<', $infile) || die "Could not open $infile: $!\n";
while(<IN>){
   chomp;
   my $query = $_;
   if ($lookup{$query}){
      print ">$query\n$lookup{$query}\n";
   } else {
      warn "Could not find $query in $fasta\n";
   }
}
close(IN);

sub store_fasta {
   my ($infile) = @_;
   my %store = ();
   open(IN, '<', $infile) || die "Could not open $infile: $!\n";
   my $current_id = '';
   while(<IN>){
      chomp;
      if (/^>(.*?)\s/){
         $current_id = $1;
         # remove version
         $current_id =~ s/\.\d+$//;
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

Where:   -i         List of IDs
         -f         FASTA file
         -h         this helpful usage message

EOF
exit();
}

__END__
