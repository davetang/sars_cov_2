#!/usr/bin/env perl
#
# Output unique sequences; FASTA identifiers are lost
#

use strict;
use warnings;
use Getopt::Std;

my %opts = ();
getopts('h:f:w:', \%opts);

if ($opts{'h'} ||
    !exists $opts{'f'}
){
   usage();
}

my $wrap_len = 60;
if (exists $opts{'w'}){
   $wrap_len = $opts{'w'};
}

my $fasta_file = $opts{'f'};
my %fasta = store_fasta($fasta_file);
my %seq_hash = ();

my $i = 0;
my $total = 0;
ACC: foreach my $id (keys %fasta){
   ++$total;
   my $seq = $fasta{$id};
   if (exists $seq_hash{$seq}){
      next ACC;
   } else {
      $seq_hash{$seq} = 1;
      ++$i;
      print ">$i\n";
      wrap_seq($seq, $wrap_len);
   }
}

warn("Unique: $i\nTotal: $total\n");

sub wrap_seq {
   my ($seq, $wrap_len) = @_;
   for(my $i = 0; $i <length($seq); $i+=$wrap_len){
      my $ss = substr($seq, $i, $wrap_len);
      print "$ss\n";
   }
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
         -w         length to wrap FASTA sequence (optional, default: 60)
         -h         this helpful usage message

EOF
exit();
}

__END__
