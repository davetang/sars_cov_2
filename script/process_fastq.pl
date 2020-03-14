#!/usr/bin/env perl

use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;

my $abs_path = abs_path($0);
$abs_path = dirname($abs_path);

my $thread = 8;
my $metadata = "$abs_path/../sra/metadata.txt";

my @header = ();
my $run;
my $layout;

open(IN, '<', $metadata) || die "Could not open $metadata $!\n";
while(<IN>){
   chomp;
   if ($. == 1){
	   @header = split(/,/);
	   for (my $i = 0; $i <scalar(@header); $i++){
		   if ($header[$i] eq 'LibraryLayout'){
			   $layout = $i;
		   }
		   if ($header[$i] eq 'Run'){
			   $run = $i;
		   }
	   }
	   if (!defined $layout){
		   die "No LibraryLayout information in $metadata\n";
	   }
	   if (!defined $run){
		   die "No Run information in $metadata\n";
	   }
	   next;
   }
   next if /^,/;
   next if /^$header[0]/;
   my @s = split(/,/);
   if ($s[$layout] eq 'PAIRED'){
	   my $fastq1 = $s[$run] . "_1.fastq";
	   my $fastq2 = $s[$run] . "_2.fastq";
	   my $fastq1_base = basename($fastq1, ".fastq");
	   my $fastq2_base = basename($fastq2, ".fastq");
	   if (-f "$abs_path/../raw/fastq/$fastq1" && -f "$abs_path/../raw/fastq/$fastq2"){
		   my $fastp = "fastp --thread $thread -i $abs_path/../raw/fastq/$fastq1 -I $abs_path/../raw/fastq/$fastq2 -o $abs_path/../raw/fastq/${fastq1_base}_fastp.fastq -O $abs_path/../raw/fastq/${fastq2_base}_fastp.fastq";
		   system($fastp);
		   my $bwa = "bwa mem -t $thread $abs_path/../raw/bwa_index/MN908947.fa $abs_path/../raw/fastq/${fastq1_base}_fastp.fastq $abs_path/../raw/fastq/${fastq2_base}_fastp.fastq | samtools sort - -o $abs_path/../result/$s[$run]_MN908947.bam";
		   system($bwa);
		   my $samtools  = "samtools view -F 0x804 -f 2 -b $abs_path/../result/$s[$run]_MN908947.bam > $abs_path/../result/$s[$run]_MN908947_mapped.bam";
		   system($samtools);
		   my $bcftools = "bcftools mpileup -f $abs_path/../raw/MN908947.fa $abs_path/../result/$s[$run]_MN908947_mapped.bam | bcftools call -mv -Ov -o $abs_path/../result/$s[$run]_MN908947_mapped.vcf";
		   system($bcftools);
	   }
   } elsif ($s[$layout] eq 'SINGLE'){
	   if (-f "$abs_path/../raw/fastq/$s[$run].fastq"){
		   my $fastp = "fastp --thread $thread -i $abs_path/../raw/fastq/$s[$run].fastq -o $abs_path/../raw/fastq/$s[$run]_fastp.fastq";
		   system($fastp);
		   my $bwa = "bwa mem -t $thread $abs_path/../raw/bwa_index/MN908947.fa $abs_path/../raw/fastq/$s[$run]_fastp.fastq | samtools sort - -o $abs_path/../result/$s[$run]_MN908947.bam";
		   system($bwa);
		   my $samtools  = "samtools view -F 0x804 -b $abs_path/../result/$s[$run]_MN908947.bam > $abs_path/../result/$s[$run]_MN908947_mapped.bam";
		   system($samtools);
		   my $bcftools = "bcftools mpileup -f $abs_path/../raw/MN908947.fa $abs_path/../result/$s[$run]_MN908947_mapped.bam | bcftools call -mv -Ov -o $abs_path/../result/$s[$run]_MN908947_mapped.vcf";
		   system($bcftools);
	   }
   } else {
	   die "Unrecognised LibraryLayout: $s[$layout]\n";
   }

}
close(IN);

__END__
