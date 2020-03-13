## README

The [Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)](https://en.wikipedia.org/wiki/Severe_acute_respiratory_syndrome_coronavirus_2) is currently causing a pandemic. This repository contains scripts and analysis code for analysing the sequence of SARS-CoV-2.

## Install

* [Entrez Direct (EDirect)](https://www.ncbi.nlm.nih.gov/books/NBK179288/)

```bash
conda env create --file environment.yml
```

### Entrez Direct Functions

Navigation functions support exploration within the Entrez databases:

* esearch performs a new Entrez search using terms in indexed fields.
* elink looks up neighbors (within a database) or links (between databases).
* efilter filters or restricts the results of a previous query.

Records can be retrieved in specified formats or as document summaries:

* efetch downloads records or reports in a designated format.

Desired fields from XML results can be extracted without writing a program:

* xtract converts EDirect XML output into a table of data values.

Several additional functions are also provided:

* einfo obtains information on indexed fields in an Entrez database.
* epost uploads unique identifiers (UIDs) or sequence accession numbers.
* nquire sends a URL request to a web page or CGI service.

## Sequences

* SARS-CoV-2 - https://www.ncbi.nlm.nih.gov/nuccore/MN908947

```bash
mkdir raw

efetch -db sequences -format genbank -id MN908947 > raw/MN908947.genbank
efetch -db sequences -format fasta -id MN908947 > raw/MN908947.fa
efetch -db sequences -format fasta_cds_aa -id MN908947
```

41874 results as of 2020/03/01.

```bash
esearch -db nuccore -query coronavirus
<ENTREZ_DIRECT>
  <Db>nuccore</Db>
  <WebEnv>NCID_1_121026101_130.14.22.33_9001_1583054896_904858778_0MetA0_S_MegaStore</WebEnv>
  <QueryKey>1</QueryKey>
  <Count>41874</Count>
  <Step>1</Step>
</ENTREZ_DIRECT>
```

Get FASTA.

```bash
esearch -db nuccore -query coronavirus | efetch -db sequences -format fasta > raw/coronavirus_20200301.fa

esearch -db nuccore -query coronavirus | efetch -db sequences -format fasta > raw/coronavirus_20200306.fa

cat raw/coronavirus_20200301.fa | grep "^>" | wc -l
   41874

cat raw/coronavirus_20200301.fa | grep MN908947
>MN908947.3 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome
```

FASTA stats.

```bash
script/fasta_stats.pl -f raw/coronavirus_20200301.fa | gzip > result/coronavirus_20200301_stat.txt.gz
```

## BLAST

Create BLAST database.

```bash
mkdir db

makeblastdb -dbtype nucl \
            -in raw/coronavirus_20200301.fa \
            -input_type fasta \
            -title coronavirus_20200301 \
            -out db/coronavirus_20200301

Building a new DB, current time: 03/01/2020 19:25:30
New DB name:   /Users/dtang/github/sars_cov_2/db/coronavirus_20200301
New DB title:  coronavirus_20200301
Sequence type: Nucleotide
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 41874 sequences in 9.91486 seconds.
```

BLAST.

```bash
# -evalue <Real> - Expectation value (E) threshold for saving hits 
# Default = `10'
blastn -outfmt 7 -query raw/MN908947.fa -db db/coronavirus_20200301 > result/MN908947_blast.txt

blastn -evalue 1 -outfmt 7 -query raw/MN908947.fa -db db/coronavirus_20200301 | wc -l

# accessions with BLAST hits
cat result/MN908947_blast.txt | grep -v "^#" | cut -f2 | sort -u | wc -l
     500

cat result/MN908947_blast.txt | grep -v "^#" | cut -f2 | sort -u > result/MN908947_matched.txt
```

Extract FASTA.

```bash
script/extract_fasta.pl -i result/MN908947_matched.txt -f raw/coronavirus_20200301.fa > result/MN908947_matched.fa

cat result/MN908947_matched.fa | grep "^>" | wc -l
     500
```

FASTA stats.

```bash
script/fasta_stats.pl -f result/MN908947_matched.fa > result/MN908947_matched_stats.txt
```

## Parse results

```bash
script/parse_outfmt7.pl -i result/MN908947_blast.txt -p 80 -l 10000 -f raw/coronavirus_20200301.fa | less
```

## ClustalW

```bash
script/extract_fasta.pl -i raw/wanted.txt -f raw/coronavirus_20200301.fa > raw/wanted.fa
clustalw -infile=raw/wanted.fa

script/extract_fasta.pl -i raw/MN908947_MN996532.txt -f raw/coronavirus_20200301.fa > raw/MN908947_MN996532.fa
clustalw -infile=raw/MN908947_MN996532.fa
```

## SRA

Download YML file with SARS-CoV-2 accessions.

```bash
wget -N https://www.ncbi.nlm.nih.gov/core/assets/genbank/files/ncov-sequences.yaml -O raw/ncov-sequences.yaml
```

Use `wget` To [obtain metadata](https://www.ncbi.nlm.nih.gov/books/NBK242621/).

```bash
mkdir sra

# sra-study: SRP242226
wget -O sra/SRP242226_info.csv 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=SRP242226'
```

Checkout metadata using [csvkit](https://csvkit.readthedocs.io/en/latest/).

```bash
csvcut -n sra/SRP242226_info.csv  | head -6
1: Run
2: ReleaseDate
3: LoadDate
4: spots
5: bases
6: spots_with_mates

csvcut -c Run,Sample,BioSample,spots,LibraryStrategy,LibrarySource,LibraryLayout,Platform,Model,SampleType sra/SRP242226_info.csv | csvlook
| Run         | Sample     | BioSample    |   spots | LibraryStrategy | LibrarySource      | LibraryLayout | Platform | Model          | SampleType |
| ----------- | ---------- | ------------ | ------- | --------------- | ------------------ | ------------- | -------- | -------------- | ---------- |
| SRR10903402 | SRS6007143 | SAMN13872786 | 676,694 | RNA-Seq         | METATRANSCRIPTOMIC | PAIRED        | ILLUMINA | Illumina MiSeq | simple     |
| SRR10903401 | SRS6007144 | SAMN13872787 | 476,632 | RNA-Seq         | METATRANSCRIPTOMIC | PAIRED        | ILLUMINA | Illumina MiSeq | simple     |
|             |            |              |         |                 |                    |               |          |                |            |
```

Download all metadata using SRPs from ncov-sequences.yaml.

```bash
for acc in `cat raw/ncov-sequences.yaml | grep sra-study | sort -u | cut -f3 -d' '`; do
   wget -O sra/${acc}_info.csv "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=$acc"
done

rm -f sra/metadata.txt
touch sra/metadata.txt
for file in `ls sra/*info.csv`; do
   echo $file;
   csvcut -c Run,Sample,BioSample,spots,LibraryStrategy,LibrarySource,LibraryLayout,Platform,Model,SampleType $file >> sra/metadata.txt
done
```

Download all Illumina data.

```bash
for acc in `cat sra/metadata.txt | grep ILLUMINA | cut -f1 -d','`; do
   echo $acc
   fasterq-dump -p --outdir raw/fastq $acc
done
```

### SRR10971381

Follow https://github.com/galaxyproject/SARS-CoV-2/blob/master/1-PreProcessing/pp_wf.png

Use [SRA Toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc) to download FASTQ sequences from the SRA. First use [prefetch](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=prefetch), a command-line for downloading SRA, dbGaP, and ADSP data, to download `sra` files. See https://www.ncbi.nlm.nih.gov/books/NBK242621/ for more information.

Download and install https://downloads.asperasoft.com/connect2/.

```bash
prefetch --output-directory raw SRR10971381

# connection keeps dying
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra46/SRR/010714/SRR10971381

md5sum SRR10971381 > SRR10971381.md5sum
cat SRR10971381.md5sum
5496488662893a836e23541b84bfb7cd  SRR10971381
```

Uploaded to https://davetang.org/file/SRR10971381, so download file from my server.

```bash
wget -c -N https://davetang.org/file/SRR10971381  
wget -c -N https://davetang.org/file/SRR10971381.md5sum

md5sum -c SRR10971381.md5sum
SRR10971381: OK
```

Next use [fastq-dump](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump) to convert SRA data into FASTQ.

```bash
fastq-dump --split-files ./SRR10971381
2020-03-10T14:51:00 fastq-dump.2.10.3 err: name not found while resolving query within virtual file system module - failed to resolve accession './SRR10971381' - Cannot resolve accession ( 404 )
Read 28282964 spots for ./SRR10971381
Written 28282964 spots for ./SRR10971381
```

[fastp](https://github.com/OpenGene/fastp).

```bash
fastp --thread 8 -i SRR10971381_1.fastq -I SRR10971381_2.fastq -o SRR10971381_1_fastp.fastq -O SRR10971381_2_fastp.fastq
Read1 before filtering:
total reads: 28282964
total bases: 4017125680
Q20 bases: 1783384314(44.3945%)
Q30 bases: 1735659038(43.2065%)

Read2 before filtering:
total reads: 28282964
total bases: 4013917534
Q20 bases: 1723725994(42.9437%)
Q30 bases: 1652844944(41.1779%)

Read1 after filtering:
total reads: 13054241
total bases: 1786633510
Q20 bases: 1671652872(93.5644%)
Q30 bases: 1634552420(91.4878%)

Read2 aftering filtering:
total reads: 13054241
total bases: 1782180210
Q20 bases: 1625652911(91.2171%)
Q30 bases: 1568126467(87.9892%)

Filtering result:
reads passed filter: 26108482
reads failed due to low quality: 30441256
reads failed due to too many N: 16190
reads failed due to too short: 0
reads with adapter trimmed: 582728
bases trimmed due to adapters: 30162896

Duplication rate: 5.57505%

Insert size peak (evaluated by paired-end reads): 150

JSON report: fastp.json
HTML report: fastp.html

fastp --thread 8 -i SRR10971381_1.fastq -I SRR10971381_2.fastq -o SRR10971381_1_fastp.fastq -O SRR10971381_2_fastp.fastq --thread 8
fastp v0.20.0, time used: 544 seconds
```

FastQC.

```bash
mkdir fastqc_out
fastqc -o fastqc_out -f fastq SRR10971381_1_fastp.fastq SRR10971381_2_fastp.fastq
```

Use BWA to map to MN908947.

```bash
mkdir bwa_index
cp MN908947.fa bwa_index
cd bwa_index
bwa index MN908947.fa

bwa mem -t 8 raw/bwa_index/MN908947.fa raw/SRR10971381/SRR10971381_1_fastp.fastq raw/SRR10971381/SRR10971381_2_fastp.fastq | samtools sort - -o result/SRR10971381_MN908947.bam
```

Stats.

```bash
samtools flagstat -@8 SRR10971381_MN908947.bam 
26137357 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
28875 + 0 supplementary
0 + 0 duplicates
174256 + 0 mapped (0.67% : N/A)
26108482 + 0 paired in sequencing
13054241 + 0 read1
13054241 + 0 read2
136034 + 0 properly paired (0.52% : N/A)
136414 + 0 with itself and mate mapped
8967 + 0 singletons (0.03% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)

samtools view -F 0x804 -f 2 -b SRR10971381_MN908947.bam > SRR10971381_MN908947_mapped.bam
```

[Variant calling](https://samtools.github.io/bcftools/howtos/variant-calling.html).

```bash
bcftools mpileup -f raw/MN908947.fa result/SRR10971381_MN908947_mapped.bam | bcftools call -mv -Ov -o result/SRR10971381_MN908947_mapped.vcf
```

## Links

* https://github.com/CSSEGISandData/COVID-19
* https://github.com/galaxyproject/SARS-CoV-2
* https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/
* https://www.sciencemag.org/news/2020/01/mining-coronavirus-genomes-clues-outbreak-s-origins
* https://www.cdc.gov/coronavirus/2019-ncov/lab/rt-pcr-panel-primer-probes.html
* https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6356540/

