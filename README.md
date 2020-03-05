## README

* [Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)](https://en.wikipedia.org/wiki/Severe_acute_respiratory_syndrome_coronavirus_2)

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
esearch -db nuccore -query coronavirus |
efetch -db sequences -format fasta > raw/coronavirus_20200301.fa

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
makeblastdb -dbtype nucl \
            -in ../raw/coronavirus_20200301.fa \
            -input_type fasta \
            -title coronavirus_20200301 \
            -out coronavirus_20200301

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
```

## Links

* https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/
* https://www.sciencemag.org/news/2020/01/mining-coronavirus-genomes-clues-outbreak-s-origins
* https://www.cdc.gov/coronavirus/2019-ncov/lab/rt-pcr-panel-primer-probes.html
* https://github.com/CSSEGISandData/COVID-19
* https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6356540/

