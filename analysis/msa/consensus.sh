#!/usr/bin/env bash

set -euo pipefail

variant=P.1
threads=8

rand=$$$RANDOM
tmpdir=/tmp/${rand}

if [[ -d ${tmpdir} ]]; then
   >&2 echo ${tmpdir} already exists!
   exit 1
else
   mkdir ${tmpdir}
fi

trap "rm -rf ${tmpdir}" SIGINT SIGTERM

if [[ ! -d ${variant} ]]; then
   mkdir ${variant}
else
   >&2 echo Directory ${variant} already exists
fi

if [[ ! -e ${variant}/genomic.fna ]]; then
   >&2 echo Unzipping ../../raw/variants/SARS-CoV-2-${variant}.20210819.zip
   unzip -p ../../raw/variants/SARS-CoV-2-${variant}.20210819.zip ncbi_dataset/data/genomic.fna > ${tmpdir}/genomic.fna
   mv ${tmpdir}/genomic.fna ${variant}
else
   >&2 echo ${variant}/genomic.fna already exists
fi

if [[ ! -e ${variant}/genomic.uniq.fa ]]; then
   >&2 echo Generating ${variant}/genomic.uniq.fa
   ../../script/unique_seq.pl -f ${variant}/genomic.fna > ${tmpdir}/genomic.uniq.fa
   mv ${tmpdir}/genomic.uniq.fa ${variant}
else
   >&2 echo ${variant}/genomic.uniq.fa already exists
fi

if [[ ! -e ${variant}/genomic.uniq.mafft.fa ]]; then
   >&2 echo Generating ${variant}/genomic.uniq.mafft.fa

   # FFT-NS-1 (very fast; recommended for >2000 sequences; progressive method with a rough guide tree):
   # mafft --retree 1 --maxiterate 0 input [> output]

   mafft --thread ${threads} --retree 1 --maxiterate 0 --inputorder ${variant}/genomic.uniq.fa > ${tmpdir}/genomic.uniq.mafft.fa
   mv ${tmpdir}/genomic.uniq.mafft.fa ${variant}
else
   >&2 echo ${variant}/genomic.uniq.mafft.fa already exists
fi

if [[ ! -e ${variant}/consensus.fa ]]; then
   >&2 echo Generating ${variant}/consensus.fa
   cons -sequence ${variant}/genomic.uniq.mafft.fa -outseq ${tmpdir}/consensus.fa
   mv ${tmpdir}/consensus.fa ${variant}
else
   >&2 echo ${variant}/consensus.fa already exists
fi

>&2 echo Done!

