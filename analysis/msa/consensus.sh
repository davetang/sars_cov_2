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
   ../../script/unique_seq.pl -f ${variant}/genomic.fna > ${tmpdir}/genomic.uniq.fa 2> ${tmpdir}/seq_num.txt
   mv ${tmpdir}/genomic.uniq.fa ${tmpdir}/seq_num.txt ${variant}
else
   >&2 echo ${variant}/genomic.uniq.fa already exists
fi

numseq=$(cat ${variant}/seq_num.txt | grep "Unique" | cut -f2 -d' ')

if [[ ${numseq} -gt 50000 ]]; then
   >&2 echo ${variant}/genomic.uniq.fa contains ${numseq} sequences and will be subsampled to 50,000 sequences
   seqtk sample -s1984 ${variant}/genomic.uniq.fa 50000 > ${tmpdir}/genomic.uniq.fa
   mv -f ${tmpdir}/genomic.uniq.fa ${variant}
fi

if [[ ! -e ${variant}/genomic.uniq.mafft.fa ]]; then
   >&2 echo Generating ${variant}/genomic.uniq.mafft.fa

   if [[ ${numseq} -lt 2000 ]]; then
      mafft --thread ${threads} --retree 2 --inputorder ${variant}/genomic.uniq.fa > ${tmpdir}/genomic.uniq.mafft.fa
   elif [[ ${numseq} -lt 10000 ]]; then
      # FFT-NS-1 (very fast; recommended for >2000 sequences; progressive method with a rough guide tree):
      mafft --thread ${threads} --inputorder --retree 1 --maxiterate 0 ${variant}/genomic.uniq.fa > ${tmpdir}/genomic.uniq.mafft.fa
   elif [[ ${numseq} -le 50000 ]]; then
      # NW-NS-PartTree-1 (recommended for ~10,000 to ~50,000 sequences; progressive method with the PartTree algorithm):
      mafft --thread ${threads} --inputorder --retree 1 --maxiterate 0 --nofft --parttree ${variant}/genomic.uniq.fa > ${tmpdir}/genomic.uniq.mafft.fa
   else
      >&2 echo ${variant}/genomic.uniq.fa contains more than 50,000 sequences!
      exit 1
   fi

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

rmdir ${tmpdir}

>&2 echo Done!

