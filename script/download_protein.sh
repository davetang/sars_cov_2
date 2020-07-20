#!/usr/bin/env bash
#
# Download SARS-CoV-2 protein sequences; run script from this directory because of relative paths
#

proteins=(
orf1ab
orf1a
nsp1
nsp2
nsp3
nsp4
nsp5
nsp6
nsp7
nsp8
nsp9
nsp10
RdRp
nsp11
nsp13
nsp14
nsp15
nsp16
S
ORF3a
E
M
ORF6
ORF7a
ORF7b
ORF8
N
ORF10
)

mkdir -p ../raw/protein/
today=$(date +%Y%m%d)

for protein in ${proteins[*]}; do
   echo $protein;
   ../bin/datasets download virus protein $protein --filename ../raw/protein/SARS2.${today}.${protein}.zip
done

exit 0

