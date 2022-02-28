#!/usr/bin/env bash
#
# Generates VCF file for a variant using NC_045512.2 as the reference
#

set -euo pipefail

usage() {
   >&2 echo "Usage: $0 [ -i accession_id ]"
   exit 1
}

num_param=1
required_param=$(bc -l<<<${num_param}*2+1)

while getopts ":i:" options; do
  case "${options}" in
    i)
      id=${OPTARG}
      ;;
    :)
      echo "Error: -${OPTARG} requires an argument."
      exit 1
      ;;
    *)
      usage ;;
  esac
done

if [[ ${OPTIND} -ne ${required_param} ]]; then
   usage
fi

check_depend(){
   tool=$1
   if [[ ! -x $(command -v ${tool}) ]]; then
     >&2 echo Could not find ${tool}
     exit 1
   fi
}

check_status(){
   step=$1
   if [[ $? -ne 0 ]]; then
     >&2 echo Something went wrong with step ${step}
     exit 1
   fi
}

dependencies=(java efetch snp-sites muscle)
for tool in ${dependencies[@]}; do
   check_depend ${tool}
done

path=$(dirname $0)/..
if [[ ! -e ${path}/raw/GCF_009858895.2_ASM985889v3_genomic.fna.gz ]]; then
   >&2 echo ${path}/raw/GCF_009858895.2_ASM985889v3_genomic.fna.gz does not exist
   exit 1
fi

if [[ ! -e ${path}/bin/snpEff/snpEff.jar ]]; then
   >&2 echo ${path}/bin/snpEff/snpEff.jar does not exist
   exit 1
fi

now(){
   date '+%Y/%m/%d %H:%M:%S'
}

SECONDS=0

>&2 printf "[ %s %s ] Fetching ${id}\n" $(now) 
efetch -db sequences -format fasta -id ${id} | gzip > ${path}/raw/${id}.fa.gz

>&2 printf "[ %s %s ] Combining FASTA files\n" $(now) 
gunzip -c ${path}/raw/GCF_009858895.2_ASM985889v3_genomic.fna.gz ${path}/raw/${id}.fa.gz | gzip > ${path}/raw/ref_vs_${id}.fa.gz

>&2 printf "[ %s %s ] Aligning with MUSCLE\n" $(now) 
gunzip -c ${path}/raw/ref_vs_${id}.fa.gz | muscle -out ${path}/raw/ref_vs_${id}.aln.fa

>&2 printf "[ %s %s ] Generating VCF\n" $(now) 
snp-sites -v -o ${path}/raw/ref_vs_${id}.vcf ${path}/raw/ref_vs_${id}.aln.fa

>&2 printf "[ %s %s ] Renaming contig\n" $(now) 
cat ${path}/raw/ref_vs_${id}.vcf | sed 's/ID=1/ID=NC_045512.2/;s/^1/NC_045512.2/' > ${path}/raw/blah
mv -f ${path}/raw/blah ${path}/raw/ref_vs_${id}.vcf

>&2 printf "[ %s %s ] Annotating variants\n" $(now) 
java -Xmx8g -jar ${path}/bin/snpEff/snpEff.jar -noStats NC_045512.2 ${path}/raw/ref_vs_${id}.vcf > ${path}/raw/ref_vs_${id}.ann.vcf

duration=$SECONDS
>&2 echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
>&2 echo Done
exit 0

