#!/usr/bin/env bash
#
# Download SARS-CoV-2 variants
#

os=$(uname -s)
arc=$(arch)
today=$(date +%Y%m%d)
path=$(dirname $0)

if [[ ${os} == Darwin ]]; then
   datasets=${path}/../bin/macos/datasets
elif [[ ${os} == Linux ]]; then
   if [[ ${arc} == x86_64 ]]; then
      datasets=${path}/../bin/ubuntu/datasets
   elif [[ ${arc} =~ arm ]]; then
      datasets=${path}/../bin/arm32/datasets
   fi
else
   >&2 echo Unrecognised operating system
   exit 1
fi

>&2 echo Using ${datasets}

lineages=(
BA.1
BA.2
BA.3
B.1.1.529
AY.1
AY.2
AY.3
B.1.351
B.1.351.2
B.1.351.3
P.1
P.1.1
P.1.2
B.1.617.2
B.1.1.7
)

if [[ -d ${path}/../raw/variants ]]; then
   >&2 echo ${path}/../raw/variants exists
else
   mkdir ${path}/../raw/variants
fi

for lineage in ${lineages[@]}; do

   ${datasets} download virus genome taxon SARS-CoV-2 --lineage ${lineage} --filename ${path}/../raw/variants/SARS-CoV-2-${lineage}.${today}.zip

   while [ $? -ne 0 ]
   do
      ${datasets} download virus genome taxon SARS-CoV-2 --lineage ${lineage} --filename ${path}/../raw/variants/SARS-CoV-2-${lineage}.${today}.zip
   done

done

>&2 echo Done!

exit 0

