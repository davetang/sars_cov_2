#!/usr/bin/env bash
#
# Download SARS-CoV-2 variants using a Raspberry Pi
#

lineages=(
B.1.1.7
B.1.351
B.1.351.2
B.1.351.3
B.1.617.2
AY.1
AY.2
AY.3
P.1
P.1.1
P.1.2
)

today=$(date +%Y%m%d)
path=$(dirname $0)

if [[ -d ${path}/../raw/variants ]]; then
   >&2 echo ${path}/../raw/variants exists
else
   mkdir ${path}/../raw/variants
fi

for lineage in ${lineages[@]}; do
   bin/arm32/datasets download virus genome taxon SARS-CoV-2 --lineage ${lineage} --filename ${path}/../raw/variants/SARS-CoV-2-${lineage}.${today}.zip
done

>&2 echo Done!

exit 0

