#!/usr/bin/env bash

set -euo pipefail

os=$(uname -s)
arc=$(arch)
path=$(dirname $0)

if [[ ${os} == Darwin ]]; then
   mkdir -p ${path}/macos
   wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/mac/datasets -O ${path}/macos/datasets && chmod 755 ${path}/macos/datasets
   wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/mac/dataformat -O ${path}/macos/dataformat && chmod 755 ${path}/macos/dataformat
elif [[ ${os} == Linux ]]; then
   if [[ ${arc} == x86_64 ]]; then
      mkdir -p ${path}/ubuntu
      wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets -O ${path}/ubuntu/datasets && chmod 755 ${path}/ubuntu/datasets
      wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/dataformat -O ${path}/ubuntu/dataformat && chmod 755 ${path}/ubuntu/dataformat
   elif [[ ${arc} =~ arm ]]; then
      mkdir -p ${path}/arm32
      wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-arm/datasets -O ${path}/arm32/datasets && chmod 755 ${path}/arm32/datasets
      wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-arm/dataformat -O ${path}/arm32/dataformat && chmod 755 ${path}/arm32/dataformat
   fi
else
   >&2 echo Unrecognised operating system
   exit 1
fi

if [[ ! -d ${path}/snpEff ]]; then
   wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
   unzip snpEff_latest_core.zip
   rm snpEff_latest_core.zip
fi

>&2 echo Done
exit 0

