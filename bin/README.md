## README

[Seqtk](https://github.com/lh3/seqtk) is a fast and lightweight tool for processing sequences in the FASTA or FASTQ format.

```bash
# on macOS
cd macos

# on Ubuntu
cd ubuntu

wget https://github.com/lh3/seqtk/archive/refs/tags/v1.3.tar.gz
tar -xzf v1.3.tar.gz
cd seqtk-1.3
make

cd ..
mv seqtk-1.3/seqtk .
rm -rf seqtk-1.3 v1.3.tar.gz
```

Use [datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/quickstarts/command-line-tools/) to download biological sequence data across all domains of life from NCBI.


```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets -O ubuntu/datasets
wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/mac/datasets -O macos/datasets
chmod 755 macos/datasets
chmod 755 ubuntu/datasets
```

