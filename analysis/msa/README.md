## README

Align with MAFFT.

```bash
mafft --retree 2 --inputorder genomic.fna > genomic.mafft.fna
```

Create consensus with cons.

```bash
cons -sequence genomic.mafft.fna -outseq consensus.fa
```

