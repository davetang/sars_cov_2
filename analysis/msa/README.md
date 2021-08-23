## README

The script `consensus.sh` will extract the genomic sequences for a lineage, create a new FASTA with only unique sequences, produce a multiple sequence alignment (MSA) with MAFFT, and finally output a consensus sequence from the MSA using Cons.

MAFFT [speed-oriented methods](https://mafft.cbrc.jp/alignment/software/manual/manual.html):

* FFT-NS-i (iterative refinement method; two cycles only):

    mafft --retree 2 --maxiterate 2 input [> output]
    fftnsi input [> output]

* FFT-NS-i (iterative refinement method; max. 1000 iterations):

    mafft --retree 2 --maxiterate 1000 input [> output]

* FFT-NS-2 (fast; progressive method):

    mafft --retree 2 --maxiterate 0 input [> output]
    fftns input [> output]

* FFT-NS-1 (very fast; recommended for >2000 sequences; progressive method with a rough guide tree):

    mafft --retree 1 --maxiterate 0 input [> output]

* NW-NS-i (iterative refinement method without FFT approximation; two cycles only):

    mafft --retree 2 --maxiterate 2 --nofft input [> output]
    nwnsi input [> output]

* NW-NS-2 (fast; progressive method without the FFT approximation):

    mafft --retree 2 --maxiterate 0 --nofft input [> output]
    nwns input [> output]

*NW-NS-PartTree-1 (recommended for ~10,000 to ~50,000 sequences; progressive method with the PartTree algorithm):

    mafft --retree 1 --maxiterate 0 --nofft --parttree input [> output]

