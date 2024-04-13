**BADGER**: **B**ivariate **A**nalysis of **D**ifferential **G**ene **E**xpression **R**eactions

**Authors**: Ellen I. Weinheimer, James B. Pease

**Description**: BADGER is used for detection of differentially reacting genes across multi-treatment gene expression time courses. The pipeline can be applied to time courses of any length or number of treatments, but requires at least two of each.

**Python Library Requirements**: sys, re, math, itertools, statistics, numpy, scipy, argparse

**Testing BADGER**: Run ./test.sh to ensure the script is working properly. 

**Basic Usage**: BADGER requires as input (1) a csv count table of normalized RNAseq reads (--input); (2) a csv file summarizing the experimental design factors with one sample per line and columns indicating treatment, timepoint, individual, and sample ID (--factors), and (3) the desired control condition label as seen in the treatment column of the experimental design factors file (--control). See files in test_files/ and command in ./test.sh as an example. 

**Questions?** Please contact Ellen Weinheimer (weine18@wfu.edu).
