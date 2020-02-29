# TreeMer
A simple tool to generate hierarchical clustering trees from nucleotide sequences using kmer spectra distance. 
Included is a large testset of SARSCOV2 genomes downloaded from https://www.nlm.nih.gov/news/coronavirus_genbank.html.

## Overview

This tool calculates the distance between a set of nucleotide sequences in FASTA format by digesting them into kmer count 
vectors (effectively kmer spectra). The pairwise distance between all pairs of vectors are calculated and clustered to build a
Hierarchical clustering tree. A number of distance metrics and clustering methods are supported (see supported methods section
below).




![SARSCOV2 Hierarchical Clustering (ward, euclidean)](https://github.com/ArthurVM/TreeMer/blob/master/SARSCOV2/HC_ed_ward.png)
