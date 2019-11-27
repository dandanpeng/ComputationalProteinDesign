# Computational Protein Design
Genetic Algorithm and Cuckoo Seach Algorithm are implemented here to work on protein design.

## Basic Idea
TERMs (TERtiary Motifs) are compact modules that recur in protein sturctures. In Mackensie et al. [1], they show that the protein structural universe is highly degenerate ("quantized" so to speak) at the level of TERMs, and a surprisingly small number of such motifs describes the majority of known structure space.

Based on this, we broke a protein structure into several TERMs and found structural matches for each TERM with MASTER[2] and TERMANAL[3]. Then we implemented genetic algorithm and cuckoo search algorithm to find the best combination of these structural matches.

## Enery Function
We wanted the structural matches' overlap parts' have the same amino acids and use as least fragments as possible to cover the whole protein sequence.
The energy function was defined as:

- First term: sum up the blosum score for all positions on the sequence
- Second term: regularization term to limit the number of fragments used

## Comparison to others
![](https://github.com/dandanpeng/ComputationalProteinDesign/blob/master/comparison.png)
