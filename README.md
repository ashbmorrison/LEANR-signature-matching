# LEANR-signature_matching
## Function to do signature matching using a modified version of the LEANR package and example results with public data

### Introduction
Signature matching is a common bioinformatic technique that quantitatively matches samples to gene sets that share a common biological annotation which can then provide molecular evidence of the association of biological pathways with diseases. We present LEANR+ORA as a potential method for signature matching using a graph networks approach. 

### Method Overview:
LEANR+ORA uses a modified version of the Local Enrichment ANalysis (LEAN) method, which implements pathway topology to find enriched protein subnetworks for signature matching (1). LEAN+ORA consists of three steps (Figure 1):
1. Differential expression analysis (DEA) is conducted on a transcriptomic/proteomic expression dataset. 
2. The full list of output markers and the associated p-values are overlayed onto a protein-protein interaction network (2) and the LEANR R package used to statistically calculate enrichment of local subnetworks (1). 
3. Significantly enriched subnetworks are individually tested for gene set enrichment using ORA (fora from FGSEA R package) (3). 

<p align="center">
  <img width="460" height="300" src="![image](https://github.com/ashbmorrison/LEANR-signature_matching/assets/16310308/0e7e80d3-a4c3-40c9-a927-d6dbcb158288)">
</p>


## References
1. Gwinner, Frederik, et al. “Network-Based Analysis of Omics Data: The LEAN Method | Bioinformatics | Oxford Academic.” OUP Academic, Oxford University Press, 6 Dec. 2016
2. Szklarczyk, Damian, et al. “STRING Database in 2023: Protein–Protein Association Networks and Functional Enrichment Analyses for Any Sequenced Genome of Interest | Nucleic Acids Research | Oxford Academic.” OUP Academic, Oxford University Press, 12 Nov. 2022
3. Korotkevich, Gennady, et al. “Fast Gene Set Enrichment Analysis | BioRxiv.” BioRxiv, bioRxiv



