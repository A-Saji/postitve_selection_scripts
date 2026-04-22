# postitve_selection_scripts
Scripts used for calculating positive selection in ortholg clusters of genes

Perl scripts written by Vivek Thakur.
Before using make sure you have installed, muscle or clustalw or hmmer for translatorX.
Create and store the scripts in a single directory for easy use. 

Download and install RaXML from "https://cme.h-its.org/exelixis/web/software/raxml/" and cite them when using it. 
"RAxML Version 8: A tool for Phylogenetic Analysis and Post-Analysis of Large Phylogenies". In Bioinformatics, 2014, open access.

Similarly DOWNLOAD TranslatorX from "http://translatorx.co.uk/", and cite it when using it 
Citation: Abascal F, Zardoya R, Telford MJ (2010)
TranslatorX: multiple alignment of nucleotide sequences guided by amino acid translations
Nucleic Acids Res. 38:W7-13.	 

Download PAML from "http://abacus.gene.ucl.ac.uk/software/paml.html" cite them if you're using it as
Ziheng Yang, PAML 4: Phylogenetic Analysis by Maximum Likelihood, Molecular Biology and Evolution, Volume 24, Issue 8, August 2007, Pages 1586–1591, https://doi.org/10.1093/molbev/msm088

You can use Convert_omega_codeml_to_tree.py to generate newick trees for PAML results. 

We can use MSA_positions.py to generate postion table for MSA, which can be used for running the scripts on plotforsubstitutionfrequency.R to generate motif substitution frequency plots. 
