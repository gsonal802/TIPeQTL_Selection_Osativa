# TIPeQTL_Selection_Osativa
Estimate and plot selection statistics

This repository contains code to perform the selection analysis on TIP-eQTLs, and plot the same

#Estimate Selection Statistics The selection script (getselstat.pl) estimates Gst between populations and nucleotide diversity (pi), and Tajima's D within the two populations. H12 was estimated using the script from Garud et al., 2012. Window size and slinding window size can be changed by editing the $window_size and $sliding variables. The input file is a tab separated Allele Count (AC) file with four columns -- Chromosome, Position, AC in Pop1, AC in Pop2 -- with header column present.

To run the script: perl getselstat.pl namePop1 namePop2 max-possible-AC-pop1 max-possible-AC-pop2 output-folder/ input-file

#Plot the results 
Rscript plotSel_haplotype.R
