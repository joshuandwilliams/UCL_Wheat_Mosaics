# UCL_Wheat_Mosaics
Contains code from my summer internship at UCL with Professor Richard Mott

ShowSNPs.sh includes downloading of appropriate files (delta files and genome level assemblies), splitting them by chromosome, and calling SNPs (requiring data reformatting).

ForDynamic.R involves reformatting the output files from ShowSNPs.sh (.delta.snp files) into a form to be inputted into Prof Mott's dynamic programming algorithm. This includes subsetting non-matching positions, and removing indels and non-bialellic SNP.

The dynamic programming algorithm itself is not made available here since it was not my work.
