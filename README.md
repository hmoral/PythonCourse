## PythonCourse

This repository contains python scripts for calculating basic population genetics statistics from pool-seq samples

# Genetic statistics: pi, Tajima's D, Fst and Dxy

# Details

The scripts take an allele frequency table containing likelihood probabilites and allele frequencies obtained with a model-based SNP caller (e.g. SNAPE or Freebayes).
Analysis are implemented in non-overlaping sliding windows
Inputs can be filtered from the frequencies table by missing data, likelihood probabilites and depth. Other filters such as allelic balance, strand bias, He, etc. should be implemented before in the vcf file

