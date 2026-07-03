# SequenceTools
Programs used in the lab for various tasks

getRBH.pl program for obtaining reciprocal best hits (RBH) using either of blastp, mmseqs, diamond, or lastal

calcANI.pl program for obtaining average nucleotide identity (ANI) by several methods.

getdNdS.pl program to obtain non-synonymous to synonymous substitutions with codeml (paml)

## Family/domain scans

scanFams.pl program for runnign domain scans using different programs. These assume you already have the domain databases built

cleanFams.pl this programs takes the output of scanFams.pl and cleans it by eliminating overlapped matches, or matches covering less than some minimal proportion of the domain model
