#!/bin/bash

# Parameter file to run PrimerPhasing
# File generated 2020-1-02

python3 /mnt/hgfs/OneDrive_UNC/Projects/Programs/NGS_PrimerPhasing/PrimerPhasing.py --options_file /mnt/hgfs/OneDrive_UNC/Projects/Programs/NGS_PrimerPhasing/run_PrimerPhasing_Template.sh

exit

--Build_PhiX_DataFrame	False
--PickleFileFolder	/mnt/hgfs/OneDrive_UNC/Projects/Programs/NGS_PrimerPhasing/

--WorkingFolder	/mnt/hgfs/Drive_D/Testing/
--RefSeq	/mnt/hgfs/OneDrive/Bioinformatics/RefSeq/GRCh38/GRCh38.p12.fa.bgz
# The amplicon can be genomic coordinates that use RefSeq or a string.
--Amplicon	chr1:225423928-225424162
# --Amplicon	ACGACCTGGTGAACACCTAGGACGCACCATTCTCACAAAGGGAGTTTTCCACACGGACACCCCCCTCCTCACCACAGCCCTGCCAGGACGGGGCTGGCTACTGGCCTTATCTCACAGGTAAAACTGACGCACGGAGGAACAATATAAATTGGGGACTAGAAAGGTGAAGAGCCAAAGTTAGAACTCAGGACCAACTTATTCTGATTTTGTTTTTCCAAACTGCTTCTCCTCTTGGGAAGTGTAAGGAAGCT
--SeqLength	150
--PairedEnd	True
# This is the sequence to insert to do the phasing.  Keep it <=8
--PhasingSeq	AGACTAAA

# % PhiX spike.  Must be integer 0 to 50
--PhiX_Fraction	30

--OutFileName	LBR2.1_30_PhiX
--ChartTitle	Phasing Test, LBR2.1; 30% PhiX
