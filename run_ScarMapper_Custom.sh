#!/bin/bash
#Parameter file to run ScarMapper
#File generated 2020-01-19

python3 /mnt/hgfs/OneDrive_UNC/Projects/Programs/ScarMapper/scarmapper.py --options_file /mnt/hgfs/Drive_D/ScarMapper_Data/run_ScarMapper_Custom.sh
exit

--IndelProcessing	True

--FASTQ1	/mnt/hgfs/Drive_D/ScarMapper_Data/RJ_JCG1_S1_L001_R1_001.fastq.gz
--FASTQ2	/mnt/hgfs/Drive_D/ScarMapper_Data/RJ_JCG1_S1_L001_R2_001.fastq.gz

--RefSeq	/mnt/hgfs/OneDrive/Bioinformatics/RefSeq/GRCm38/GRCm38.p6.fa.gz
--Master_Index_File	/mnt/hgfs/OneDrive/Bioinformatics/Indices/Master_Lists/Ramsden_Rosa26a_Oligos.csv
--SampleManifest	/mnt/hgfs/Drive_D/ScarMapper_Data/RamsdenSampleManifest.csv
--TargetFile	/mnt/hgfs/OneDrive_UNC/Projects/ScarMapper_Runs/ScarMapper_Targets.txt
--PrimerPhasingFile	/mnt/hgfs/OneDrive_UNC/Projects/ScarMapper_Runs/PrimerPhasing.txt
--WorkingFolder	/mnt/hgfs/Drive_D/ScarMapper_Data/


--Verbose	DEBUG
--Job_Name	Testing
--Spawn	5
--Demultiplex	False
--Species	Mouse
--Cell_Line	hTERT_RPE1
--Platform	Ramsden

--N_Limit	0.1
--Minimum_Length	100	# Length after trimming
--OutputRawData	False

