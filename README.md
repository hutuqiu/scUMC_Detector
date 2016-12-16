1. Collect 31 WGBS data
	The data information for all the WGBS data is in:                                                                         	  https://www.dropbox.com/s/mdu5vxdmu5hr2r8/31WGBS_Sample_Information.xlsx?dl=0
	
2. WGBS data pre-processing 
	2.1 Mapping with BSMAP 
	2.2 QC with BSeQC
	2.3 Mlevel Calling with MOABS
	2.4 UMR (>=4CpG) detection
	
3. Collect all the UMCs in 31 WGBS data
	The file with all the UMCs in 31 WGBS data is in:																			https://www.dropbox.com/s/vl3xrkujhnpczxd/1.%2031WGBS_UMC_All.bed.zip?dl=0

4. Remove UMRs (>=4CpG; detected in 31 WGBS) and SNP (SNP142 downloaded from UCSC), then get the Candidate UMC: 
	https://www.dropbox.com/s/1ost397bpui6nun/2.%2031WGBS_UMC_Candidate_For_scUMC.bed.zip?dl=0
	
5. Collect the sample occupancy for the Candidate UMC
	
6. Using the "scUMC_detector.py" to detect the sparse conserved under-methylated CpG
	python scUMC_detector.py C31WGBS_UMC_Candidate_For_scUMC.bed 0.01 31WGBS_UMC_scUMC_Fdr0.01.txt
