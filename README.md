##1. Collect 31 WGBS data
The data information for all the WGBS data is in: https://www.dropbox.com/s/mdu5vxdmu5hr2r8/31WGBS_Sample_Information.xlsx?dl=0
	
##2. WGBS data pre-processing 
  #2.1 Mapping with BSMAP 
  #2.2 QC with BSeQC
  2.3 Mlevel Calling with MOABS
  2.4 UMR (>=4CpG) detection
	
##3. Collect all the UMCs in 31 WGBS data
The file with all the UMCs in 31 WGBS data is in: https://www.dropbox.com/s/vl3xrkujhnpczxd/1.%2031WGBS_UMC_All.bed.zip?dl=0

##4. Remove UMRs (>=4CpG; detected in 31 WGBS) and SNP (SNP142 downloaded from UCSC), then get the Candidate UMC
	
	
##5. Integrate the sample occupancy for the Candidate UMC
   Input: the file name of outfile; list of Mlevel .wig file (samne chr .wig file);
   Output:  UMC with occupancy in 31 WGBS 
   		row: UMC coordilation; 
		col: chr; start; end; #UMC; #noUMC; #sample; UMC state in 31 WGBS [value: <=0.1 -> 1 | >0.1 -> -1 | nan -> 0]
   Example: 
   python ~/script/intergrate_UMC.py 31WGBS_UMC_Occupancy_chr1.txt ../WGBS_wig/cellwig/CD14_cells/chr1.wig ../WGBS_wig/cellwig/CD19cells/chr1.wig (other 29 .wig files) CD14_cells CD19cells (the names of other 29 WGBS samples)
	
##6. Using the "scUMC_detector.py" to calculate the Conservation Score and detect the sparse conserved under-methylated CpG
   ##Input: Candidate UMC with occupancy in 31 WGBS (https://www.dropbox.com/s/1ost397bpui6nun/2.%2031WGBS_UMC_Candidate_For_scUMC.bed.zip?dl=0)
   ##Output: scUMC detected in 31 WGBS based on 0.01 P-value in Chebyshev’s inequality method
	python scUMC_detector.py 31WGBS_UMC_Candidate_For_scUMC.bed 0.01 31WGBS_UMC_scUMC_Pvalue0.01.txt
