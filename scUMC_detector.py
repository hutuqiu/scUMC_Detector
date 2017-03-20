import os
import sys
import numpy as np

def get_umc_dict(umc_file):
	umc_dict = {}
	for line in open(umc_file):
		cols = line.rstrip().split()
		chrom = cols[0]
		coord = int(cols[1])
		sample_number = int(cols[3])
		if not umc_dict.has_key(chrom):
			umc_dict[chrom] = {}
		umc_dict[chrom][coord] = sample_number
	#coord_list = [] 
	value_list = []
	for chrom in umc_dict.keys():
		#coord_list = coord_list + umc_dict[chrom].keys()
		value_list = value_list + umc_dict[chrom].values()
	#coord_list.sort()
	#smean = np.mean(umc_dict.values())
	#ssd = np.std(umc_dict.values())
	smean = np.mean(value_list)
	ssd = np.std(value_list)
	return umc_dict,smean,ssd

def scan_con(umc_dict,  smean, ssd, h, sample_cutoff, outfile):
	out = open(outfile,'w')
	umc_con = {}
	if sample_cutoff != 0:
		h = (sample_cutoff - smean)/ssd
		print sample_cutoff, h
	else: 
		print h*ssd+smean, h
	for chrom in umc_dict.keys():
		for coord in umc_dict[chrom]:	
			hmax = (umc_dict[chrom][coord] - smean)/ssd
			if hmax < h:
				continue
			else:
				if not umc_con.has_key(chrom):
					umc_con[chrom] = {}
				umc_con[chrom][coord] = umc_dict[chrom][coord]
				start = coord
				end = coord
				out.write("%s\t%d\t%d\t%d\t%d\n" %(chrom, start-1, end, end-start+1, umc_dict[chrom][coord]))	
	out.close()
	return umc_con			


def main():
	'''
	Input: umc file
	Output: Conserved UMR bed file (coordilation, name,conserved score, strand, average meth level, #CpG, including CpGs; 	 
	'''	
	umc_file = sys.argv[1]
	pvalue_cutoff = float(sys.argv[2])
	name = sys.argv[3]
	if len(sys.argv) > 4:
		sample_cutoff = int(sys.argv[4])
	else:
		sample_cutoff = 0
	h=(1/pvalue_cutoff)**0.5
	#outfile = name + "_conserved_solo_UMC.bed"
	outfile = name
	umc_dict,smean,ssd = get_umc_dict(umc_file)
	print umc_file
	print pvalue_cutoff, h,  
	umc_con = scan_con(umc_dict,  smean, ssd, h, sample_cutoff,outfile)


if __name__ == '__main__':
	main()
