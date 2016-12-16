import os
import sys
import numpy as np


def scan_cpg(wigfile, wigfile_index,  wigfile_num, cpg_dict, cpg_ml_dict):
	for line in open(wigfile):
		cols = line.rstrip().split()
		if cols[0] == "track":
			continue
 		elif cols[0] == "variableStep":
			chrom = cols[1].split("=")[1]
		else:
			coord = cols[0]	
			value = cols[1]
			if not cpg_dict.has_key(coord):
				cpg_dict[coord] = [0] * wigfile_num
				cpg_ml_dict[coord] = ["nan"] * wigfile_num
			cpg_ml_dict[coord][wigfile_index] = float(value)
			if float(value) <= 0.1:
				cpg_dict[coord][wigfile_index] = 1
			else:
				cpg_dict[coord][wigfile_index] = -1
	return cpg_dict,cpg_ml_dict,chrom

def detect_umc(cpg_dict, cpg_ml_dict, wigfile_num, wigfile_label, name, chrom):
	'''
	output the UMC; discard CpG > 0.1 
	get the UMC occuring sample numbers
	'''
	umc_occurance = {}
	umc_file_name = name + "_all_umc.txt"	
	umc_ml_file_name = name + "_all_umc_mratio.txt"	
	umc_file_out = open(umc_file_name, "w")
	umc_ml_file_out = open(umc_ml_file_name, "w")
	umc_file_out.write("chr\tstart\tend\t#UMC\t#noUMC\t#sample") 
	umc_ml_file_out.write("chr\tstart\tend") 
	for l in wigfile_label:
		umc_file_out.write("\t%s" %l)
		umc_ml_file_out.write("\t%s" %l)
	umc_file_out.write("\n")
	umc_ml_file_out.write("\n")
	for coord in cpg_dict.keys():
		if 1 in cpg_dict[coord]:
			umc_num =  len([i for i,val in enumerate(cpg_dict[coord]) if val == 1])
			umc_nonum =  len([i for i,val in enumerate(cpg_dict[coord]) if val == -1])
			umc_file_out.write("%s\t%d\t%d\t%d\t%d\t%d" %(chrom, int(coord), int(coord)+1, umc_num, umc_nonum,wigfile_num))
			umc_ml_file_out.write("%s\t%d\t%d" %(chrom, int(coord), int(coord)+1 ))
			for i in range(0,len(cpg_dict[coord])):
				umc_file_out.write("\t%d" %cpg_dict[coord][i])
				umc_ml_file_out.write("\t%s" %str(cpg_ml_dict[coord][i]))
			umc_file_out.write("\n")
			umc_ml_file_out.write("\n")
			umc_occurance[coord] = umc_num
	del cpg_dict
	del cpg_ml_dict
	umc_file_out.close()
	umc_ml_file_out.close()
	return umc_occurance

def main():

	'''
	Input: name of outfile; list of wig.file (samne chr wig.file); 
	Output: 1. UMC with occupancy in 31 WGBS (row: UMC coordilation; col: chr; start; end; #UMC; #noUMC; #sample; UMC state in 31 WGBS [value: <=0.1 -> 1 | >0.1 -> -1 | nan -> 0])
			#2. Conserved UMR bed file (coordilation, name,conserved score, strand, average meth level, #CpG, including CpGs; 	
	'''
	parameters_num = len(sys.argv)
	name = sys.argv[1]
	wigfile_list = sys.argv[2: (2+(parameters_num-2)/2)]
	wigfile_label = sys.argv[(2+(parameters_num-2)/2):parameters_num]
	wigfile_num = len(wigfile_list)
	cpg_dict = {} ###key: coordilation value; list:  1 -> UMC,-1 -> noUMC, 0 -> noCovered)
	cpg_ml_dict = {} ###key: coordilation value; list:  methylation level)
	for w in range(wigfile_num):
		print("processing %s..." %wigfile_list[w])
		cpg_dict,cpg_ml_dict, chrom = scan_cpg(wigfile_list[w], w, wigfile_num, cpg_dict, cpg_ml_dict)	
	umc_occurance = detect_umc(cpg_dict, cpg_ml_dict, wigfile_num, wigfile_label, name, chrom)
	print np.mean(umc_occurance.values())
	print np.std(umc_occurance.values())

if __name__ == '__main__':
	main()
