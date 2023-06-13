#!/usr/bin/python3
#-----------------------------------------------------#
# MOBFinder
# Random Select Training and Test Program
# Written by FengTao, in China (SMU)
#---------------------------------------

import argparse
import random
import os
import numpy as np

def parameter_passing():
	usage = 'python3 %(prog)s [-h] [-i input_fasta] [-b input_contig_bin] [-o output_dir]'
	main = 'Mobilizable Plasmid Relaxase Annotation'
	parser = argparse.ArgumentParser(usage = usage, description = main)
	parser.add_argument('-i', type = str, metavar = '', default = '', help = "input contig fasta")
	parser.add_argument('-b', type = str, metavar = '', default = '', help = "input contig binning result")
	parser.add_argument('-o', type = str, metavar = '', default = '', help = "output directory")
	parser.add_argument('-v', action="version", version = "%(prog)s version 1.0", help = "show script's version and exit")
	parameter = parser.parse_args()
	return parameter.i, parameter.b, parameter.o

def dir_check(check_dir):
	check_dir = check_dir + '/' if not check_dir.endswith('/') else check_dir
	if not os.path.exists(check_dir):
		os.mkdir(check_dir)
	return check_dir

def vector_read(vector):
	vector_dir = {}
	with open(vector, 'r') as vectors:
		for line in vectors:
			wordvec = line.strip().split()
			if len(wordvec) == 101:
				vector_dir[wordvec[0]] = [float(num) for num in wordvec[1:]]
			else:
				continue
	return vector_dir

def complement_chain(chain):
	base_pair = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
	rev_chain = ''
	for base in chain:
		if base in base_pair.keys():
			rev_chain = base_pair[base] + rev_chain
		else:
			rev_chain = base + rev_chain
	return rev_chain

def average_vector(vector_dic, chain):
	tag_num = 0
	total = np.array([0]*100)
	for num in range(0, len(chain)-3):
		tag = chain[num : num+4]
		if set(tag) <= set(['A', 'T', 'C', 'G']):
			tag_arr = np.array(vector_dic[tag])
			total = total + tag_arr
			tag_num += 1
		else:
			continue
	if tag_num != 0:
		total = total / tag_num
	return total

def quotient_remainder(num1, num2):
	num_1 = num1 // num2
	num_2 = num1 % num2
	return num_1, num_2

def vector_output(seq_id, seq_len, contig_len, contig_vec, output):
	output.write("{}\t".format(seq_id))
	output.write("{}\t".format(seq_len))
	output.write("{}\t".format(contig_len))
	for vec in range(0, len(contig_vec)-1):
		output.write("{}\t".format(contig_vec[vec]))
	output.write("{}\n".format(contig_vec[len(contig_vec)-1]))

def vector_calculate(seq_id, seq, vector_dic, output):
	quot, rema = quotient_remainder(len(seq), 1600)
	if quot != 0:
		for num in range(0, quot):
			seq_start = num * 1600
			seq_end = seq_start + 1600
			seq_tmp = seq[seq_start:seq_end]
			seq_tmp_rev = complement_chain(seq_tmp)
			seq_vec = average_vector(vector_dic, seq_tmp)
			seq_vec_rev = average_vector(vector_dic, seq_tmp_rev)
			seq_vec_ave = (seq_vec + seq_vec_rev) / 2
			vector_output(seq_id, len(seq), 1600, seq_vec_ave, output)
	if rema >= 4:
		seq_start = quot * 1600
		seq_end = len(seq)
		seq_tmp = seq[seq_start:seq_end]
		seq_tmp_rev = complement_chain(seq_tmp)
		seq_vec = average_vector(vector_dic, seq_tmp)
		seq_vec_rev = average_vector(vector_dic, seq_tmp_rev)
		seq_vec_ave = (seq_vec + seq_vec_rev) / 2
		vector_output(seq_id, len(seq), rema, seq_vec_ave, output)

def filename_extract(inputfile):
	inputfile = inputfile.split('/')[-1] if '/' in inputfile else inputfile
	filename = inputfile.split('.fa')[0]
	return(filename)

def contig_bin_extract(contig_bin):
	bin_result = open(contig_bin, 'r')
	bin_dic = {}
	for line in bin_result:
		line = line.strip().split(',')
		if line[1] in bin_dic.keys():
			bin_dic[line[1]].append(line[0])
		else:
			bin_dic[line[1]] = [line[0]]
	bin_result.close()
	return(bin_dic)

def contig_bin_vector_convert(vector, fasta, fasta_name, fasta_bin, bin_list, outputdir):
	vectors = vector_read(vector)
	output_list = open(bin_list, 'w+')
	for binid in fasta_bin.keys():
		binname = "{}{}_{}.tsv".format(outputdir, fasta_name, binid)
		output_list.write("{}\n".format(binname))
		output = open(binname, 'w+')
		number = 1
		fastaseq = open(fasta, 'r')
		for line in fastaseq:
			if '>' in line:
				line = line.strip()
				if number != 1:
					if contig_id in fasta_bin[binid]:
						vector_calculate(binid, contig, vectors, output)
				contig_id = line.split()[0][1:]
				contig = ''
				number += 1
			else:
				contig += line.strip()
		if contig_id in fasta_bin[binid]:
			vector_calculate(binid, contig, vectors, output)
		fastaseq.close()
		output.close()
	output_list.close()

def vector_remove(bin_list):
	with open(bin_list) as binlist:
		for line in binlist:
			os.system(("rm {}".format(line.strip())))
	os.system(("rm {}".format(bin_list)))

def mobfinder_predict(fasta, contig_bin, outputdir):
	outputdir = dir_check(outputdir)
	vector = './MOBFinder_vector.w2v'
	fasta_name = filename_extract(fasta)
	fasta_bin = contig_bin_extract(contig_bin)
	fasta_bin_list = "{}bin_list.tsv".format(outputdir)
	contig_bin_vector_convert(vector, fasta, fasta_name, fasta_bin, fasta_bin_list, outputdir)
	result = "{}{}_mobfinder_binning_result.tsv".format(outputdir, fasta_name)
	os.system(("Rscript MOBFinder.R {} {}".format(fasta_bin_list, result)))
	vector_remove(fasta_bin_list)

def main():
	input_fasta, input_fasta_bin, output_dir = parameter_passing()
	mobfinder_predict(input_fasta, input_fasta_bin, output_dir)

if __name__ == '__main__':
	main()
