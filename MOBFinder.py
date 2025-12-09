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
	usage = 'python3 %(prog)s [-h] [-i input_fasta] [-o output_dir]'
	main = 'Mobilizable Plasmid Relaxase Annotation'
	parser = argparse.ArgumentParser(usage = usage, description = main)
	parser.add_argument('-i', type = str, metavar = '', default = '', help = "input assembly fasta fragments")
	parser.add_argument('-o', type = str, metavar = '', default = '', help = "output directory")
	parser.add_argument('-v', action="version", version = "%(prog)s version 1.0", help = "show script's version and exit")
	parameter = parser.parse_args()
	return parameter.i, parameter.o

def dir_check(check_dir):
	check_dir = check_dir + '/' if not check_dir.endswith('/') else check_dir
	return check_dir

def dir_create(create_dir):
	if not os.path.exists(create_dir):
		os.mkdir(create_dir)

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

def fasta_vector_convert(vector, fasta, fasta_name, batch_list, outputdir):
	vectors = vector_read(vector)
	output_list = open(batch_list, 'w+')
	number = 0
	with open(fasta) as fastas:
		for line in fastas:
			if '>' in line:
				line = line.strip()
				if number != 0:
					index1 = int(number / 10000)
					index2 = number % 10000
					if index2 == 1:
						if index1 != 0:
							output.close()
						batch_id = "{}{}_{}.tsv".format(outputdir, fasta_name, index1)
						batch_result = "{}{}_{}.result".format(outputdir, fasta_name, index1)
						output_list.write("{}\t{}\n".format(batch_id, batch_result))
						output = open(batch_id, 'w+')
					vector_calculate(contig_id, contig, vectors, output)
				contig_id = line.split()[0][1:]
				contig = ''
				number += 1
			else:
				contig += line.strip()
		if number == 1:
			batch_id = "{}{}_1.tsv".format(outputdir, fasta_name)
			batch_result = "{}{}_1.result".format(outputdir, fasta_name)
			output_list.write("{}\t{}\n".format(batch_id, batch_result))
			output = open(batch_id, 'w+')
		vector_calculate(contig_id, contig, vectors, output)
	output.close()
	output_list.close()

def batch_remove(batch_list, output):
	batchlist = open(batch_list, 'r')
	index = 1
	for line in batchlist:
		line = line.strip().split()
		if index == 1:
			os.system(("cat {} >> {}".format(line[1], output)))
		else:
			os.system(("sed '1d' {} >> {} ".format(line[1], output)))
		index += 1
		os.system(("rm {}".format(line[0])))
		os.system(("rm {}".format(line[1])))
	batchlist.close()
	os.system(("rm {}".format(batch_list)))

def mobfinder_predict(vector, fasta, outputdir):
	outputdir = dir_check(outputdir)
	dir_create(outputdir)
	fasta_name = filename_extract(fasta)
	batch_list = "{}batch_list.tsv".format(outputdir)
	fasta_vector_convert(vector, fasta, fasta_name, batch_list, outputdir)
	result = "{}{}_mobfinder_result.tsv".format(outputdir, fasta_name)
	os.system(("Rscript MOBFinder.R {}".format(batch_list)))
	batch_remove(batch_list, result)

def main():
	input_fasta, output_dir = parameter_passing()
	input_vector = './MOBFinder_vector.w2v'
	mobfinder_predict(input_vector, input_fasta, output_dir)

if __name__ == '__main__':
	main()
