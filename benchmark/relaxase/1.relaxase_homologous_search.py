#!/usr/bin/python3
#-------------------------------------------
# Generate Plasmid Relaxase database
# Written by FengTao, based on blast 2.12.0+
#-------------------------------------------

# Usage: python3 1.relaxase_homologous_db.py -r ./01.relaxase_index/ -d /home/Database/nr/nr_db/nr -t 20 -c 70 -i 70 -e 1e-10 -o .

import argparse
import os

def arg_parser_relaxasedb():

	mob_usage = 'python3 %(prog)s [-h] [-r relaxase_index_dir] [-d protein_db] [-t threads] [-c query_cov] [-i identity] [-e evalue] [-o output_dir]'
	mob_main = 'Mobilizable Plasmid Relaxase Database Generating'

	parser = argparse.ArgumentParser(usage = mob_usage, description = mob_main)
	parser.add_argument('-r', type = str, metavar = '', default = '', help = "input directory of relaxase index sequence")
	parser.add_argument('-d', type = str, metavar = '', default = '', help = "protein database path of homology blast")
	parser.add_argument('-t', type = str, metavar = '', default = '', help = "threads number of homology blast")
	parser.add_argument('-c', type = str, metavar = '', default = '', help = "query coverage cut-off of homology blast")
	parser.add_argument('-i', type = str, metavar = '', default = '', help = "identity cut-off of homology blast")
	parser.add_argument('-e', type = str, metavar = '', default = '', help = "evalue cut-off ofhomology blast")
	parser.add_argument('-o', type = str, metavar = '', default = '', help = "output directory of relaxase blast result")
	parser.add_argument('-v', action="version", version = "%(prog)s version 1.0", help = "show script's version and exit")
	parameter = parser.parse_args()

	return parameter.r, parameter.d, parameter.t, parameter.c, parameter.i, parameter.e, parameter.o


def blastp_homology_search(input_dir, input_db, input_threads, input_qcov, input_identity, input_evalue, output_dir):

	relaxase_list = ['MOBB', 'MOBC', 'MOBF', 'MOBH', 'MOBL', 'MOBM', 'MOBP', 'MOBQ', 'MOBT', 'MOBV']

	if not input_dir.endswith('/'):
		input_dir = input_dir + '/'

	if not output_dir.endswith('/'):
		output_dir = output_dir + '/'

	blastp_result = output_dir + '2.relaxase_blastp'

	os.mkdir(blastp_result)

	for relaxase in relaxase_list:

		print("%s homologous protein blast begins!" % (relaxase))

		relaxase_index = input_dir + relaxase + '.fasta'
		relaxase_blastp = blastp_result + '/' + relaxase + '_blastp.result'

		relaxase_blast = 'blastp -query ' + relaxase_index + \
				 ' -db ' + input_db + \
				 ' -num_threads ' + input_threads + \
				 ' -max_target_seqs 50000 -outfmt ' + "'6 qseqid qlen qstart qend sseqid slen sstart send length pident positive qcovhsp bitscore evalue'" + \
				 ' -out ' + relaxase_blastp

		os.system(relaxase_blast)

		homologous_id = blastp_result + '/' + relaxase + '_homologous_id.tsv'
		homologous_fasta = blastp_result + '/' + relaxase + '_homologous.fasta'

		output = open(homologous_id, 'a')

		with open(relaxase_blastp) as relaxase_blastp_result:
			for line in relaxase_blastp_result:
				line = line.strip().split()
				if bool(int(line[11]) >= int(input_qcov))==True & bool(float(line[9]) >= float(input_identity))==True & bool(float(line[13]) < float(input_evalue))==True:
					line = line[4].split('|')
					output.write("%s\n" % (line[1]))
				else:
					continue
		output.close()

		blastdbcmd_extract = 'blastdbcmd -db ' + input_db + \
				   ' -entry_batch ' + homologous_id + \
				   ' -dbtype prot -out ' + homologous_fasta

		os.system(blastdbcmd_extract)
		os.system(('rm ' + relaxase_blastp))
		os.system(('rm ' + homologous_id))


if __name__ == '__main__':
	blastp_input, blastp_db, blastp_threads, blastp_qcov, blastp_identity, blastp_evalue, blastp_dir = arg_parser_relaxasedb()
	blastp_homology_search(blastp_input, blastp_db, blastp_threads, blastp_qcov, blastp_identity, blastp_evalue, blastp_dir)
