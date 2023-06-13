#!/usr/bin/python3
#---------------------------------------
# Mobilizable Plasmid Classification
# Written by FengTao, based on blastx
#---------------------------------------

# Usage: python3 mob_classification.py -h

from Bio import SeqIO
import argparse
import re
import os
import math
import numpy as np

def parameter_passing():
	usage = 'python3 %(prog)s [-h] [-d relaxase_db] [-i input_fasta] [-g input_genebank] [-e evalue_cutoff] [-o output_dir]'
	main = 'Mobilizable Plasmid Relaxase Annotation'
	parser = argparse.ArgumentParser(usage = usage, description = main)
	parser.add_argument('-d', type = str, metavar = '', default = '', help = "local path of relaxase blast database")
	parser.add_argument('-i', type = str, metavar = '', default = '', help = "input plasmid fasta file")
	parser.add_argument('-g', type = str, metavar = '', default = '', help = "input plasmid genebank file")
	parser.add_argument('-e', type = str, metavar = '', default = '', help = "blast evalue cutoff")
	parser.add_argument('-o', type = str, metavar = '', default = '', help = "blast output directory")
	parser.add_argument('-v', action="version", version = "%(prog)s version 1.0", help = "show script's version and exit")
	parameter = parser.parse_args()
	return parameter.d, parameter.i, parameter.g, parameter.e, parameter.o


def dir_check(input_dir):
	input_dir = input_dir + '/' if not input_dir.endswith('/') else input_dir
	return input_dir


def feature_check(check_tag, check_dic):
	check_value = check_dic.get(check_tag, [''])[0] if check_tag in check_dic else 'na'
	return check_value


def genebank_cds_extract(genebank, input_dir):
	for gb_record in SeqIO.parse(genebank, 'gb'):
		gb_id = gb_record.id
		gb_output = open("{}{}_cds.fasta".format(input_dir, gb_id), 'a')
		cd_num = 0
		for feature in gb_record.features:
			if feature.type == "CDS":
				cd_num += 1
				if "parts" in feature.location.__dict__:
					location = ''
					for num in range(len(feature.location.parts)):
						location_parts = "[{},{}]_({})".format(feature.location.parts[num].start, feature.location.parts[num].end, feature.location.parts[num].strand)
						location = "{}-{}".format(location, location_parts) if location else location_parts
				else:
					location = "[{},{}]_({})".format(feature.location.start, feature.location.end, feature.location.strand)
				gene = feature_check('gene', feature.qualifiers)
				note = feature_check('note', feature.qualifiers)
				product = feature_check('product', feature.qualifiers)
				protein_id = feature_check('protein_id', feature.qualifiers)
				translation = feature_check('translation', feature.qualifiers)
				if "translation" in feature.qualifiers:
					gb_output.write(">{}_cds_{}_{}|".format(gb_id, cd_num, protein_id))
					gb_output.write("{} | ".format(location))
					gb_output.write("{} | ".format(gene))
					gb_output.write("{} | ".format(product))
					gb_output.write("{}\n".format(note))
					gb_output.write("{}\n".format(translation[:-1]))
		gb_output.close()


def blast_extract(blast_output):
	if not os.path.getsize(blast_output):
		qseqid, location, sseqid, bitscore, evalue, pident, qcovhsp = '', '', '', 0, 0, 0, 0
	else:
		with open(blast_output) as blast:
			line = next(blast).strip().split()

			qseqid, location, sseqid = line[0].split('|')[0], line[0].split('|')[1], line[4]
			bitscore, evalue, pident, qcovhsp = float(line[8]), float(line[9]), float(line[15]), float(line[18])

	return qseqid, location, sseqid, bitscore, evalue, pident, qcovhsp


def relaxase_search(db_dir, fasta, genebank, evalue, output_folder):
	"""blastp -query input.fasta -db ./4.relaxase_db/MOBB/MOBB_DB -out mobf_blastx_results.tsv -outfmt 6"""
	"""qseqid qlen qstart qend sseqid slen sstart send bitscore evalue length nident positive mismatch gaps pident ppos qcovs qcovhsp"""
	relaxase_list = ['MOBB', 'MOBC', 'MOBF', 'MOBH', 'MOBL', 'MOBM', 'MOBP', 'MOBQ', 'MOBT', 'MOBV']
	gb_id = fasta.split('/')[-1] if re.search('/', fasta) else fasta
	gb_id = gb_id.split('.fasta')[0]

	db_dir = dir_check(db_dir)
	output_folder = dir_check(output_folder)
	mob_output = open("{}{}_mob.tsv".format(output_folder, gb_id), 'a')
	mob_output.write("ID\tCategory\tProteinID\tLocation\tMOBscore\tConfidence\t")
	mob_output.write('MOBB\tMOBC\tMOBF\tMOBH\tMOBL\tMOBM\tMOBP\tMOBQ\tMOBT\tMOBV\n')

	qseqid_max, location_max, sseqid_max, bitscore_max, evalue_max, pident_max, qcovhsp_max = '', '', '', 0, 0, 0, 0
	mob_score, plasmid_category, confidence, relaxase_record = 0, '', '', []

	genebank_cds_extract(genebank, output_folder)

	for relaxase in relaxase_list:
		blast_outfmt = "-outfmt '6 qseqid qlen qstart qend sseqid slen sstart send bitscore evalue length nident positive mismatch gaps pident ppos qcovs qcovhsp'"
		blast_query = "{0}{1}_cds.fasta".format(output_folder, gb_id)
		blast_db = "{0}{1}/{1}_DB".format(db_dir, relaxase)
		blast_output = "{0}{1}_{2}_blast.tsv".format(output_folder, gb_id, relaxase)
		blast_result = "{0}{1}_{2}_blast.result".format(output_folder, gb_id, relaxase)
		blast_search = "blastp -query {0} -db {1} -evalue {2} -out {3} {4}".format(blast_query, blast_db, evalue, blast_output, blast_outfmt)
		blast_sort = "sort -n -r -k 9 {0} >> {1}".format(blast_output, blast_result)

		os.system(blast_search)
		os.system(blast_sort)
		os.system(("rm {0}".format(blast_output)))

		mob_qseqid, mob_location, mob_sseqid, mob_bitscore, mob_evalue, mob_pident, mob_qcovhsp = blast_extract(blast_result)

		if mob_bitscore > bitscore_max or mob_bitscore == bitscore_max and mob_qcovhsp > qcovhsp_max:
			qseqid_max, location_max, sseqid_max, bitscore_max, evalue_max, pident_max, qcovhsp_max = blast_extract(blast_result)
			plasmid_category = relaxase

		if mob_qseqid:
			relaxase_record.append((mob_qseqid, mob_location, mob_sseqid, mob_bitscore, mob_evalue, mob_pident, mob_qcovhsp))
		else:
			relaxase_record.append(('na'))

		os.system(("rm {0}".format(blast_result)))

	os.system(("rm {0}".format(blast_query)))

	if qseqid_max:
		mob_score = round(math.sqrt(0.01 * qcovhsp_max * (1 - 1 / np.log10(float(bitscore_max)))), 4)
		if mob_score >= 0.5 and evalue_max < float(1e-10):
			confidence = 'sure'
		else:
			confidence = 'possible'
	else:
		plasmid_category, confidence = 'non-mob', 'sure'

	mob_output.write("%s\t%s\t%s\t%s\t%s\t%s\t" % (gb_id, plasmid_category, qseqid_max.split('_')[-1], location_max, mob_score, confidence))
	mob_output.write("%s\t%s\t%s\t%s\t%s\t" % (relaxase_record[0], relaxase_record[1], relaxase_record[2], relaxase_record[3], relaxase_record[4]))
	mob_output.write("%s\t%s\t%s\t%s\t%s\n" % (relaxase_record[5], relaxase_record[6], relaxase_record[7], relaxase_record[8], relaxase_record[9]))

	mob_output.close()

def main():
	relaxasedb_dir, input_fasta, input_genebank, input_evalue, output_dir = parameter_passing()
	relaxase_search(relaxasedb_dir, input_fasta, input_genebank, input_evalue, output_dir)

if __name__ == '__main__':
	main()
