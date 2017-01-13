#!/usr/bin/env python
from subprocess import call, Popen
#from Bio.Blast.Applications import NcbiblastnCommandline
import argparse
from glob import glob
import os
import collections
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import sys
from ete3 import Tree, ClusterTree
import random

def main ():

	'''
	Screen assemblies
	Can use conitgs (nuc), CDS or protein seqs
	Can use six frame contig translation if don't trust the CDS caller or you're looking for sub seqs
	If assemblies and queries are nucleotide sequence then blastn used.
	If 'assemblies' are accually proteins sequences and queries are protein sequence then blastp used.
	If assemblies are nucleotide sequence and query is protein then tblastn is used.
	'''

 	#Args
	parser = argparse.ArgumentParser(description='Screens assemblies')
	#parser.add_argument("-i", "--input", type=str, help="Directory with assemblies") 
	#(rax requires that the call be in the dir)
	parser.add_argument("-q", "--query", type=str, help="Query")
	parser.add_argument("-ph", "--phylogeny", type=str, help="Rax tree")
	parser.add_argument("-t", "--threads", type=str, help="How many threads to give blast, muscle, raxml etc", default = '4')
	parser.add_argument("-pi", "--percent_identiy", type=str, help="percent_identiy", default='80')
	parser.add_argument("-pl", "--percent_length", type=str, help="percent_length", default='80')
	parser.add_argument("-d", "--direct_str", type=bool, help="Use direct string comparison instead of blast - faster but only get hits at 100% len and sequence identity", default=False)
	parser.add_argument("-r", "--regex", type=str, help="regular expression; i.e., '.fa' will concatenate all files ending with .fa ")
	args = parser.parse_args()
	
	#run
	number_of_aassemblies = cat(args)
	blast_type = blast(args)
	csv(args)
	fasta(args)
	align(args, blast_type)
	gene_tree(args, blast_type)
	genome_tree(args)
	plot_vars(args)
	box(args, number_of_aassemblies)
	genome_tree(args)
	itol(args)

def randy_colors():
	
	'''
	Returns random colors for itol
	'''

	colors = set([])
	with open('/Users/lmcintyre/Documents/canvas/itol_templates/500hex_colors.txt','r') as fin:
		for i in fin:
   			if '#' in i and len(i.strip()) == 7:
   				colors.add(i.strip())
    	return random.sample(colors,552)

def convert_id():
	
	'''
	Removes # from sanger ids
	Oppitionally adds 'LD_' to id name (for regex)
	'''
	pass

def get_seq_type(seqs):
	
	#better with regex?
	amino_acids = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']	
	nucleotides = ['A','T','C','G','N','R','Y','S','W','K','M','B','D','H','V','N']
	# from http://www.bioinformatics.org/sms/iupac.html
	protein = False
	DNA = True
	for record in SeqIO.parse(seqs,'fasta'):
		for pos in str(record.seq).upper():
			if pos not in nucleotides:
				print pos
				DNA  = False
				protein = True
			elif pos not in nucleotides and pos not in amino_acids:
				print 'Error, not DNA or protein: ', pos, ' in ', record
				protein = False
				sys.exit(0)
		
	if DNA:
		seq_type = 'DNA'
	elif protein:
		seq_type = 'prot'
 		
	return seq_type


def cat(args):

	'''
	Concatenates assemblies; takes n assemblies and makes one file with all contigs
	'''
	
	fout = open('concatenated_assmblies.fasta', 'w')
	for i, assembly in enumerate(glob('*' + args.regex)):
		with open(assembly, 'r') as fin:
			for line in fin:
				if '>' in line: # to get all contigs with just sanger name and contig with all _
					if '.' in line: 
						line = '>' + '_'.join(line.strip().split('.')[1:]) + '\n'
				fout.write(line)
		fout.write('\n')
	fout.close()
	print str(i+1), ' assemblies concatenated.'
	return i+1


def call_blast(args, blast_type):

 
	subjects = 'concatenated_assmblies.fasta'
	#cant use qcov_hsp_perc with cline (HSP = High- scoring Segment Pair )
	cmd = [blast_type,
	'-query', args.query,
	'-db', subjects,
	'-out', 'blast_out.txt',
	'-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp',
	'-max_target_seqs', '500000',
	'-qcov_hsp_perc', args.percent_length,
	'-num_threads', args.threads]

	if blast_type == 'blastn':
		cmd += ['-perc_identity', args.percent_identiy]
	print cmd
	call(cmd)	


def blast(args):

	'''
	Blasts (Nuc or Prot) seq/s of interest (Query) against db of concatenated 
	assemblies(contigs)/CDS/proteins
	'''
 	print 'testing query seqs'	
	seq_type_query = get_seq_type(args.query)
 	print 'fine -> ', seq_type_query
	print 'testing db seqs'	
	seq_type_db = get_seq_type('concatenated_assmblies.fasta')
	print 'fine -> ', seq_type_db
	#db
	#index if not already done
	if not os.path.exists('concatenated_assmblies.fasta.nhr'):
		if seq_type_db == 'DNA':
			seq_type = 'nucl'
		else:
			seq_type = 'prot'
		try:
			call(['makeblastdb',
			'-in', 'concatenated_assmblies.fasta',
			'-parse_seqids',
			'-dbtype', seq_type])
		except:
			print 'you need to install makeblastdb or fix paths to it'

	#blast
	if seq_type_query == 'DNA' and seq_type_db == 'DNA':
		blast_type = 'blastn'
		call_blast(args, 'blastn')
		print 'Running blastn...'
	elif seq_type_query == 'prot' and seq_type_db == 'prot':
		blast_type = 'blastp'
		call_blast(args, 'blastp')
		print 'Running blastp...'
	elif seq_type_query == 'prot' and seq_type_db == 'DNA':
		blast_type = 'tblastn'
		call_blast(args, 'tblastn')
		print 'Running tblastn...'
	else:
		print "Error with sequence type, can't run blastn, blastp or tblastx"

	return blast_type

def direct_str():

	'''
	Exact match search - faster than blast, same output.
	Required for peptides < 10 aa long
	'''
	pass


def get_query_seqs(args):

	query_seqs = {} 
	for record in SeqIO.parse(args.query, 'fasta'):
		query_seqs[record.id] = record

	return query_seqs

def csv(args, binary = True):
	
	'''
	Makes to CSVs; a binary presence absence and a redundant total hits
	'''
	
	assemblies = []
	for assembly in glob('*' + args.regex):
		assemblies.append(assembly.strip().split('/')[-1])
	
	hits_dict, query_seqs = parse_blast(args, dict_key = 'assembly', dict_value = 'query')
	fout = open('hits.csv','w')
	header = '#Names,' + ','.join(query_seqs) +'\n'
	fout.write(header)
	for ass in hits_dict:
		fout.write(ass +',')
		for query in query_seqs:
			if query in hits_dict.get(ass):
				if binary:
					fout.write('1,')
				else:
					fout.write(str(hits_dict.get(ass).count(query))+',')
			else:
				fout.write('0,')
		fout.write('\n')
	fout.close()

def fasta(args):

	'''
	Makes a fasta file of all hits for each query
	'''
	query_seqs = get_query_seqs(args)
	for query in query_seqs:
		
		#make oonfig file from blast output
		current = query + '_config.txt'
		fout = open(current,'w')
		with open('blast_out.txt', 'r') as fin:
			for line in fin:
				bits = line.strip().split()
				if bits[0] == query:
					contig = bits[1] 
					seq_start = bits[8]
					seq_end = bits[9]
					if int(seq_start) > int(seq_end):#minus/plus is right, although seems backwards..
						line = ' '.join([contig, seq_end + '-' + seq_start, 'minus'])
						fout.write(line + '\n')
					else:
						line = ' '.join([contig, seq_start + '-' + seq_end, 'plus'])
						fout.write(line + '\n')
		fout.close()

		#Cut seq of interest out of contig
		call(['blastdbcmd',
			'-db', 'concatenated_assmblies.fasta',
			'-entry_batch', current,
			'-out', query + '_seqs.fasta'])

def align(args, blast_type):

	'''
	Aligns fastas with muscle
	'''
	query_seqs = get_query_seqs(args)
	for query in query_seqs:
		#add ref
		fout = open(query + '_seqs_and_ref.fasta','w')
		fout.write('>ref\n')
		record = query_seqs.get(query)
		fout.write(str(record.seq) + '\n')
		if blast_type == 'tblastn':#11 = bacterial
		        for record in SeqIO.parse(query + '_seqs.fasta', 'fasta'):
				record.seq = record.seq.translate(table = 11)
				SeqIO.write(record,fout,'fasta')
		else:
			with open(query + '_seqs.fasta','r') as fin:
				for line in fin:
					fout.write(line)
		fout.close()
		call(['muscle',
			'-in', query + '_seqs_and_ref.fasta',
			'-out', query + '_seqs.aln'])

def gene_tree(args, blast_type):

	'''
	Makes a Raxml tree of alignments
	'''
	#set model. need to make this accessible one day, not hard coded...
	if blast_type == 'blastn':
		model = 'GTRGAMMA'
	else:
		model = 'PROTGAMMAWAG'
	query_seqs = get_query_seqs(args)
	for query in query_seqs:
		cmd = ['raxmlHPC-PTHREADS',
			'-s', query + '_seqs.aln',
			'-n', query + '_seqs.tre',
			'-m', model,
			'-p', '12345',
			'-T', args.threads]
		call(cmd)

def genome_tree(args):

	'''
	Takes a tre file and anotes it with with gene variants
	'''

	percent_dict, query_seqs = parse_blast(args,dict_key = 'assembly', dict_value = 'percent')
	t = Tree(args.phylogeny)
	fout = open('percent.tsv','w')
	header = '#Names\t' + '\t'.join(query_seqs) +'\n'
	fout.write(header)
	for ass in percent_dict:
		fout.write(ass +'\t' + '\t'.join(percent_dict.get(ass)) + '\n')
	fout.close()

	t = ClusterTree(args.phylogeny, 'percent.tsv')
	t.render('heatmap.svg','heatmap')

def itol(args):

	'''
	Makes an itol colortrip template 
	'''
	color_dict = {}

	#Group by percent similarity
	percent_dict, query_seqs = parse_blast(args, dict_key = 'assembly', dict_value = 'percent')
	for i, query in enumerate(query_seqs):
		fout = open(query+'_itol.txt','w')
		randy = randy_colors()
		#header
		fout.write('DATASET_COLORSTRIP\n')
		fout.write('SEPARATOR SPACE\n')
		fout.write('DATASET_LABEL ' + query + '\n')
		fout.write('COLOR #ff0000\n')
		fout.write('DATA\n')
	
		for ass in percent_dict:
			#lazy round to int, need better way to divide up seqs, no of snps?
			percent = str(int(float(percent_dict.get(ass)[i])))
			if ass in color_dict:
				color = color_dict.get(ass)
			else:
				color = randy.pop()
				color_dict[ass] = color
			fout.write(' '.join([ass,color,percent]) + '\n')
	fout.close()


def plot_vars(args):

	'''
	Plots the distribution of vars (SNPs/indels) on the seq and saves as SVG
	'''
	
	query_seqs = get_query_seqs(args)
	for query in query_seqs:
		if not os.path.exists(query + '_seqs.aln'):
			print 'Skipping ' + query + '_seqs.aln'
			continue
		
		#ref
		ref_dik = collections.OrderedDict()
		
		for record in SeqIO.parse(query + '_seqs.aln', 'fasta'):
			if record.id == 'ref':
				ref_seq = str(record.seq)
				for i, nuc in enumerate(ref_seq):
					i+=1
					ref_dik[str(i)] = [nuc,0,0,0]#[ref_nuc, count_of_snps, count_of_ins, count_of_del]
		#hits
		number_hits = 0
		for record in SeqIO.parse(query + '_seqs.aln', 'fasta'):
			if record.id == 'ref':
				continue
			number_hits += 1 
			seq = str(record.seq)
			try:
				assert len(seq) == len(ref_seq)
			except:
				print 'alignments need to be same length!'
				print len(seq),  len(ref_seq)
			for i, nuc in enumerate(seq):
				i+=1
				if ref_dik.get(str(i))[0] == '-' and nuc == '-':
					continue
				elif ref_dik.get(str(i))[0] == '-':#ref is - == ins
					ref_dik.get(str(i))[2] += 1
				elif nuc == '-': #current seq is - == del
					ref_dik.get(str(i))[3] += 1
				else:
					if ref_dik.get(str(i))[0] != nuc:
						ref_dik.get(str(i))[1] += 1
		if number_hits == 0:
			continue
		for i,pos in enumerate(ref_dik):
			print query, ref_dik.get(pos)
			if i ==6:
				break
		#plot
		#don't really need this unless using indels...
		snp_dik = collections.OrderedDict()
		for pos in ref_dik:
			nuc, snp, ins, del_count = ref_dik.get(pos)	
			gaps = 0
			if nuc == '-':
				gaps +=1
			else:
				pos_in_seq = int(pos)-gaps
				snp_dik[pos_in_seq] = (float(snp)/float(number_hits))*100
		if snp_dik:		
			df = pd.DataFrame([snp_dik], index=['snp'])
			df_bar = df.transpose()
			#matplotlib.use('Agg')
			n = 100 #ticks
			ax = df_bar.plot.bar(figsize=(15,8), edgecolor = "none", width = 2, stacked=True)
			ticks = ax.xaxis.get_ticklocs()
			ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
			ax.xaxis.set_ticks(ticks[::n])
			ax.xaxis.set_ticklabels(ticklabels[::n])
			plt.xlabel('ref aa seq')
			plt.title(query + ' (hits in samples = ' + str(number_hits) +')')
			plt.ylabel('Number of samples with SNPs %')
			plt.savefig(query+'.svg')

def parse_blast(args, dict_key = 'assembly', dict_value = 'percent'):
	
	'''
	Parses blast output. How to do this cleaner? a Class? At least its in one place
	'''
	query_seqs = get_query_seqs(args)
	query_seqs = query_seqs.keys()

	if dict_key == 'query' and dict_value == 'percent':
		hits_dict = collections.defaultdict(list)
		with open('blast_out.txt', 'r') as fin:
			for line in fin:
				query, _, percent = line.strip().split()[:3]
				hits_dict[query].append(float(percent))
	else:
		hits_dict = {}
		with open('blast_out.txt', 'r') as fin:
			for line in fin:
				query, hit, percent = line.strip().split()[:3]
				ass = '_'.join(hit.split('_')[:-1])
				hits_dict[ass] = [str(0.0) for query in query_seqs]
		with open('blast_out.txt', 'r') as fin:
			for line in fin:
				query, hit, percent = line.strip().split()[:3]
				ass = '_'.join(hit.split('_')[:-1])
				index = query_seqs.index(query)
				if dict_key == 'assembly' and dict_value == 'percent':
					hits_dict[ass][index] = percent
				elif dict_key == 'assembly' and dict_value == 'query':
					hits_dict[ass][index] = query
				else:
					print "Can't parse Blast"
	return hits_dict, query_seqs	

def box(args, number_of_assemblies):

	'''
	Plot variation (box) and carriage (bar) on separate axis, save as svg 
	'''
	
	#Get data	
	percent_dict, labels = parse_blast(args, dict_key = 'query', dict_value = 'percent')
	variation_box = [[float(no) for no in percent_dict.get(x)] for x in percent_dict]
	carriage_bar = [100.0*(len(percent_dict.get(x))/float(number_of_assemblies)) for x in percent_dict]
	print 'variation_box,carriage_bar', variation_box,carriage_bar
	
	#plot
	plt.figure()
	fig, ax = plt.subplots()
	plt.boxplot(variation_box)
	plt.xlabel('Antigens')
	plt.xlabel('Percentage')
	plt.title('Blast')
	ax.set_xticklabels(labels,rotation='vertical')

	ax2 = ax.twinx()
	ax2.bar(ax.get_xticks(),carriage_bar,color='#eeefff',align='center')
	ax2.set_zorder(1)
	ax.set_zorder(2)
	ax2.set_alpha(0.1)
	ax.patch.set_facecolor('None')
	plt.savefig("box_plots.svg", figsize=(22,22))
	
def six_frame_translation():

	'''
	Translate nucleotide sequence to all six frames of protein
	'''
	pass

if __name__ == "__main__":
	main()
