#!/usr/bin/env python3
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
import shutil

def main ():

	'''
	Screen assemblies
	Can use conitgs (nuc), CDS or protein seqs
	Can use six frame contig translation if don't trust the CDS caller or you're looking for sub seqs
	If assemblies and queries are nucleotide sequence then blastn used.
	If 'assemblies' are accually proteins sequences and queries are protein sequence then blastp used.
	If assemblies are nucleotide sequence and query is protein then tblastn is used.
	'''
	#2do - make cated assemblies a tmp file on external hdd, too big for mac hdd
 	#Args
	parser = argparse.ArgumentParser(description='Screens assemblies. No file will be overwiten, please manually delete to start fresh. Example usage: screen_assembly.py -q prots.fasta -i db_seqs -ph microreact-project-EkqTIX5U-tree.nwk.tre')
	#parser.add_argument("-i", "--input", type=str, help="Directory with assemblies") 
	#(rax requires that the call be in the dir)
	parser.add_argument("-q", "--query", type=str, help="Query")
	parser.add_argument("-t", "--threads", type=str, help="How many threads to give blast, muscle, raxml etc", default = '4')
	parser.add_argument("-pi", "--percent_identity", type=str, help="percent_identiy", default='80')
	parser.add_argument("-pl", "--percent_length", type=str, help="percent_length", default='80')
	parser.add_argument("-i", "--input_folder", type=str, help="Folder with ONLY db seqs in it. Can be one file or many. Must be the same sequence type.")
	parser.add_argument('-f', "--fast_mode", action='store_true', default=False, help='Disables slow stuff like the raxml run')
	parser.add_argument('-o', "--operon", action='store_true', default=False, help='Seq is operon or kmer, something where stops are expected')

	args = parser.parse_args()
	print ('args',args)	
	#run
	assemblies = cat(args)
	
	blast_type = blast(args)
	fasta(args)
	if not args.fast_mode:
		align(args, blast_type)
		gene_tree(args, blast_type)
		#if blast_type == 'blastn':
			#sum_snps(args) #- neeed to make an emble file to get syn..?
		variant_types(args, assemblies)
		plot_vars(args)
		var_pos_csv(args)
		hits_per_query_dik = box(args, assemblies)
		itol(args, assemblies)
	hits_per_query_dik2 = csv(args, assemblies)
	if not args.fast_mode:
		try:
			assert hits_per_query_dik == hits_per_query_dik2 #check the hits in the csv are teh same as in the plot
		except:
			print ('Hits per csv dont match hits in the box plot!')
			print ('from box',hits_per_query_dik)
			print ('from csv', hits_per_query_dik2)
			sys.exit(0)
	csv(args, assemblies, False) #non binary hits
	#Tidy
	for query in get_query_seqs(args):
		if not os.path.exists(query):
		    	os.makedirs(query)
		for query_output in glob(query+'*'):
			try:
				shutil.move(query_output, query + '/' + query_output)
			except:
				pass#lazy dir into itself handle
	for i in glob('*_seqs.aln.reduced'):
		os.remove(i)
	#for i in glob('*_config.txt'):
		#os.remove(i)
	directory = 'muscle_and_raxml'
	for rax_file in glob('RAxML*'):
		shutil.move(rax_file,directory + '/' + rax_file)
	for rax_file in glob('*seqs.aln'):
		shutil.move(rax_file,directory + '/' + rax_file)
	for rax_file in glob('*.fasta'):
		shutil.move(rax_file,directory + '/' + rax_file)

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
	b = False
	for record in SeqIO.parse(seqs,'fasta'):
		for pos in str(record.seq).upper():
			if pos not in nucleotides:
				DNA  = False
				protein = True
				b = True
				break #speed it up..
			elif pos not in nucleotides and pos not in amino_acids:
				print ('Error, not DNA or protein: ', pos, ' in ', record)
				protein = False
				sys.exit(0)
		break
		
	if DNA:
		seq_type = 'DNA'
	elif protein:
		seq_type = 'prot'
 		
	return seq_type


def cat(args):

	'''
	Concatenates assemblies; takes n assemblies and makes one file with all contigs
	'''
	print ('cating assemblies...')
	if os.path.isfile(args.input_folder + '/concatenated_assmblies/concatenated_assmblies.fasta'):
		make = False	
	else:
		if not os.path.exists(args.input_folder + '/concatenated_assmblies'):
			os.makedirs(args.input_folder + '/concatenated_assmblies')
		make = True
		fout = open(args.input_folder + '/concatenated_assmblies/concatenated_assmblies.fasta', 'w')
	assemblies = set([])
	#test if allready cated
	assemblies_names = glob(args.input_folder + '/*')
	if make:
		for assembly in glob(args.input_folder + '/*'):
			ass_file_name = assembly.strip().split('/')[-1].split('.')
			if ass_file_name == ['concatenated_assmblies']:
				continue	
			try:
				assert len(ass_file_name) == 2
			except:
				print ('len(ass_file_name) ==', len(ass_file_name), ass_file_name)
				sys.exit(0)
			ass_file_name = ass_file_name[0]
			assemblies.add(ass_file_name)
			for i, record in enumerate(SeqIO.parse(assembly, 'fasta')):
				record.id = ass_file_name + '_' + str(i)#can't handle all the variablity any other way
				SeqIO.write(record,fout,'fasta')
	if make:
		fout.close()
	#If pre cated
	print ('Getting assembly names...')
	for record in SeqIO.parse(args.input_folder + '/concatenated_assmblies/concatenated_assmblies.fasta', 'fasta'):
		ass = '_'.join(record.id.split('_')[:-1])
		if make:
			assert ass in assemblies
		assemblies.add(ass)
	print (len(assemblies), 'assemblies used')
	return assemblies


def call_blast(args, blast_type, cat):

 
	
	#cant use qcov_hsp_perc with cline (HSP = High- scoring Segment Pair )
	cmd = [blast_type,
	'-query', args.query,
	'-db', cat,
	'-out', 'blast_out.txt',
	'-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp',
	'-max_target_seqs', '500000',
	'-qcov_hsp_perc', args.percent_length,
	'-num_threads', args.threads]

	if blast_type == 'blastn':
		cmd += ['-perc_identity', args.percent_identity]
	call(cmd)	


def blast(args):

	'''
	Blasts (Nuc or Prot) seq/s of interest (Query) against db of concatenated 
	assemblies(contigs)/CDS/proteins
	'''
	print ('testing query seqs')
	seq_type_query = get_seq_type(args.query)
	print ('fine -> ', seq_type_query)
	print ('testing db seqs')
	cat = args.input_folder + '/concatenated_assmblies/concatenated_assmblies.fasta'
	seq_type_db = get_seq_type(cat)
	print ('fine -> ', seq_type_db)
	#db
	if seq_type_db == 'DNA':
		seq_type = 'nucl'
	else:
		seq_type = 'prot'
	try:
		call(['makeblastdb',
		'-in', cat,
		'-parse_seqids',
		'-dbtype', seq_type])
	except:
		print ('you need to install makeblastdb or fix paths to it')

	#blast
	if seq_type_query == 'DNA' and seq_type_db == 'DNA':
		blast_type = 'blastn'
		call_blast(args, 'blastn', cat)
		print ('Running blastn...')
	elif seq_type_query == 'prot' and seq_type_db == 'prot':
		blast_type = 'blastp'
		call_blast(args, 'blastp', cat)
		print ('Running blastp...')
	elif seq_type_query == 'prot' and seq_type_db == 'DNA':
		blast_type = 'tblastn'
		call_blast(args, 'tblastn', cat)
		print ('Running tblastn...')
	else:
		print ("Error with sequence type, can't run blastn, blastp or tblastx")
	
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

def csv(args, assemblies, binary = True):
	
	'''
	Makes to CSVs; a binary presence absence and a redundant total hits
	'''
	
	print ('Making csv ...')
	_, omit_dik = rejects(args)

	hits_dict, query_seqs = parse_blast(args, assemblies, dict_key = 'assembly', dict_value = 'query')
	if binary:
		fout = open('binary_hits.csv','w')
	else:
		fout = open('total_hits.csv','w')
	header = '#Names,' + ','.join(query_seqs) +'\n'
	fout.write(header)
	if binary:
		hits_per_query = {}
		for query in query_seqs:
			hits_per_query[query] = 0
	for ass in assemblies:
		fout.write(ass +',')
		for query in query_seqs:
			if query in hits_dict.get(ass, []):
				if binary:
					hits_per_query[query] += 1
					fout.write('1,')
				else:
					fout.write(str(hits_dict.get(ass).get(query))+',')
			else:
				if omit_dik[ass][query]:
					fout.write(omit_dik[ass][query][0] + ',')
				else:
					fout.write('0,')
		fout.write('\n')
	fout.close()
	if binary:
		return hits_per_query

def fasta(args,  blast_folder = 'blast'):

	'''
	Makes a fasta file of all hits for each query
	'''

	print ('Cutting seqs from db for fasta ...')

	query_seqs = get_query_seqs(args)
	for query in query_seqs:
		
		#make oonfig file from blast output
		current = query + '_config.txt'
		fout = open(current,'w')
		with open('blast_out.txt', 'r') as fin:
			for line in fin:
				bits = line.strip().split()
				if float(bits[2]) < (float(args.percent_identity) -1.0):
					continue
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
		cat = args.input_folder + '/concatenated_assmblies/concatenated_assmblies.fasta'
		cmd = ['blastdbcmd',
			'-db', cat,
			'-entry_batch', current,
			'-out', query + '_seqs_with_contig_breaks.fasta']
		call(cmd)
		boot_hits_at_contig_breaks(args, query, cat)

def boot_functionality(args, fout, fout2, direction, contig, query, query_seqs, record):

	last_pos = len(query) + contig.index(query)
	if last_pos in range(len(contig)-1, len(contig)+1):#lil wriggle room...
		fout2.write('>' + record.id + ' ' + direction +' strand end of contig\n')
		fout2.write(query + '\n')
	elif contig.index(query) in [0,1]:
		fout2.write('>' + record.id + ' ' + direction +' strand start of contig\n')
		fout2.write(query + '\n')
	else:
		record_query = query_seqs.get(record.id)
		SeqIO.write(record_query, fout,'fasta')


def boot_hits_at_contig_breaks(args, query, cat):

	'''
	Identifies and removes hits that touch a contig break
	'''
	seq_type_db = get_seq_type(cat)
	if seq_type_db != 'prot':

		print ('Identiying contig breaks...')
		query_seqs = {}
		for record in SeqIO.parse(query + '_seqs_with_contig_breaks.fasta','fasta'):
			query_seqs[record.id] = record
		fout = open(query + '_seqs.fasta','w')#Don't touch contig break
		fout2 = open(query + '_nuc_seqs_excluded_due_to_contig_break.fasta','w')
		for record in SeqIO.parse(cat,'fasta'):
			if record.id in query_seqs:
				found = False
				query = str(query_seqs.get(record.id).seq).upper()
				contig = str(record.seq).upper()#contig
				if query in contig:
					boot_functionality(args, fout, fout2, 'Fwd', contig, query, query_seqs, record)
					found = True
				else:
					record.seq = record.seq.reverse_complement()
					contig = str(record.seq).upper()
					boot_functionality(args, fout, fout2, 'Rev', contig, query, query_seqs, record)
					found = True
				assert found
		fout.close()
		fout2.close()


def make_fasta_non_redundant(args, query, seq_type):

	'''
	Boil each multi fasta down to unique seqs for every query - at aa level
	'''
	print ('Eliminating redundant seqs...')	
	redundant_map = {}
	for i, record in enumerate(SeqIO.parse(query + '_seqs_and_ref_' + seq_type + '.fasta','fasta')):
		seq_str = str(record.seq)
		if seq_str not in redundant_map:
			redundant_map[seq_str] = [record.id]
		else:
			redundant_map[seq_str].append(record.id)
	print (query,'total seqs = ',i+1, ' .Non redundant = ',len(redundant_map))
	
	#write non redundant fasta
	fout= open(query + '_non_redundant.fasta','w')
	for i, seq in enumerate(redundant_map):
		if 'ref' in redundant_map.get(seq):
			fout.write('>0 ' + ','.join(redundant_map.get(seq))+'\n')
		else:
			#random_name = redundant_map.get(seq).pop()
			fout.write('>'+ str(i + 1) + ' ' + ','.join(redundant_map.get(seq))+'\n')
		fout.write(seq+'\n')
	fout.close()

def add_ref(query_seqs, query, seq_type, blast_type):

	'''
	Adds to ref for alignment
	'''
	fout = open(query + '_seqs_and_ref_' + seq_type + '.fasta','w')
	fout.write('>ref\n')
	record = query_seqs.get(query)
	print (query_seqs)
	if seq_type == 'aa':
		try:
			record.seq = record.seq.translate(table = 11)
		except:
			pass #hack to not crash if is aa
	fout.write(str(record.seq) + '\n')
	
	return fout

def gap_plot(args, query, seq_type):

	'''
	Line plot of aln to see where gaps are
	'''
	d= collections.defaultdict(int)
	for i, record in enumerate(SeqIO.parse(query + '_' + seq_type + '_non_redundan_seqs.aln','fasta')):
		seq_str = str(record.seq)
		for j, aa in enumerate(seq_str):
			if aa != '-':
				d[j+1] += 1

	plt.plot(d.values())
	plt.axis([0, len(d), 0, max(d.values())+1])
	plt.ylabel('Unique sequences')
	plt.xlabel('Alignment length including gaps')
	plt.title('Position of indels in ' + query + ' with ' + str(max(d.values())) + ' sequences')
	plt.savefig(query +'_gaps.svg')
	plt.close()

def muscle(args, query, seq_type):
	
	call(['muscle',
		'-in', query + '_non_redundant.fasta',
		'-out', query + '_' + seq_type + '_non_redundan_seqs.aln'])
	call(['muscle',
		'-in', query + '_seqs_and_ref_' + seq_type + '.fasta',
		'-out', query + '_' + seq_type + '_seqs.aln'])
	if seq_type == 'aa':
		gap_plot(args, query, seq_type)

def translated(args, query_seqs, query, seq_type, blast_type):
	
	'''
	Handles seq conversion if needed, excludes rubbish and multiple hits from further analysis
	'''
	fout = add_ref(query_seqs, query, seq_type, blast_type)
	foutN = open(query + '_seqs_and_ref_' + seq_type + '_Ns.fasta','w')
	if seq_type == 'aa':
		fout2 = open(query + '_seqs_and_ref_' + seq_type + '_multiple_stops.fasta','w')
	catch_multiple_hits = set([])
	for record in SeqIO.parse(query + '_seqs.fasta', 'fasta'):
		ass = '_'.join(record.id.split('_')[:-1])
		if ass in catch_multiple_hits:
			print ('record excluded due to multple hits:', record)
			continue
		else:
			catch_multiple_hits.add(ass)
		seq = str(record.seq).upper()
		if seq_type == 'nuc':
			if 'N' in seq:
				SeqIO.write(record,foutN,'fasta')
			else:
				SeqIO.write(record,fout,'fasta')
		if seq_type == 'aa':
			if 'N' in seq:
				SeqIO.write(record,foutN,'fasta')
				continue#got to double up for when have prot query...
			record.seq = record.seq.translate(table = 11)
			if args.operon:#keep seqs with stop if operon or kmer
				print ('not omiting seqs with stops')
				SeqIO.write(record,fout,'fasta')
			else:
				if record.seq.count('*') > 1:
					SeqIO.write(record,fout2,'fasta')
					#punt frame shifts and other errors to file for manual inspection
				else:
					SeqIO.write(record,fout,'fasta')
	fout.close()
	foutN.close()
	if seq_type =='aa':
		fout2.close()
		make_fasta_non_redundant(args, query, seq_type)
	muscle(args, query, seq_type)

	
def align(args, blast_type, directory = 'muscle_and_raxml'):

	'''
	Aligns fastas with muscle
	'''

	if not os.path.exists(directory):
		    os.makedirs(directory)

	query_seqs = get_query_seqs(args)
	for query in query_seqs:
		if blast_type == 'blastp':#no nuc output if prot query used
			fout = add_ref(query_seqs, query, 'aa', blast_type)
			with open(query + '_seqs.fasta','r') as fin:
				for line in fin:
					fout.write(line)
			fout.close()
			make_fasta_non_redundant(args, query, 'aa')
			muscle(args, query, 'aa')
		elif blast_type == 'tblastn':#this dup needs to be a function - if it works
			seq_type = 'aa'
			#add ref - should this be separate file?
			translated(args, query_seqs, query, seq_type, blast_type)	
		else:
			for seq_type in ['nuc','aa']:
				translated(args, query_seqs, query, seq_type, blast_type)


def gene_tree(args, directory = 'muscle_and_raxml'):

	'''
	Makes a Raxml tree of alignments
	'''
	model = 'PROTGAMMAWAG'
	query_seqs = get_query_seqs(args)
	for query in query_seqs:
		cmd = ['raxmlHPC-PTHREADS',
			'-s', query + '_aa_non_redundan_seqs.aln',
			'-n', query + '.tre',
			'-m', model,
			'-p', '12345',
			'-T', args.threads]
		call(cmd)

def sum_snps(args):

	from Bio.Alphabet import generic_dna, generic_protein
	from Bio import SeqFeature
	from Bio.SeqFeature import FeatureLocation
	'''
	Calls sum_snps.py
	'''
	print ('Running sum snps...')	
	#can't get teh synonomous calulation to work..	
	query_seqs = get_query_seqs(args)
	for query in query_seqs:
		from Bio import SeqFeature
		for record in SeqIO.parse(query + '_nuc_seqs.aln','fasta'):
			if record.id == 'ref':
				my_start_pos = SeqFeature.ExactPosition(1)
				my_end_pos = SeqFeature.ExactPosition(len(record.seq)-1)
				my_feature_location = FeatureLocation(my_start_pos,my_end_pos)
				my_feature_type = "CDS"
				record.seq.alphabet = generic_dna
				from Bio.SeqFeature import SeqFeature#pretty sure this is bio pythons fault...
				my_feature = SeqFeature(my_feature_location, type=my_feature_type,qualifiers={'codon_start':1,'transl_table':11})
				record.features.append(my_feature)
				with open(query + '.gbk','w') as fout:
					SeqIO.write(record,fout,'genbank')

	for query in query_seqs:
		cmd = ['summarise_snps_local.py',
			'-i', query + '_nuc_seqs.aln',
			'-r', 'ref',
			'-e', query + '.gbk',
			'-t',
			'-o', query + '_sum_snps']
		call(cmd)
	
def genome_tree(args):

	'''
	Takes a tre file and anotes it with with gene variants
	'''
	#superceded by itol
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

def variant_types(args, assemblies):

	'''
	makes csv and itol of different seq variants
	'''
	print ('writing variant types csv')
	query_seqs = get_query_seqs(args)
	query_seqs = list(query_seqs)#order
	variants_dik = collections.defaultdict(dict)
	print (query_seqs)
	for query in query_seqs:
		for i, record in enumerate(SeqIO.parse(query + '_non_redundant.fasta', 'fasta')):
			if record.id == '0': #skip record.id, or ref twice
				for sample in record.description.split(' ')[1].split(','):
					if sample =='ref':
						variants_dik['Ref'][query] = '0'
					else:
						ass = '_'.join(sample.split('_')[:-1])
						variants_dik[ass][query] = '0'
			else:
				for sample in record.description.split(' ')[1].split(','):
					ass = '_'.join(sample.split('_')[:-1])
					variants_dik[ass][query] = str(i + 1)
				
	print (variants_dik)	
	fout = open('sequence_variants.csv', 'w')
	header = 'Sample,'+','.join(query_seqs)
	fout.write(header + '\n')
	for sample in assemblies:
		tw = []
		for query in query_seqs:
			if variants_dik[sample]:
				tw.append(variants_dik.get(sample, 'NA').get(query, 'NA'))
			else:
				tw.append('NA')
		tw = sample + ',' + ','.join(tw)
		fout.write(tw+'\n')
	fout.close()



def itol(args, assemblies):

	'''
	Makes an itol heatmap template 
	'''
	color_dict = {}
	fout = open(args.query[:-3]+'_itol.txt','w')
	#Group by percent similarity
	percent_dict, query_seqs = parse_blast(args, assemblies, dict_key = 'assembly', dict_value = 'percent')
	#Header			
	fout.write('DATASET_HEATMAP\n')
	fout.write('SEPARATOR SPACE\n')
	fout.write('DATASET_LABEL ' + args.query[:-3] + '\n')
	fout.write('COLOR #ff0000\n')
	FIELD_LABELS = ' '.join(query_seqs)
	fout.write('FIELD_LABELS ' + FIELD_LABELS + '\n')
	fout.write('DATA\n')
	for ass in percent_dict:
		fout.write(ass+' ')
		for query_percent in percent_dict.get(ass):
			fout.write(str(query_percent) + ' ')
		fout.write('\n')
	fout.close()

def var_pos_csv(args):

	'''
	Csv of var pos for aa and nuc alns
	'''
	query_seqs = get_query_seqs(args)
	for query in query_seqs:
		for seq_type in ['nuc', 'aa']:
			if not os.path.exists(query + '_' + seq_type + '_seqs.aln'):
				continue # skip nuc if using aa and or any seqs without hits
			number_hits, ref_dik, length = var_pos(args, seq_type, query)
			fout = open(query + '_' + seq_type + '_pos.csv','w')#checked 
			for pos in ref_dik:
				fout.write(ref_dik.get(pos)[0] +',')
				for var in ref_dik.get(pos)[4]:
					count = str(ref_dik.get(pos)[4].get(var))
					fout.write(var + ',' + count + ',')
				fout.write('\n')
			fout.close()

def var_pos(args, seq_type, query):

	'''
	Get positions of vars in aa or nuc seqs
	'''
	print ('Get positions of vars: ', seq_type, query)	
	#ref
	ref_dik = collections.OrderedDict()
	
	for record in SeqIO.parse(query + '_' + seq_type + '_seqs.aln', 'fasta'):
		if record.id == 'ref':
			ref_seq = str(record.seq)
			for i, nuc in enumerate(ref_seq):
				i+=1
				ref_dik[str(i)] = [nuc,0,0,0,collections.defaultdict(int)]#[ref_nuc, count_of_snps, count_of_ins, count_of_del]
	#hits
	number_hits = 0
	length = 10 # hack so doesn't fall over if no hits
	for record in SeqIO.parse(query + '_' + seq_type + '_seqs.aln', 'fasta'):
		if record.id == 'ref':
			continue
		number_hits += 1 
		seq = str(record.seq)
		try:
			assert len(seq) == len(ref_seq)
			length = len(seq)
		except:
			print ('alignments need to be same length!')
			print (len(seq),  len(ref_seq))
		for i, nuc in enumerate(seq):
			i+=1
			if ref_dik.get(str(i))[0] == '-' and nuc == '-':
				continue
			elif ref_dik.get(str(i))[0] == '-':#ref is - & nuc != - == ins
				ref_dik.get(str(i))[2] += 1
			elif nuc == '-': #current seq is - == del
				ref_dik.get(str(i))[3] += 1
			else:
				if ref_dik.get(str(i))[0] != nuc:
					ref_dik.get(str(i))[1] += 1
					ref_dik.get(str(i))[4][nuc] += 1
	return number_hits, ref_dik, length

def plot_vars_func(args, list_of_diks, list_of_dik_names, length, total, number_hits, query, seq_types):

	'''
	PLot vars v ref seq
	'''

	df = pd.DataFrame(list_of_diks, index=list_of_dik_names)
	df_bar = df.transpose()
	#matplotlib.use('Agg')
	n = 100 #ticks
	ax = df_bar.plot.bar(figsize=(15,8), edgecolor = "none", width = 1, stacked=True)
	ticks = ax.xaxis.get_ticklocs()
	ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
	yticks = [i for i in range(101)]
	ax.yaxis.set_ticks(yticks[0:101:20])
	ax.xaxis.set_ticks(ticks[::n])
	ax.xaxis.set_ticklabels(ticklabels[::n])
	plt.xlabel('Reference lenght: ' + str(length) + ' positions with var: ' + str(total))
	plt.title(query + ' (samples = ' + str(number_hits) +')')
	plt.ylabel('Number of samples with ' + seq_types + ' %')
	plt.savefig(query+'_'+seq_types+'.svg')



def plot_vars(args, seq_type = 'aa'):

	'''
	Plots the distribution of vars (SNPs/indels) on the seq and saves as SVG
	'''
	print ('Plottings vars....')
	query_seqs = get_query_seqs(args)
	for query in query_seqs:
		if not os.path.exists(query + '_' + seq_type + '_seqs.aln'):
			print ('Skipping ' + query + '_seqs.aln')
			continue
		number_hits, ref_dik, length = var_pos(args, seq_type, query)	
		if number_hits == 0:
			continue
		#Get SNP/indel positions in gap free ref seq
		snp_dik = collections.OrderedDict()
		ins_dik = collections.OrderedDict()
		del_dik = collections.OrderedDict()
		total = 0
		gaps = 0 # Want to plot ref seq without '-'s the pos csv has them so will be out of sunc at end
		for pos in ref_dik:
			pos_in_seq = int(pos)-gaps
			nuc, snp, ins, del_count, _ = ref_dik.get(pos)	
			if int(snp) != 0:
				total +=1
			if nuc == '-':#not perfect - colapses ins and overwrites with last count... but how else..?
				ins_dik[pos_in_seq] = (float(ins)/float(number_hits))*100 
				gaps +=1
			else:
				#ins_dik[pos_in_seq] = 0.0
				if del_count == 0:
					del_dik[pos_in_seq] = 0.0
				else:
					del_dik[pos_in_seq] = (float(del_count)/float(number_hits))*100
				if snp == 0:
					snp_dik[pos_in_seq] = 0.0
				else:
					snp_dik[pos_in_seq] = (float(snp)/float(number_hits))*100
		#plot
		plot_vars_func(args, [snp_dik], ['SNPs'], length, total, number_hits, query, 'aa_variants')	
		plot_vars_func(args, [ins_dik, del_dik], ['Insertions', 'Deletions'], length, total, number_hits, query, 'aa_indels')


def parse_blast(args, assemblies, dict_key = 'assembly', dict_value = 'percent'):
	
	'''
	Parses blast output. How to do this cleaner? a Class? At least its in one place
	'''
	query_seqs = get_query_seqs(args)
	query_seqs = list(query_seqs.keys())#list for fixed order
	omit, _ = rejects(args)

	percent_identity = float(args.percent_identity) - 1.0 # Blast regularly calls 99% where its def 100% by direct str
	#query as key
	if dict_key == 'query' and dict_value == 'percent':
		tmp = collections.defaultdict(lambda: collections.defaultdict(list)) # need tmp to get max hit v db seq
		hits_dict = collections.defaultdict(list)
		with open('blast_out.txt', 'r') as fin:
			for line in fin:
				query, hit, percent = line.strip().split()[:3]
				if hit not in omit.get(query, 'NA'):
					if float(percent) >= percent_identity:
						ass = '_'.join(hit.split('_')[:-1])
						tmp[query][ass].append(float(percent))
		for query in tmp:
			#print 'tmp.get(query)',tmp.get(query)
			for ass in tmp.get(query):
				biggest_hit = max(tmp[query][ass])
				assert biggest_hit >= percent_identity
				hits_dict[query].append(biggest_hit)
	#ass as key
	elif dict_key == 'assembly' and dict_value == 'query': #for csv
		hits_dict = collections.defaultdict(lambda: collections.defaultdict(int))
		with open('blast_out.txt', 'r') as fin:
			for line in fin: 
				query, hit, percent = line.strip().split()[:3]
				if hit not in omit.get(query, 'NA'):
					if float(percent) >= percent_identity:
						ass = '_'.join(hit.split('_')[:-1])
						hits_dict[ass][query] += 1
	#heatmap itol
	elif dict_key == 'assembly' and dict_value == 'percent':
		hits_dict = {}
		#fill dict with 0.0
		for ass in assemblies:
			hits_dict[ass] = [str(0.0) for query in query_seqs]
		#Positionally replace 0.0 where there's hits
		with open('blast_out.txt', 'r') as fin:
			for line in fin:
				query, hit, percent = line.strip().split()[:3]
				if hit not in omit.get(query, 'NA'):
					if float(percent) >= percent_identity:
						ass = '_'.join(hit.split('_')[:-1])
						index = query_seqs.index(query)
						if float(percent) > float(hits_dict[ass][index]):#if 2 + use largest
							hits_dict[ass][index] = percent#itol heatmap
	else:
		print ("Can't parse Blast")
		
	return hits_dict, query_seqs	

def box(args, assemblies):

	'''
	Plot variation (box) and carriage (bar) on separate axis, save as svg 
	'''
	number_of_assemblies = len(assemblies)
	#Get data	
	percent_dict, query_seqs = parse_blast(args, assemblies, dict_key = 'query', dict_value = 'percent')
	print ('percent_dict.keys()', percent_dict.keys(), percent_dict)

	for i in percent_dict:
		print (i, len(percent_dict.get(i)))
	#See how many are left
	labels = {}#try use ordered dict if still no order
	for query in query_seqs:
		for i, record in enumerate(SeqIO.parse(query + '_aa_seqs.aln', 'fasta')):
			#print record.id 
			remaining_seqs = i #not i + 1 cuz ref in there
		labels[query] = remaining_seqs
		if query not in percent_dict:
			percent_dict[query] = [0.0]
		print (query, remaining_seqs)
	print ('percent_dict.keys()', percent_dict.keys(), percent_dict)
	keys = list(labels.copy().keys())#future proof incase very multi thread it
	#- https://blog.labix.org/2008/06/27/watch-out-for-listdictkeys-in-python-3
	for box_number, query_seq_names in enumerate([keys[i:i + 10] for i in range(0, len(keys), 10)]):
		print ('query_seq_names',query_seq_names)
		variation_box = [[float(no) for no in percent_dict.get(query)] for query in percent_dict if query in query_seq_names]
		carriage_bar = []
		for query in percent_dict:
			if query in query_seq_names:
				if percent_dict.get(query) == [0.0]:
					carriage_bar.append(0.0)
				else:
					carriage_bar.append(100.0*(len(percent_dict.get(query))/float(number_of_assemblies)))
		#carriage_bar = [100.0*(len(percent_dict.get(query))/float(number_of_assemblies)) for query in percent_dict if query in query_seq_names]
		print ('carriage_bar',carriage_bar)
		print ('variation_box',variation_box)
		#plot
		plt.figure()
		fig, ax = plt.subplots()
		plt.boxplot(variation_box)
		plt.title('Blast screen of ' + str(number_of_assemblies) + ' seqs')
		ax.set_xticklabels(query_seq_names,rotation=45)
		ax.set_yticks([0.0,20.0,40.0,60.,80.0,100.0])
		plt.ylabel('Percent variation (box & wisker)')

		ax2 = ax.twinx()
		#ax2.set_yticks([0.0,20.0,40.0,60.,80.0,100.0])
		ax2.bar(ax.get_xticks(),carriage_bar,color='#eeefff',align='center')
		ax2.set_zorder(1)
		ax.set_zorder(2)
		ax2.set_alpha(0.1)
		ax.patch.set_facecolor('None')
		rects = ax2.patches
		seq_counts = [str(labels.get(x)) for x in query_seq_names]
		print ('seq_counts',seq_counts)
		for rect, label in zip(rects, seq_counts):
			height = rect.get_height()
			ax2.text(rect.get_x() + rect.get_width()/2, height/2, label, ha='center', va='bottom')

		ax2.set_yticks([0.0,20.0,40.0,60.,80.0,100.0])
		plt.ylabel('Percent carriage (bar)')
		plt.tight_layout()
		plt.savefig("box_and_carriage_plot" + str(box_number + 1) + ".svg", figsize=(22,22))
	return labels

def six_frame_translation():

	'''
	Translate nucleotide sequence to all six frames of protein
	'''
	pass

def reg (name, reject_set, reject_dik, query, reason):

	if os.path.exists(name):
		for record in SeqIO.parse(name,'fasta'):
			ass = '_'.join(record.id.split('_')[:-1])
			reject_dik[ass][query].append(reason)
			reject_set[query].add(record.id)
	return reject_set, reject_dik

def rejects(args):

	'''
	Get info pertiaing to seqs removed due to hitting a contig break or have multiple stops in aa seq
	'''
	reject_set = collections.defaultdict(set)
	reject_dik = collections.defaultdict(lambda: collections.defaultdict(list))
	query_seqs = get_query_seqs(args)
	for query in query_seqs:
		reject_set, reject_dik = reg(query + '_seqs_and_ref_nuc_Ns.fasta', reject_set, reject_dik, query, 'Ns')
		reject_set, reject_dik = reg(query + '_seqs_and_ref_aa_Ns.fasta', reject_set, reject_dik, query, 'Ns')
		reject_set, reject_dik = reg(query + '_nuc_seqs_excluded_due_to_contig_break.fasta', reject_set, reject_dik, query, 'Contig')
		#reject_set, reject_dik = reg(query + '_seqs_with_contig_breaks.fasta', reject_set, reject_dik, query, 'Contig')
		if not args.operon:
			print ('omiting stops...')
			reject_set, reject_dik = reg(query + '_seqs_and_ref_aa_multiple_stops.fasta', reject_set, reject_dik, query, 'Stops')
	return reject_set, reject_dik

if __name__ == "__main__":
	main()
