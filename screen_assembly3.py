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
import random
import shutil
sys.path.append('../common_modules')#Set this specific to you 
from lab_modules import *
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MuscleCommandline
from Bio.Alphabet import generic_dna
from Bio.codonalign.codonseq import cal_dn_ds
from Bio import codonalign
import multiprocessing
from multiprocessing import Pool, TimeoutError
import matplotlib as mpl 
import itertools

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
    parser = argparse.ArgumentParser(description='Screens assemblies for given gene/genes/proteins/operons.')
    
    required = parser.add_argument_group('Required', 'Need query fastas and a database to look in')
    required.add_argument("-q", "--query", type=str, help="Query sequences fasta. Must only have one . in name. ie sample_1.fa, NOT sample.1.fa")
    required.add_argument("-i", "--input_folder", type=str, help="Folder with ONLY db seqs in it. Can be one file or many. Must be the same sequence type.")

    group = parser.add_argument_group('Options', 'Alignment required for SNP plots, this will increase run time')
    group.add_argument("-t", "--threads", type=str, help="How many threads to give blast, muscle, raxml etc.", default = '4') 
    group.add_argument("-d", "--dNdS", action='store_true', help="Calculate dNdS. Off by default.", default=False)
    group.add_argument("-a", "--aln", action='store_false', help="Make alignments with muscle. On by default.", default=True)
    group.add_argument("-x", "--plots", action='store_false', help="Make plots. On by default.", default=True)
    group.add_argument("-p", "--percent_identity", type=str, help="Percent_identity. Default 80.", default='80')
    group.add_argument("-l", "--percent_length", type=str, help="Percent_length. Default 80.", default='80')
    group.add_argument('-o', "--operon", action='store_true', default=False, help='Seq is operon or kmer, something where stops are expected')
    group.add_argument("-k", "--keep_flagged", action='store_true', help="Include hits flagged with Ns, premature stops or contig breaks. Default off.", default=False)
   
    phy = parser.add_argument_group('Phylogeny', 'Can take a long time')
    phy.add_argument("-m", "--raxml_model", type=str, help="Model for raxml. Default is GTRGAMMA", default='GTRGAMMA')
    phy.add_argument("-b", "--bootstraps", type=int, help="Number of bootstraps for raxml. Default 100.", default=100)
    phy.add_argument("-r", "--raxml", action='store_true', help="Run raxml, requires alignments. Off by default.", default=False)
    phy.add_argument("-z", "--raxml_executable", type=str, help="Default is raxmlHPC-PTHREADS-SSE3", default='raxmlHPC-PTHREADS-SSE3')

    box_wisker = parser.add_argument_group('Required', 'Box and wisker plot')
    box_wisker.add_argument("-c", "--label_rotation", type=int, help="Labels on box plot. Default 90.", default=90)
    box_wisker.add_argument("-s", "--style", action='store_false', help="Box plot style. Default classic", default=True)
    box_wisker.add_argument("-g", "--print_count", action='store_false', help="Overlay total count on percent carriage bars. Default True", default=True)
    box_wisker.add_argument("-n", "--number_samples_on_box_plot", type=int, help="number_samples_on_box_plot. Default 10", default = 10)
  
    args = parser.parse_args()
    print ('args',args)	
    if args.style:
        mpl.style.use('classic')
        
        
    #run

    #put stuff back if re-running
    query_seqs = list(get_query_seqs(args))
    for query in query_seqs:
        if os.path.exists(query):
            for gene_file in glob(query+'/*'):
                try: shutil.move(gene_file, gene_file.split('/')[-1])
                except: pass
    for log_file in glob('logs/*'):
        try: shutil.move(log_file, log_file.split('/')[-1])
        except: pass

    #start pipe    
    assemblies = cat(args)
    try: assert '' not in assemblies
    except: print ('issue with assembly name', list(sorted(assemblies)))
    print ('Starting Blast...')
    blast_type = blast(args)
    print ('Snipping sequences from db and checking for contigs breaks...')
        #Get fastas
    with Pool(processes=int(args.threads)) as pool:
        tmp = [(args, query) for query in query_seqs]
        pool.map(fasta, tmp)
        
    with Pool(processes=int(args.threads)) as pool:
        tmp = [(args, query) for query in query_seqs]
        pool.map(boot_hits_at_contig_breaks, tmp)
        #Process seqs
    with Pool(processes=int(args.threads)) as pool:
        tmp = [(args, blast_type, query) for query in query_seqs]
        pool.map(process_seqs, tmp)
        
    print ('Making an itol template...')
    itol(args, assemblies)
    print ('Making a CSV summary of results...')
    hits_per_query_dik2 = csv(args, assemblies)
    box(args, assemblies, blast_type)
    if args.aln:
        print ('Making a Box plot')
        hits_per_query_dik = sanity_check(args, blast_type)
        try:
            assert hits_per_query_dik == hits_per_query_dik2 #check the hits in the csv are teh same as in the plot
        except:
            print ('Hits per csv dont match hits in the box plot!')
            print ('from box plot (which is from query_aa_seqs.aln)',hits_per_query_dik)
            print ('from csv binary hits (from blast)', hits_per_query_dik2)
            sys.exit(0)
        if args.plots:
            print ('Plotting...')
            if args.operon:
                if blast_type != 'tblastn':#no ref in nuc
                    plot_vars(args, 'nuc')
            else:
                variant_types(args, assemblies)
                plot_vars(args)   
            var_pos_csv(args, blast_type)
    if args.raxml:
        print ('Running RAXML... please be patient!')
        gene_tree(args, blast_type)
    if not args.operon:
        if args.dNdS:
            print ('Calculating dNdS...')
            DNDS(args)
    #Tidy
    for query in get_query_seqs(args):
        if not os.path.exists(query):
            os.makedirs(query)
        for query_output in glob(query+'*'):
            try:
                shutil.move(query_output, query + '/' + query_output)
            except:
                pass#lazy dir into itself handle
    for log in glob('*RAxML_*'):
        shutil.move(log, 'logs/'+log)
    versions(args)
    boil_down(args)
    print ('Done!')

def boil_down(args):

    '''
    just report hits
    '''
    with open('binary_hits_boiled_down.csv','w') as fout:
        with open('binary_hits.csv','r') as fin:
            genes = fin.readline().strip().split(',')
            for line in fin:
                bits = line.strip().split(',')
                for i, bit in enumerate(bits):
                    if i ==0:
                        fout.write(bit+',')
                    if bit == '1':
                        fout.write(genes[i]+',')
                fout.write('\n')

def versions(args):
    
    '''
    Save wrapped program version to text
    '''

    with open('versions.txt', 'w') as fout:
        fout.write(str(args)+'\n')
        call(['blastn', '-version'], stdout=fout)
        if args.aln:
            call(['muscle', '-version'], stdout=fout)
        if args.raxml:
            call([args.raxml_executable, '-v'], stdout=fout)

def helper(csv, hits_per_query, query, fout, hits_dict, ass, omited = []):

    if csv == 'binary_hits':
        hits_per_query[query] += 1
        fout.write('1,')
    if csv == 'total_hits':
        fout.write(str(len(hits_dict.get(ass).get(query)))+',')
    if csv == 'length_and_sequence_identity':
        max_tup = 0
        biggest_tup = ()
        for tup in hits_dict.get(ass).get(query):
            length, percent, coords = tup    
            if coords not in omited:
                if sum([length, percent]) > max_tup:
                    max_tup = sum([length, percent])
                    biggest_tup = tup
        length, percent, coords = biggest_tup
        fout.write('percent_length='+str(length)+' percent_identity='+str(percent)+',')

    return hits_per_query

def csv(args, assemblies, binary = True):
    
    '''
    Makes to CSVs; a binary presence absence and a redundant total hits
    '''
    
    print ('Making csv ...')
    omit_set, omit_dik = rejects(args)#todo make this work with specific multi hits
    print ('gggggggg',omit_set, omit_dik)
    hits_dict, query_seqs = parse_blast(args, assemblies, dict_key = 'assembly', dict_value = 'query')
    for csv in ['binary_hits', 'total_hits', 'length_and_sequence_identity']:
        with open(csv + '.csv', 'w') as fout:
            header = '#Names,' + ','.join(query_seqs) +'\n'
            fout.write(header)
            if csv == 'binary_hits':
                hits_per_query = {}
                for query in query_seqs:
                    hits_per_query[query] = 0
            for ass in assemblies:
                fout.write(ass +',')
                for query in query_seqs:
                    if query in hits_dict.get(ass, []):
                        if args.keep_flagged:#don't omit flagged
                            print ('keeping all hits')
                            hits_per_query = helper(csv, hits_per_query, query, fout, hits_dict, ass)
                        else:#omit flagged, use other hits if there are any
                            number_of_hits = len(hits_dict.get(ass).get(query))
                            print (ass, query, number_of_hits)
                            if number_of_hits == 1:
                                length, percent, coords = hits_dict.get(ass).get(query)[0]
                                if coords in omit_set:
                                    print('omit', coords)
                                    fout.write('-'.join(list(itertools.chain.from_iterable([omit_dik[ass][query][coord] for coord in omit_dik[ass][query]]))) + ',')
                                else:
                                    print('dont ommit', coords)
                                    hits_per_query = helper(csv, hits_per_query, query, fout, hits_dict, ass)
                            else:
                                omited = []
                                not_omited = []
                                for tup in hits_dict.get(ass).get(query):
                                    length, percent, coords = tup
                                    if coords in omit_set:
                                        omited.append(coords)
                                    else:
                                        not_omited.append(coords)
                                print ('multi hit!!!', 'omited',omited,'not',not_omited)
                                if not_omited:#use the best one that passed qc
                                    hits_per_query = helper(csv, hits_per_query, query, fout, hits_dict, ass, omited)
                                else:#report all reasons for failing qc
                                    fout.write('-'.join(list(itertools.chain.from_iterable([omit_dik[ass][query][coord] for coord in omit_dik[ass][query]]))) + ',')
                                        
                    else:
                        fout.write('0,')
                fout.write('\n')
    fout.close()
    
    return hits_per_query


def boot_functionality(args, fout, fout2, direction, contig, query, query_seqs, record):

    last_pos = len(query) + contig.index(query)
    record_query = query_seqs.get(record.id)
    if last_pos in range(len(contig)-1, len(contig)+1):#lil wriggle room...
        SeqIO.write(record_query, fout2,'fasta')
    elif contig.index(query) in [0,1]:
        SeqIO.write(record_query, fout2,'fasta')
    else:
        SeqIO.write(record_query, fout,'fasta')


def boot_hits_at_contig_breaks(tup):

    '''
    Identifies and removes hits that touch a contig break
    '''
    args, query = tup
    cat = args.input_folder + '/concatenated_assemblies/concatenated_assemblies.fasta'
    seq_type_db = get_seq_type(cat)
    if seq_type_db != 'prot':
        query_seqs = {}
        if not os.path.exists(query + '_seqs_without_contig_break.fasta'): #don't redo if already done when re-running
            for record in SeqIO.parse(query + '_all_nuc_seqs.fasta','fasta'):
                if ':' in record.id:#sometimes randomly adds coordinates? wtf
                    id_, pos = record.id.strip().split(':')
                    record.id = id_
                    record.decription = pos
                query_seqs[record.id] = record
            fout = open(query + '_seqs_without_contig_break.fasta','w')#Don't touch contig break
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
    redundant_map = collections.defaultdict(list)
    for i, record in enumerate(SeqIO.parse(query + '_seqs_and_ref_' + seq_type + '.fasta','fasta')):
        seq_str = str(record.seq)
        redundant_map[seq_str].append(record.id)
    
     
    #write non redundant fasta
    fout= open(query + '_non_redundant.fasta','w')
    if i > 1:
        for i, seq in enumerate(redundant_map):
            if 'ref' in redundant_map.get(seq):
                fout.write('>0 ' + ','.join(redundant_map.get(seq))+'\n')
            else:
                fout.write('>'+ str(i + 1) + ' ' + ','.join(redundant_map.get(seq))+'\n')
            fout.write(seq+'\n')
    fout.close()

def add_ref(query_seqs, query, seq_type, blast_type):

    '''
    Adds to ref for alignment
    '''
    
    fout = open(query + '_seqs_and_ref_' + seq_type + '.fasta','w')
    record = query_seqs.get(query)
    if blast_type == 'tblastn' and seq_type == 'nuc':
        return fout #don't want aa ref with nuc hits 
    else:
        fout.write('>ref\n')
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
    if os.path.exists(query + '_' + seq_type + '_non_redundan_seqs.aln'):
        for i, record in enumerate(SeqIO.parse(query + '_' + seq_type + '_non_redundan_seqs.aln','fasta')):
            seq_str = str(record.seq)
            for j, aa in enumerate(seq_str):
                    if aa != '-':
                            d[j+1] += 1

        if i > 1:
            plt.plot(d.values())
            plt.axis([0, len(d), 0, max(d.values())+1])
            plt.ylabel('Unique sequences')
            plt.xlabel('Alignment length including gaps')
            plt.title('Position of indels in ' + query + ' with ' + str(max(d.values())) + ' sequences')
            plt.savefig(query +'_gaps.svg')
            plt.close()

def clustal(args, query, seq_type):
	
    print ('Clustal Omega call', query, seq_type)
    if not os.path.exists(query + '_' + seq_type + '_non_redundan_seqs.aln'):	
        call(['clustalo',
            '-i', query + '_non_redundant.fasta',
            '-o', query + '_' + seq_type + '_non_redundan_seqs.aln',
            '--threads', '1'])#seems faster than just giving clustal muli threads
    if not os.path.exists(query + '_' + seq_type + '_seqs.aln'):	
        call(['clustalo',
            '-i', query + '_seqs_and_ref_' + seq_type + '.fasta',
            '-o', query + '_' + seq_type + '_seqs.aln',
            '--threads', '1'])
    #broken
    #if seq_type == 'aa':
    #    print ('Start gap plot', query)
    #    gap_plot(args, query, seq_type)
    #    print ('Finishd gap plot', query)
    #print ('Finishd aln and gap plot', query)

   
def translated(args, query_seqs, query, seq_type, blast_type):
    
    '''
    Handles seq conversion if needed, excludes rubbish and multiple hits from further analysis
    '''
    fout = add_ref(query_seqs, query, seq_type, blast_type)
    foutN = open(query + '_seqs_and_ref_' + seq_type + '_Ns.fasta','w')
    if seq_type == 'aa':
        fout2 = open(query + '_seqs_and_ref_aa_multiple_stops.fasta','w')
    catch_multiple_hits = defaultdict(list)
    for record in SeqIO.parse(query + '_seqs_without_contig_break.fasta', 'fasta'):
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
                SeqIO.write(record,fout,'fasta')
            else:
                if record.seq.count('*') > 1:
                    SeqIO.write(record,fout2,'fasta')
                    #punt frame shifts and other errors to file for manual inspection
                else:
                    SeqIO.write(record,fout,'fasta')
    fout.close()
    foutN.close()
    if seq_type == 'aa':
        fout2.close()
 
def multi(args, query_seqs, query, seq_type, blast_type):

    catch_multiple_hits = defaultdict(list)

    qc = QC_fails(args, query)
    for record in SeqIO.parse(query + '_seqs_without_contig_break.fasta', 'fasta'):
        if record.id+':'+record.description not in qc:   
            ass = '_'.join(record.id.split('_')[:-1])
            catch_multiple_hits[ass].append(record)

    for ass in catch_multiple_hits:
        #use longest (this might mean a low quality one if pl and pi are very low)
        record_len = 0
        if len(catch_multiple_hits.get(ass)) > 1:
            print ('Sample ' + ass + ' has multiple hits for query ' + query)
        for var in catch_multiple_hits.get(ass):
            if len(str(var.seq)) > record_len:
                record_len = len(str(var.seq))
                record = var
    
    if seq_type =='aa':
        make_fasta_non_redundant(args, query, seq_type)
    if args.aln:
        clustal(args, query, seq_type)
    
	
def process_seqs(tup):

    '''
    Aligns fastas with muscle
    '''

    args, blast_type, query = tup
    print ('start aln', query)
    query_seqs = get_query_seqs(args)
    if blast_type == 'blastp':#no nuc output if prot query used
        fout = add_ref(query_seqs, query, 'aa', blast_type)
        with open(query + '_seqs_without_contig_break.fasta','r') as fin:
            for line in fin:
                fout.write(line)
        fout.close()
        make_fasta_non_redundant(args, query, 'aa')
        clustal(args, query, 'aa')
    #elif blast_type == 'tblastn':
    #    seq_type = 'aa'
    #    translated(args, query_seqs, query, seq_type, blast_type)	
    else:
        for seq_type in ['aa','nuc']:#aa first so have info re frame shifts
            if args.operon and seq_type == 'aa' and blast_type == 'blastp':
                continue
            translated(args, query_seqs, query, seq_type, blast_type)
        for seq_type in ['aa','nuc']:
            multi(args, query_seqs, query, seq_type, blast_type)
            

def gene_tree(args, directory = 'muscle_and_raxml'):

    query_seqs = get_query_seqs(args)
    for query in query_seqs:
        for i, record in enumerate(SeqIO.parse(query + '_nuc_seqs.aln', 'fasta')):
            pass
        if i > 1:#make sure it not just ref seq
            tree(args, query,  '_nuc_seqs.aln', args.raxml_executable, args.bootstraps, args.raxml_model)


def variant_types(args, assemblies):

	'''
	makes csv and itol of different seq variants
	'''
	print ('writing variant types csv')
	query_seqs = get_query_seqs(args)
	query_seqs = list(query_seqs)#order
	variants_dik = collections.defaultdict(dict)
	for query in query_seqs:
		for i, record in enumerate(SeqIO.parse(query + '_non_redundant.fasta', 'fasta')):
			if record.id == '0': #skip record.id, or ref twice
				for sample in record.description.split(' ')[1].split(','):
					if sample =='ref':
						variants_dik['Ref'][query] = '0'
					else:
						variants_dik[ass][query] = '0'
			else:
				for sample in record.description.split(' ')[1].split(','):
					ass = '_'.join(sample.split('_')[:-1])
					variants_dik[ass][query] = str(i + 1)
				
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
    
    name = args.query[:-3]+'_itol.txt'.replace('../','')
    fout = open(name,'w')
    #Group by percent similarity
    percent_dict, query_seqs = parse_blast(args, assemblies, dict_key = 'assembly', dict_value = 'percent')
    df = pd.DataFrame.from_dict(percent_dict, orient='index')
    df = df.reindex(df.mean().sort_values(ascending=False).index, axis=1)#sort so most hits and highest homology first
    #Header			
    fout.write('DATASET_HEATMAP\n')
    fout.write('SEPARATOR COMMA\n')
    fout.write('DATASET_LABEL,' + args.query[:-3] + '\n')
    fout.write('COLOR,#ff0000\n')
    FIELD_LABELS = ','.join(query_seqs)
    fout.write('FIELD_LABELS,' + FIELD_LABELS + '\n')
    '''
    fout.write('LEGEND_TITLE,Dataset_legend')
    FIELD_LABELS = ','.join(['1' for i in range(len(query_seqs))])
    fout.write('LEGEND_SHAPES,' + FIELD_LABELS + '\n')
    FIELD_LABELS = ','.join(['#ff0000' for i in range(len(query_seqs))])
    fout.write('LEGEND_COLORS,' + FIELD_LABELS + '\n')
    LEGEND_LABELS = ','.join(query_seqs)
    fout.write('LEGEND_LABELS,' + LEGEND_LABELS + '\n') 
    '''
    fout.write('DATA\n')
    fout.close()
    
    with open(name, 'a') as fout:
        df.to_csv(fout, header=False)

def var_pos_csv(args, blast_type):

    '''
    Csv of var pos for aa and nuc alns
    '''
    query_seqs = get_query_seqs(args)
    for query in query_seqs:
        for seq_type in ['nuc', 'aa']:
            if blast_type == 'tblastn' and seq_type == 'nuc':
                continue#no ref
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
    #ref
    ref_dik = collections.OrderedDict()
    print ('var posssss', seq_type, query)        
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
		plot_vars_func(args, [snp_dik], ['SNPs'], length, total, number_hits, query, seq_type + '_variants')	
		plot_vars_func(args, [ins_dik, del_dik], ['Insertions', 'Deletions'], length, total, number_hits, query,
                        seq_type + '_indels')


def parse_blast(args, assemblies, dict_key = 'assembly', dict_value = 'percent'):
	
    '''
    Parses blast output. How to do this cleaner? a Class? At least its in one place
    '''
    query_seqs = get_query_seqs(args)
    query_seqs = list(query_seqs.keys())#list for fixed order
 
    percent_length = float(args.percent_length)
    percent_identity = float(args.percent_identity) # Blast regularly calls 99% where its def 100% by direct str
    
    qc, _ = rejects(args)
    #query as key
    if dict_key == 'query' and dict_value == 'percent':#for box n wisker
        tmp = collections.defaultdict(lambda: collections.defaultdict(list)) # need tmp to get max hit v db seq
        hits_dict = collections.defaultdict(list)
        with open('blast_out.txt', 'r') as fin:
            for line in fin:
                bits = line.strip().split()
                query, hit, percent = bits[:3]
                if not args.keep_flagged:
                    if hit+':'+bits[8]+'-'+bits[9] in qc:
                        continue
                if float(percent) >= percent_identity:
                    ass = '_'.join(hit.split('_')[:-1])
                    tmp[query][ass].append(float(percent))
        for query in tmp:
            for ass in tmp.get(query):
                biggest_hit = max(tmp[query][ass])
                assert biggest_hit >= percent_identity
                hits_dict[query].append(biggest_hit)
    #ass as key
    elif dict_key == 'assembly' and dict_value == 'query': #for csv
        hits_dict = collections.defaultdict(lambda: collections.defaultdict(list))
        with open('blast_out.txt', 'r') as fin:
            for line in fin: 
                bits=line.strip().split()
                query, hit, percent = bits[:3]
                length = float(line.strip().split()[-1])
                percent = float(percent)
                if percent >= percent_identity and length >= percent_length:
                    ass = '_'.join(hit.split('_')[:-1])
                    tup = length, percent, hit+':'+bits[8]+'-'+bits[9]
                    hits_dict[ass][query].append(tup)
    #heatmap itol
    elif dict_key == 'assembly' and dict_value == 'percent':
        hits_dict = collections.defaultdict(lambda: collections.defaultdict(float))
        for ass in assemblies:#prefil
            for query in query_seqs:
                hits_dict[ass][query] = 0.0
        with open('blast_out.txt', 'r') as fin:
            for line in fin:
                bits = line.strip().split()
                query, hit, percent = bits[:3]
                percent = float(percent)
                if not args.keep_flagged:
                    if hit+':'+bits[8]+'-'+bits[9] in qc:
                        continue
                if not 'pdb' in hit and not '|' in hit:
                    if percent >= percent_identity:
                        ass = '_'.join(hit.split('_')[:-1])
                        try: assert ass in assemblies
                        except: print (ass, 'not in assemblies')
                        if percent > hits_dict.get(ass).get(query):#if 2 + use largest
                            hits_dict[ass][query] = percent#itol heatmap
                else:
                    print('skipping hit', hit, 'blast has removed contig information')
    else:
        print ("Can't parse Blast")
            
    return hits_dict, query_seqs	

def sanity_check(args, blast_type):
    
    '''
    Make sure everything lines up
    '''

    query_seqs = get_query_seqs(args)
    query_seqs = list(query_seqs.keys())

    labels = {}
    
    for query in query_seqs:
        if args.operon:
            for i, record in enumerate(SeqIO.parse(query + '_nuc_seqs.aln', 'fasta')):
                  remaining_seqs = i #not i + 1 cuz ref in there
            if blast_type == 'tblastn':
                remaining_seqs += 1 #no ref
        else:
            for i, record in enumerate(SeqIO.parse(query + '_aa_seqs.aln', 'fasta')):
                remaining_seqs = i #not i + 1 cuz ref in there

        labels[query] = remaining_seqs

    return labels

def box(args, assemblies, blast_type):

    '''
    Plot variation (box) and carriage (bar) on separate axis, save as svg 
    '''
    number_of_assemblies = len(assemblies)
    #Get data	
    percent_dict, query_seqs = parse_blast(args, assemblies, dict_key = 'query', dict_value = 'percent')

    #See how many are left
    labels = {}
    df = pd.read_csv('total_hits.csv',index_col=0,dtype='str')
    if '#Names' in df.columns:
        df.columns = list(df.columns)[1:] +['na']
    for query in query_seqs:
        labels[query] = list(df[query]).count('1')
        labels[query] += list(df[query]).count('2')
        labels[query] += list(df[query]).count('3')
        labels[query] += list(df[query]).count('4')
        labels[query] += list(df[query]).count('5')#should be enough, will come up below if not
        if query not in percent_dict:
            percent_dict[query] = [0.0]

    keys = list(labels.copy().keys())#future proof incase very multi thread it
    #- https://blog.labix.org/2008/06/27/watch-out-for-listdictkeys-in-python-3
    number_samples = args.number_samples_on_box_plot
    for box_number, query_seq_names in enumerate([keys[i:i + number_samples] for i in range(0, len(keys), number_samples)]):
        variation_box = [[float(no) for no in percent_dict.get(query)] for query in query_seq_names]
        carriage_bar = []
        for query in query_seq_names:
            try: assert labels.get(query) == len(percent_dict.get(query))  #probably not needed but I like to double check
            except: print('assert fail', query, labels.get(query), len(percent_dict.get(query)))
            if percent_dict.get(query) == [0.0]:
                carriage_bar.append(0.0)
            else:
                if args.keep_flagged:
                    labels[query] = len(percent_dict.get(query)) #overwrite above
                    carriage_bar.append(100.0*(len(percent_dict.get(query))/number_of_assemblies))
                else:
                    carriage_bar.append(100.0*(labels.get(query)/number_of_assemblies))
        #plot
        plt.figure()
        fig, ax = plt.subplots()
        plt.boxplot(variation_box, patch_artist=True)
        plt.title('Blast screen of ' + str(number_of_assemblies) + ' seqs')
        ax.set_xticklabels(query_seq_names, rotation=args.label_rotation)#,rotation=45 maybe make this optional
        plt.ylabel('Percent variation (box & wisker)')

        ax2 = ax.twinx()
        ax2.bar(ax.get_xticks(),carriage_bar,color='#eeefff',align='center')
        ax2.set_zorder(1)
        ax.set_zorder(2)
        ax2.set_alpha(0.1)
        ax.patch.set_facecolor('None')
        rects = ax2.patches
        
        if args.print_count:
            seq_counts = [str(labels.get(x)) for x in query_seq_names]
            print ('seq_counts',seq_counts)
            for rect, label in zip(rects, seq_counts):
                height = rect.get_height()
                ax2.text(rect.get_x() + rect.get_width()/2, height/2, label, ha='center', va='bottom', rotation=90)
        
        ax2.set_ylim(ax.get_ylim())
        ax.set_yticks([0.0,20.0,40.0,60.,80.0,100.0])
        ax2.set_yticks([0.0,20.0,40.0,60.,80.0,100.0])
        plt.ylabel('Percent carriage (bar)')
        plt.tight_layout()
        plt.savefig("box_and_carriage_plot" + str(box_number + 1) + ".svg", figsize=(22,22))
    print ('box labels',labels)

def DNDS(args):
    
    query_seqs = get_query_seqs(args)
    query_seqs = list(query_seqs)
    d = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(list)))
    for gene in query_seqs:
        dnds_in = gene + '_seqs_and_ref_nuc.fasta'
        try:#does all V all incase ref is an outlier
            for ref in SeqIO.parse(dnds_in,'fasta'):
                for record in SeqIO.parse(dnds_in,'fasta'):
                    if ref.id == record.id:
                        continue
                    #trim to be divisable by 3
                    min_len = min([len(str(ref.seq)), len(str(record.seq))])
                    while min_len%3!=0:
                        min_len -= 1
                    seq1 = SeqRecord(Seq(str(ref.seq)[:min_len], alphabet=IUPAC.IUPACUnambiguousDNA()), id='pro1')
                    seq2 = SeqRecord(Seq(str(record.seq)[:min_len], alphabet=IUPAC.IUPACUnambiguousDNA()), id='pro2')
                    #aln prot
                    tmp_aln = 'tmp.aln'
                    tmp_fa = 'tmp.fa'
                    with open('tmp.fa','w') as fout:
                        for seqxx in [ref, record]:
                            record_prot = SeqRecord(Seq(str(seqxx.seq),generic_dna), id=seqxx.id)
                            record_prot.seq = record_prot.seq.translate(table=11)
                            SeqIO.write(record_prot, fout, 'fasta')
                    cline = MuscleCommandline(input=tmp_fa,out=tmp_aln)
                    stdout1,stderr1 = cline()
                    #add prot
                    for i, aa_record in enumerate(SeqIO.parse('tmp.aln','fasta')):
                        if aa_record.id == ref.id:
                            pro1 = SeqRecord(Seq(str(aa_record.seq), alphabet=IUPAC.protein),id='pro1')
                        else:
                            pro2 = SeqRecord(Seq(str(aa_record.seq), alphabet=IUPAC.protein),id='pro2')
                    #make aln object
                    aln = MultipleSeqAlignment([pro1, pro2])
                    codon_aln = codonalign.build(aln, [seq1, seq2])
                    #get dnds
                    try:
                        dN, dS = cal_dn_ds(codon_aln[0], codon_aln[1], method='NG86')  
                        try: dNdS = dN/dS
                        except: dNdS = np.nan
                    except: dNdS = np.nan
                    d[gene][ref.id][record.id] = dNdS
        except:
            print ('DNDS failed for gene:', gene, 'between the seqs', ref.id, 'and', record.id, 'from ' + dnds_in)
            #print (dnds_in ,os.path.exits(dnds_in))#hash seems to break this 
    with open('dNdS_median_all_samples.csv','w') as fout:
        for gene in d:
            df = pd.DataFrame.from_dict(d.get(gene), orient='index')
            df.to_csv(gene + 'dNdS_raw.csv')
            fout.write(gene+','+str(np.median(df.median()))+'\n')

def reg(name, reject_set, reject_dik, query, reason):

    if os.path.exists(name):
        for record in SeqIO.parse(name,'fasta'):
            ass = '_'.join(record.id.split('_')[:-1])
            coords = record.description.replace(record.id,'').strip()
            exact_hit = record.id+':'+coords
            reject_dik[ass][query][exact_hit].append(reason)
            reject_set.add(exact_hit)
    return reject_set, reject_dik

def rejects(args):

    '''
    Get info pertiaing to seqs removed due to hitting a contig break or have multiple stops in aa seq
    '''
    reject_set = set([])
    reject_dik = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(list)))#specific to multiple hits
    query_seqs = get_query_seqs(args)
    for query in query_seqs:
        reject_set, reject_dik = reg(query + '_seqs_and_ref_nuc_Ns.fasta', reject_set, reject_dik, query, 'Ns')
        reject_set, reject_dik = reg(query + '_seqs_and_ref_aa_Ns.fasta', reject_set, reject_dik, query, 'Ns')
        reject_set, reject_dik = reg(query + '_nuc_seqs_excluded_due_to_contig_break.fasta', reject_set, reject_dik, query, 'Contig')
        if not args.operon:
            assert os.path.exists(query + '_seqs_and_ref_aa_multiple_stops.fasta') #make sure file has been made before using it!
            reject_set, reject_dik = reg(query + '_seqs_and_ref_aa_multiple_stops.fasta', reject_set, reject_dik, query, 'Stops')
    
    return reject_set, reject_dik

def QC_fails(args, query):

    qc = set([])
    for fasta in [query + '_seqs_and_ref_nuc_Ns.fasta',
                 query + '_seqs_and_ref_aa_Ns.fasta',
                 query + '_nuc_seqs_excluded_due_to_contig_break.fasta']:
        for record in SeqIO.parse(fasta, 'fasta'):
            coords = record.description.replace(record.id,'').strip()
            exact_hit = record.id+':'+coords
            qc.add(exact_hit)
            
    if not args.operon:
        assert os.path.exists(query + '_seqs_and_ref_aa_multiple_stops.fasta') #make sure file has been made before using it!
        for record in SeqIO.parse(query + '_seqs_and_ref_aa_multiple_stops.fasta', 'fasta'):
            coords = record.description.replace(record.id,'').strip()
            exact_hit = record.id+':'+coords
            qc.add(exact_hit)
                 
    return qc

if __name__ == "__main__":
    main()
