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
sys.path.append('/home/lmcintyre/code/github/common_modules')#Set this specific to you 
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
    group.add_argument("-x", "--plots", action='store_false', help="Make gap and aa variation plots. Some versions of MAC OS are incompatable with these plots and cause crashes. On by default.", default=True)
    group.add_argument("-p", "--percent_identity", type=str, help="Percent_identity. Default 80.", default='80')
    group.add_argument("-l", "--percent_length", type=str, help="Percent_length. Default 80.", default='80')
    group.add_argument('-o', "--operon", action='store_true', default=False, help='Seq is operon or kmer, something where stops are expected')
    group.add_argument("-k", "--keep_flagged", action='store_true', help="Include hits flagged with Ns, premature stops or contig breaks. Default off.", default=False)
    group.add_argument("-y", "--ignore_warnings", action='store_true', help="Ignore warnings", default=False)
   
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
    print ('Snipping sequences from db...')
        #Get fastas
    with Pool(processes=int(args.threads)) as pool:
        tmp = [(args, query) for query in query_seqs]
        pool.map(fasta, tmp)
    print ('Looking for hits at contig breaks...')
    with Pool(processes=int(args.threads)) as pool:
        tmp = [(args, query) for query in query_seqs]
        pool.map(boot_hits_at_contig_breaks, tmp)
        #Process seqs
    print ('Processing seqs (running QC (if enabled)) and translating (if required)...')
    with Pool(processes=int(args.threads)) as pool:
        tmp = [(args, blast_type, query) for query in query_seqs]
        pool.map(process_seqs, tmp)
    
    with Pool(processes=int(args.threads)) as pool:
        tmp = [(args, blast_type, query) for query in query_seqs]
        pool.map(process_seqs2, tmp)     
    
    print ('Making an itol template...')
    itol(args, assemblies)
    print ('Making a CSV summary of results...')
    hits_per_query_dik2 = csv(args, assemblies)
    try: box(args, assemblies, blast_type)
    except Exception as e: print(e, 'cant make box plot')
    if args.aln:
        print ('Making a Box plot')
        hits_per_query_dik = sanity_check(args, blast_type)
        for query in hits_per_query_dik:
            try:
                assert hits_per_query_dik.get(query) == hits_per_query_dik2.get(query)
            except:
                if hits_per_query_dik.get(query) and hits_per_query_dik2.get(query):
                    print ('Problem',query,
                    set(hits_per_query_dik.get(query)).difference(set(hits_per_query_dik2.get(query))),
                    set(hits_per_query_dik2.get(query)).difference(set(hits_per_query_dik.get(query))))
                elif hits_per_query_dik.get(query) == [] and hits_per_query_dik2.get(query) == None:
                    continue
                else:
                    print (1, hits_per_query_dik.get(query) == None, 2, hits_per_query_dik2.get(query) == None)  
                if not args.ignore_warnings:
                    print ('Hits per csv dont match hits in the box plot!', query)
                    sys.exit(0)
        if args.plots:
            print ('Plotting...')
            if args.operon:
                with Pool(processes=int(args.threads)) as pool:
                    tmp = [(args, query, 'nuc') for query in query_seqs]
                    pool.map(plot_vars, tmp)
            else:
                variant_types(args, assemblies) 
                with Pool(processes=int(args.threads)) as pool:
                    tmp = [(args, query, 'aa') for query in query_seqs]
                    pool.map(plot_vars, tmp)   
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
    for tmp in glob('*tmp*'):
        os.remove(tmp)
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


def pick_top_hit(query, hits_dict, ass, omited):

    max_tup = 0
    biggest_tup = ()
    if ass =='ref':
        length, percent, coords = 100.0, 100.0, 'ref'  
    else:
        try:
            for tup in hits_dict.get(ass).get(query):
                length, percent, coords = tup
                if coords not in omited:
                    if sum([length, percent]) > max_tup:
                        max_tup = sum([length, percent])
                        biggest_tup = tup
            length, percent, coords = biggest_tup
        except:
            length, percent, coords = 0.0, 0.0, 'NA'
    return length, percent, coords

def helper(csv, hits_per_query, query, fout, hits_dict, ass, omited = set([])):

    if csv == 'binary_hits':
        hits_per_query[query].add(ass)
        fout.write('1,')
    if csv == 'total_hits':
        fout.write(str(len(hits_dict.get(ass).get(query)))+',')
    if csv == 'length_and_sequence_identity':
        length, percent, coords = pick_top_hit(query, hits_dict, ass, omited) 
        fout.write('percent_length='+str(length)+' percent_identity='+str(percent)+',')

    return hits_per_query

def csv(args, assemblies, binary = True):
    
    '''
    Makes to CSVs; a binary presence absence and a redundant total hits
    '''
    
    print ('Making csv ...')
    omit_set, omit_dik = rejects(args)#todo make this work with specific multi hits
    hits_dict, query_seqs = parse_blast(args, assemblies, dict_key = 'assembly', dict_value = 'query')
    for csv in ['binary_hits', 'total_hits', 'length_and_sequence_identity']:
        with open(csv + '.csv', 'w') as fout:
            header = 'Names,' + ','.join(query_seqs) +'\n'
            fout.write(header)
            if csv == 'binary_hits':
                hits_per_query = collections.defaultdict(set)
            for ass in assemblies:
                ass = ass.replace('gnl|MYDB|','')
                if ass == 'ref':
                    continue
                fout.write(ass +',')
                for query in query_seqs:
                    if query in hits_dict.get(ass, []):
                        if args.keep_flagged:#don't omit flagged
                            hits_per_query = helper(csv, hits_per_query, query, fout, hits_dict, ass)
                        else:#omit flagged, use other hits if there are any
                            number_of_hits = len(hits_dict.get(ass).get(query))
                            if number_of_hits == 1:
                                length, percent, coords = hits_dict.get(ass).get(query)[0]
                                if coords in omit_set:
                                    fout.write('-'.join(list(itertools.chain.from_iterable([omit_dik[ass][query][coord] for coord in omit_dik[ass][query]]))) + ',')
                                else:
                                    hits_per_query = helper(csv, hits_per_query, query, fout, hits_dict, ass)
                            else:
                                not_omited = []
                                try:
                                    for tup in hits_dict.get(ass).get(query):
                                        length, percent, coords = tup
                                        if coords not in omit_set:
                                            not_omited.append(coords)
                                    if not_omited:#use the best one that passed qc
                                        hits_per_query = helper(csv, hits_per_query, query, fout, hits_dict, ass, omit_set)
                                    else:#report all reasons for failing qc
                                        fout.write('-'.join(list(itertools.chain.from_iterable([omit_dik[ass][query][coord] for coord in omit_dik[ass][query]]))) + ',')
                                except:
                                    print ('issue with', ass, query, hit_dict)        
                    else:
                        fout.write('0,')
                fout.write('\n')
    fout.close()
    
    return hits_per_query


def boot_functionality(args, fout, fout2,  contig, query, seqs, record_hit):

    last_pos = len(query) + contig.index(query)
    
    if last_pos in range(len(contig)-1, len(contig)+1):#lil wriggle room...
        SeqIO.write(record_hit, fout2,'fasta')
    elif contig.index(query) in [0,1]:
        SeqIO.write(record_hit, fout2,'fasta')
    else:
        SeqIO.write(record_hit, fout,'fasta')


def boot_hits_at_contig_breaks(tup):

    '''
    Identifies and removes hits that touch a contig break
    '''
    args, query = tup
    cat = args.input_folder + '/concatenated_assemblies/concatenated_assemblies.fasta'
    seq_type_db = get_seq_type(cat)
    if seq_type_db == 'prot' or args.keep_flagged:
        shutil.copy(query + '_all_nuc_seqs.fasta',query + '_seqs_without_contig_break.fasta')
    else:
        seqs = collections.defaultdict(list)
        if not os.path.exists(query + '_seqs_without_contig_break.fasta'): #don't redo if already done when re-running
            for record in SeqIO.parse(query + '_all_nuc_seqs.fasta','fasta'):
                seqs[record.id].append(record)
            with open(query + '_seqs_without_contig_break.fasta','w') as fout:#Don't touch contig break
                with open(query + '_nuc_seqs_excluded_due_to_contig_break.fasta','w') as fout2:
                    for record in SeqIO.parse(cat,'fasta'):
                        seq_id = str(record.id).replace('gnl|MYDB|','')
                        if seq_id in seqs:
                            for hit in seqs.get(seq_id):
                                coords = hit.description.replace(record.id,'').strip()
                                exact_hit = hit.id+':'+coords
                                found = False
                                tmp = str(hit.seq).upper()
                                contig = str(record.seq).upper()#contig
                                if tmp in contig:
                                    boot_functionality(args, fout, fout2,  contig, tmp, seqs, hit)
                                    found = True
                                else:
                                    tmp = str(hit.seq.reverse_complement()).upper()
                                    boot_functionality(args, fout, fout2, contig, tmp, seqs, hit)
                                    found = True
                                assert found


def make_fasta_non_redundant(args, query, seq_type):

    '''
    Boil each multi fasta down to unique seqs for every query - at aa level
    '''
    try:
        redundant_map = collections.defaultdict(list)
        for i, record in enumerate(SeqIO.parse(query + '_seqs_and_ref_' + seq_type + '.fasta','fasta')):
            seq_str = str(record.seq)
            redundant_map[seq_str].append(record.id)
    
     
        #write non redundant fasta
        fout= open(query + '_' +seq_type+'_unique_seqs.fasta','w')
        if i > 1:
            for i, seq in enumerate(redundant_map):
                if 'ref' in redundant_map.get(seq):
                    fout.write('>0 ' + ','.join(redundant_map.get(seq))+'\n')
                else:
                    fout.write('>'+ str(i + 1) + ' ' + ','.join(redundant_map.get(seq))+'\n')
                fout.write(seq+'\n')
        fout.close()
    except:
        print ('Not making unique seq fasta for',query + '_seqs_and_ref_' + seq_type)

def gap_plot(args, query, seq_type):

    '''
    Line plot of aln to see where gaps are
    '''
    d= collections.defaultdict(int)
    if os.path.exists(query + '_' + seq_type + '_non_redundant.aln'):
        if not os.path.exists(query +'_gaps.svg'):
            for i, record in enumerate(SeqIO.parse(query + '_' + seq_type + '_non_redundant.aln','fasta')):
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
                plt.close('all')

def clustal(args, query, seq_type):
	
    if not os.path.exists(query + '_' + seq_type + '_non_redundant.aln'):	
        call(['clustalo',
            '-i', query + '_' + seq_type + '_non_redundant.fasta',
            '-o', query + '_' + seq_type + '_non_redundant.aln',
            '--threads', '1'])#seems faster than just giving clustal muli threads
    
    if args.plots:
        if seq_type == 'aa':
            gap_plot(args, query, seq_type)

def add_ref(args, query, seq_type, blast_type, fout):

    '''
    Adds to ref for alignment
    '''
    query_seqs = get_query_seqs(args)
    record = query_seqs.get(query)
    if blast_type == 'tblastn' and seq_type == 'nuc':
        return fout #don't want aa ref with nuc hits 
    else:
        record.id = 'ref'
        if seq_type == 'aa':
            try:
                record.seq = record.seq.translate(table = 11)
                SeqIO.write(record,fout,'fasta')
            except:
                pass #hack to not crash if is aa
        else:
            SeqIO.write(record,fout,'fasta')
    
        return fout

def flag_stops(args, fout, fout2, record):

    if args.keep_flagged:#keep seqs with stop if operon or kmer
        SeqIO.write(record,fout,'fasta')
    else:
        if str(record.seq).count('*') > 1:
            SeqIO.write(record,fout2,'fasta')
            #punt frame shifts and other errors to file for manual inspection
        else:
            SeqIO.write(record,fout,'fasta')

    return fout, fout2

def translated(args, query_seqs, query, seq_type, blast_type):
    
    '''
    Handles seq conversion if needed, excludes rubbish from further analysis
    '''
    with open(query + '_seqs_and_ref_' + seq_type + '.fasta','w') as fout:
        fout = add_ref(args, query, seq_type, blast_type, fout)
        foutN = open(query + '_seqs_and_ref_' + seq_type + '_Ns.fasta','w')
        if seq_type == 'aa':
            fout2 = open(query + '_seqs_and_ref_aa_multiple_stops.fasta','w')
        for record in SeqIO.parse(query + '_seqs_without_contig_break.fasta', 'fasta'):
            seq = str(record.seq).upper()
            if seq_type == 'nuc':
                if 'N' in seq and not args.keep_flagged:
                    SeqIO.write(record,foutN,'fasta')
                else:
                    SeqIO.write(record,fout,'fasta')
            if seq_type == 'aa':
                record.seq = record.seq.translate(table = 11)
                fout, fout2 = flag_stops(args, fout, fout2, record)
        foutN.close()
        if seq_type == 'aa':
            fout2.close()
 
def multi(args, query_seqs, query, seq_type, blast_type):

    catch_multiple_hits = collections.defaultdict(lambda:collections.defaultdict(str))
    hits_dict, query_seqs = parse_blast(args, dict_key = 'assembly', dict_value = 'query')
    qc = QC_fails(args, query)
    for record in SeqIO.parse(query + '_seqs_and_ref_' + seq_type + '.fasta', 'fasta'):
        coords = record.description.replace(record.id,'').strip()
        exact_hit = record.id+':'+coords
        if exact_hit not in qc:   
            if record.id == 'ref':
                catch_multiple_hits['ref']['na'] = record
            else:
                ass = '_'.join(record.id.split('_')[:-1])
                catch_multiple_hits[ass][exact_hit] = record
    with open(query + '_' + seq_type + '_non_redundant.fasta', 'w') as fout:
        for ass in catch_multiple_hits:
            if ass == 'ref':
                record = catch_multiple_hits['ref']['na']
            else:
                length, percent, coords = pick_top_hit(query, hits_dict, ass, qc)
                record = catch_multiple_hits[ass][coords]
            SeqIO.write(record, fout,'fasta')
    make_fasta_non_redundant(args, query, seq_type)
    if args.aln:
        clustal(args, query, seq_type)
    
	
def process_seqs(tup):

    '''
    Aligns fastas with ClustalO. Translates seqs. Removes redundancy.
    '''

    args, blast_type, query = tup
    query_seqs = get_query_seqs(args)
    if blast_type == 'blastp':#no nuc output if prot query used
        with open(query + '_seqs_and_ref_aa.fasta','w') as fout:
            fout = add_ref(args, query, 'aa', blast_type, fout)
            fout2 = open(query + '_seqs_and_ref_aa_multiple_stops.fasta','w')
            for record in SeqIO.parse(query + '_seqs_without_contig_break.fasta','fasta'):
                fout, fout2 = flag_stops(args, fout, fout2, record)
            fout2.close()
        #multi(args, query_seqs, query, 'aa', blast_type) 
    else:
        for seq_type in ['aa','nuc']:#aa first so have info re frame shifts
            if args.operon and seq_type =='aa':
                continue
            else:
                translated(args, query_seqs, query, seq_type, blast_type)
        
def process_seqs2(tup):

    args, blast_type, query = tup
    query_seqs = get_query_seqs(args)
    if blast_type == 'blastp':#no nuc output if prot query used
        multi(args, query_seqs, query, 'aa', blast_type)
    else:    
        for seq_type in ['aa','nuc']:
            if args.operon and seq_type =='aa':
                continue
            else:
                multi(args, query_seqs, query, seq_type, blast_type)
            

def gene_tree(args, directory = 'muscle_and_raxml'):

    query_seqs = get_query_seqs(args)
    for query in query_seqs:
        for i, record in enumerate(SeqIO.parse(query + '_nuc_non_redundant.aln', 'fasta')):
            pass
        if i > 1:#make sure it not just ref seq
            tree(args, query,  '_nuc_non_redundant.aln', args.raxml_executable, args.bootstraps, args.raxml_model)


def variant_types(args, assemblies):

    '''
    makes csv and itol of different seq variants
    '''
    print ('writing variant types csv')
    query_seqs = get_query_seqs(args)
    query_seqs = list(query_seqs)#order
    variants_dik = collections.defaultdict(dict)
    for query in query_seqs:
        for i, record in enumerate(SeqIO.parse(query + '_aa_unique_seqs.fasta', 'fasta')):
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
    #df = df.reindex(df.mean().sort_values(ascending=False).index, axis=1)#sort so most hits and highest homology first
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
            if not os.path.exists(query + '_' + seq_type + '_non_redundant.aln'):
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
    for record in SeqIO.parse(query + '_' + seq_type + '_non_redundant.aln', 'fasta'):
        if record.id == 'ref':
            ref_seq = str(record.seq)
            for i, nuc in enumerate(ref_seq):
                i+=1
                ref_dik[str(i)] = [nuc,0,0,0,collections.defaultdict(int)]#[ref_nuc, count_of_snps, count_of_ins, count_of_del]
    #hits
    number_hits = 0
    length = 10 # hack so doesn't fall over if no hits
    for record in SeqIO.parse(query + '_' + seq_type + '_non_redundant.aln', 'fasta'):
        if record.id == 'ref':
            continue
        number_hits += 1 
        seq = str(record.seq)
        try:
            assert len(seq) == len(ref_seq)
            length = len(seq)
        except:
            print ('alignments need to be same length!')
            try:print (len(seq),  len(ref_seq))
            except:print ('No ref_seq')
            continue
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



def plot_vars(tup):

    '''
    Plots the distribution of vars (SNPs/indels) on the seq and saves as SVG
    '''
    args, query, seq_type = tup
    print ('Plottings vars....')

    if os.path.exists(query + '_' + seq_type + '_non_redundant.aln'):
        number_hits, ref_dik, length = var_pos(args, seq_type, query)	
        if number_hits != 0 and ref_dik:
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
        else:
            print ('Skipping ' + query + '_non_redundant.aln')

def parse_blast(args, assemblies = 'na', dict_key = 'assembly', dict_value = 'percent'):
	
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
                hit = hit.replace('gnl|MYDB|','')
                if not 'pdb' in hit:
                    if not args.keep_flagged:
                        if hit+':'+bits[8]+'-'+bits[9] in qc:
                            continue
                    if float(percent) >= percent_identity:
                        ass = '_'.join(hit.split('_')[:-1])
                        tmp[query][ass].append(float(percent))
                else:
                    print('csv: skipping hit', hit, 'blast has removed contig information')
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
                hit = hit.replace('gnl|MYDB|','')
                length = float(line.strip().split()[-1])
                percent = float(percent)
                if not 'pdb' in hit:
                    if percent >= percent_identity and length >= percent_length:
                        ass = '_'.join(hit.split('_')[:-1])
                        tup = length, percent, hit+':'+bits[8]+'-'+bits[9]
                        hits_dict[ass][query].append(tup)
                else:
                    print('csv: skipping hit', hit, 'blast has removed contig information')
    #heatmap itol
    elif dict_key == 'assembly' and dict_value == 'percent':
        hits_dict = collections.defaultdict(lambda: collections.defaultdict(float))
        for ass in assemblies:#prefil
            ass = ass.replace('gnl|MYDB|','')
            for query in query_seqs:
                hits_dict[ass][query] = 0.0
        with open('blast_out.txt', 'r') as fin:
            for line in fin:
                bits = line.strip().split()
                query, hit, percent = bits[:3]
                hit = hit.replace('gnl|MYDB|','')
                percent = float(percent)
                if not args.keep_flagged:
                    if hit+':'+bits[8]+'-'+bits[9] in qc:
                        continue
                if not 'pdb' in hit:#wtf is this pdb?
                    if percent >= percent_identity:
                        ass = '_'.join(hit.split('_')[:-1])
                        try: assert ass in assemblies
                        except: print (ass, 'not in assemblies', assemblies)
                        if percent > hits_dict.get(ass).get(query):#if 2 + use largest
                            hits_dict[ass][query] = percent#itol heatmap
                else:
                    print('itol: skipping hit', hit, 'blast has removed contig information')
    else:
        print ("Can't parse Blast")
            
    return hits_dict, query_seqs	

def sanity_check(args, blast_type):
    
    '''
    Make sure everything lines up
    '''

    query_seqs = get_query_seqs(args)
    query_seqs = list(query_seqs.keys())

    labels = collections.defaultdict(set)
    
    for query in query_seqs:
        try:
            if args.operon:
                for i, record in enumerate(SeqIO.parse(query + '_nuc_non_redundant.aln', 'fasta')):
                    ass = '_'.join(record.id.split('_')[:-1])
                    if record.id != 'ref':
                        labels[query].add(ass)
            else:
                for i, record in enumerate(SeqIO.parse(query + '_aa_non_redundant.aln', 'fasta')):
                    ass = '_'.join(record.id.split('_')[:-1])
                    if record.id != 'ref':
                         labels[query].add(ass)
        except:
            labels[query] = []
            print ('No aln found for', query)

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
    if 'Names' in df.columns:
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
            try: 
                assert labels.get(query) == len(percent_dict.get(query))  #probably not needed but I like to double check
            except: 
                if labels.get(query) == 0 and percent_dict.get(query) == [0.0]:
                    pass
                else:
                    print('assert fail!!iii', query, labels.get(query), len(percent_dict.get(query)), percent_dict.get(query))
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
            for rect, label in zip(rects, seq_counts):
                height = rect.get_height()
                ax2.text(rect.get_x() + rect.get_width()/2, height/2, label, ha='center', va='bottom', rotation=90)
        
        ax2.set_ylim(ax.get_ylim())
        ax.set_yticks([0.0,20.0,40.0,60.,80.0,100.0])
        ax2.set_yticks([0.0,20.0,40.0,60.,80.0,100.0])
        plt.ylabel('Percent carriage (bar)')
        plt.tight_layout()
        plt.savefig("box_and_carriage_plot" + str(box_number + 1) + ".svg", figsize=(22,22))

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
            try: assert os.path.exists(query + '_seqs_and_ref_aa_multiple_stops.fasta') #make sure file has been made before using it!
            except: print (query, 'fail check 1: ', query + '_seqs_and_ref_aa_multiple_stops.fasta')
            reject_set, reject_dik = reg(query + '_seqs_and_ref_aa_multiple_stops.fasta', reject_set, reject_dik, query, 'Stops')
    
    return reject_set, reject_dik

def QC_fails(args, query):

    qc = set([])
    for fasta in [query + '_seqs_and_ref_nuc_Ns.fasta',
                 query + '_seqs_and_ref_aa_Ns.fasta',
                 query + '_nuc_seqs_excluded_due_to_contig_break.fasta']:
        if os.path.exists(fasta):
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
