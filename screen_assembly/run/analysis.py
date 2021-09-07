#!/usr/bin/env python3

import collections
from glob import glob
from Bio import SeqIO
import os
from subprocess import call, Popen
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import sys
import random
import matplotlib as mpl 
import shutil

def main ():

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
    # if os.path.isfile(args.input_folder + '/concatenated_assemblies/concatenated_assemblies.fasta'):
    #     make = False
    # else:
    if not os.path.exists(args.input_folder + '/concatenated_assemblies'):
        os.makedirs(args.input_folder + '/concatenated_assemblies')
    make = True
    fout = open(args.input_folder + '/concatenated_assemblies/concatenated_assemblies.fasta', 'w')
    assemblies = set([])
    if make:
        with open('names_map.csv', 'w') as fout_map:
            for j, assembly in enumerate(glob(args.input_folder + '/*')):
                ass_file_name = assembly.strip().split('/')[-1].split('.')
                if ass_file_name == ['concatenated_assemblies']:
                    continue
                try:
                    assert len(ass_file_name) == 2
                except:
                    print ('Please ensure that your database file name only has one fullstop i.e., sample_1.fa NOT sample.1.fa')
                    sys.exit(0)
                ass_file_name = ass_file_name[0].replace('#','_')#Hashes break stuff
                #assemblies.add(ass_file_name)
                tig_ID = f'sample_{str(j)}'
                assemblies.add(tig_ID)
                fout_map.write(','.join([tig_ID, ass_file_name]) +'\n')
                for i, record in enumerate(SeqIO.parse(assembly, 'fasta')):
                    record.description += ' ' + ass_file_name
                    record.id = f'gnl|MYDB|sample_{str(j)}_{str(i)}'#can't handle all the variablity any other way
                    SeqIO.write(record, fout, 'fasta')
    if make:
        fout.close()
    #If pre cated
    print ('Getting assembly names...')
    with open('assembly_names.txt','w') as fout:
        for record in SeqIO.parse(args.input_folder + '/concatenated_assemblies/concatenated_assemblies.fasta', 'fasta'):
            #ass = record.id.split('|')[-1].split('_contig')[0]
            #ass = record.description.split(' ')[-1]
            ass = '_'.join(record.id.split('_')[:-1]).replace('gnl|MYDB|','')#just takes contig number off the end
            if make:
                assert ass in assemblies
            else:
                assemblies.add(ass)
            fout.write(ass+'\n')
    print (len(assemblies), 'assemblies used')

    return assemblies

def fasta(tup):

    '''
    Makes a fasta file of all hits for each query
    '''
    args, query = tup
    #make oonfig file from blast output
    current = query + '_config.txt'
    fout = open(current,'w')
    coords = []
    with open('blast_out.txt', 'r') as fin:
        for line in fin:
            bits = line.strip().split()
            if float(bits[2]) < (float(args.percent_identity) -1.0):
                continue
            if bits[0] == query:
                contig = bits[1]
                seq_start = bits[8]
                seq_end = bits[9]
                coords.append(seq_start + '-' + seq_end) #only want in one direction as just an id
                if int(seq_start) > int(seq_end):#minus/plus is right, although seems backwards..
                    line = ' '.join([contig, seq_end + '-' + seq_start, 'minus'])
                    fout.write(line + '\n')
                else:
                    line = ' '.join([contig, seq_start + '-' + seq_end, 'plus'])
                    fout.write(line + '\n')
    fout.close()
    #Cut seq of interest out of contig
    cat = args.input_folder + '/concatenated_assemblies/concatenated_assemblies.fasta'
    cmd = ['blastdbcmd',
        '-db', cat,
        '-entry_batch', current,
        '-out', query + '_all_nuc_seqs_tmp.fasta']
    call(cmd)
    #add coords to fasta as unique id for multiple hits in same contig
    with open(query + '_all_nuc_seqs.fasta','w') as fout:
        for i, record in enumerate(SeqIO.parse(query + '_all_nuc_seqs_tmp.fasta','fasta')):
            record.description = coords[i]
            record.id = str(record.id).split(':')[1]
            SeqIO.write(record,fout,'fasta')

def get_query_seqs(args):

    query_seqs = {}
    for record in SeqIO.parse(args.query, 'fasta'):
        query_seqs[record.id] = record

    return query_seqs


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
    except Exception as e:
        print ('Not making unique seq fasta for',query + '_seqs_and_ref_' + seq_type, 'because there are no hits for this query')

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
        if seq_type == 'aa' and blast_type != 'tblastn':
            try:
                record.seq = record.seq.translate(table = 11)
                SeqIO.write(record,fout,'fasta')
            except Exception as e:
                print ('Cant add ref!', e, query, seq_type, blast_type)
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

def count_invariant_sites(args, query):

    d = collections.defaultdict(lambda: collections.defaultdict(int))

    print ('Counting invar sites...')
    for record in SeqIO.parse(query + '_nuc_non_redundant.aln','fasta'):
        for pos, nuc in enumerate(str(record.seq).upper()):
            d[pos][nuc] += 1
    d_count = collections.defaultdict(int)
    for pos in d:
        if len(d.get(pos)) == 1:
            d_count[list(d.get(pos).keys())[0]] += 1
        if 'N' in list(d.get(pos).keys()):
            d_count['N'] += 1
        if '-' in list(d.get(pos).keys()):
            d_count['-'] += 1
    print (d_count)

    with open(query + '_trimmed_SNP.aln','w') as fout:
        for record in SeqIO.parse(query + '_nuc_non_redundant.aln','fasta'):
            seq = ''
            for pos, nuc in enumerate(str(record.seq).upper()):
                 if len(d.get(pos)) != 1:
                     if 'N' not in list(d.get(pos).keys()) and '-' not in list(d.get(pos).keys()):
                         seq+=nuc
            fout.write('>'+record.id+'\n')
            for i in range(0, len(seq), 60):
                fout.write(seq[i:i+60] + '\n')

    call_iqtree_with_SNP_aln(args, d_count, query)

def gene_tree(args):

    query_seqs = get_query_seqs(args)
    for query in query_seqs:
        try:
            for i, record in enumerate(SeqIO.parse(query + '_nuc_non_redundant.aln', 'fasta')):
                pass
            if i > 1:#make sure it not just ref seq
                count_invariant_sites(args, query)
            else:
                print ('Cant make tree for ', query)
        except:
            print ('Cant make tree for ', query)

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
    ax = df_bar.plot.bar(edgecolor = "none", width = 1, stacked=True)
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
                        #ass = '_'.join(hit.split('_')[:-2])
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

def call_iqtree_with_SNP_aln(args, d_count, query):

    constants=','.join([str(d_count.get('A')), str(d_count.get('C')),
           str(d_count.get('G')), str(d_count.get('T'))])
    print ('constants', query, constants)
    call(['nice', 'iqtree', '-s', query + '_trimmed_SNP.aln', '-nt', str(args.threads),
          '-m', 'GTR+G', '-bb', '1000', '-czb', '-fconst', constants])


def blast(args):

    '''
    Blasts (Nuc or Prot) seq/s of interest (Query) against db of concatenated
    assemblies(contigs)/CDS/proteins
    '''
    print ('testing query seqs')
    seq_type_query = get_seq_type(args.query)
    print ('fine -> ', seq_type_query)
    print ('testing db seqs')
    cat = args.input_folder + '/concatenated_assemblies/concatenated_assemblies.fasta'
    seq_type_db = get_seq_type(cat)
    print ('fine -> ', seq_type_db)
    #db
    if seq_type_db == 'DNA':
        seq_type = 'nucl'
    else:
        seq_type = 'prot'
    try:
        command = ['makeblastdb',
        '-in', cat,
        '-parse_seqids',
        '-dbtype', seq_type]
        
        print ('makeblastdb', ' '.join(command))
        call(command)#'-parse_seqids',
        print ('makeblastdb done')
    except:
        print ('you need to install makeblastdb or fix paths to it')

    #blast
    if seq_type_query == 'DNA' and seq_type_db == 'DNA':
        blast_type = 'blastn'
        call_blast(args, cat, 'blastn')
        print ('Running blastn...')
    elif seq_type_query == 'prot' and seq_type_db == 'prot':
        blast_type = 'blastp'
        call_blast(args, cat, 'blastp')
        print ('Running blastp...')
    elif seq_type_query == 'prot' and seq_type_db == 'DNA':
        blast_type = 'tblastn'
        call_blast(args, cat, 'tblastn')
        print ('Running tblastn...')
    else:
        print ("Error with sequence type, can't run blastn, blastp or tblastx")

    return blast_type

def call_blast(args, db, blast_type = 'blastn'):

    #cant use qcov_hsp_perc with cline (HSP = High- scoring Segment Pair )
    cmd = [blast_type,
    '-query', args.query,
    '-db', db,
    '-out', 'blast_out.txt',
    '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp',
    '-max_target_seqs', '500000',
    '-qcov_hsp_perc', args.percent_length,
    '-num_threads', args.threads]

    if blast_type == 'blastn':
        cmd += ['-perc_identity', args.percent_identity]
    print(' '.join(cmd))
    call(cmd)

def versions(args):

    '''
    Save wrapped program version to text
    '''

    with open('versions.txt', 'w') as fout:
        fout.write(str(args)+'\n')
        call(['blastn', '-version'], stdout=fout)
        if args.IQtree:
            call(['iqtree', '-version'], stdout=fout)
            
    #todo - add tot his the command that run run - ie, write out args

def clustal(args, query, seq_type):

    if not os.path.exists(query + '_' + seq_type + '_non_redundant.aln'):
        call(['clustalo',
            '-i', query + '_' + seq_type + '_non_redundant.fasta',
            '-o', query + '_' + seq_type + '_non_redundant.aln',
            '--threads', '1'])#seems faster than just giving clustal muli threads

    if args.plots:
        if seq_type == 'aa':
            gap_plot(args, query, seq_type)

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


def box(args, assemblies, blast_type):

    '''
    Plot variation (box) and carriage (bar) on separate axis, save as svg
    '''
    number_of_assemblies = len(assemblies)
    #Get data
    percent_dict, query_seqs = parse_blast(args, assemblies, dict_key = 'query', dict_value = 'percent')
    print('percent_dict', percent_dict)
    #See how many are left
    labels = collections.defaultdict(int)
    df = pd.read_csv('total_hits.csv',index_col=0,dtype='str')
    if 'Names' in df.columns:
        df.columns = list(df.columns)[1:] +['na']
    
    for i in range(111):#should be enough, will come up below if not
        for query in query_seqs:
            labels[query] += list(df[query]).count(str(i+1))
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
                    print('assert fail!!iii', query, labels.get(query), len(percent_dict.get(query)), percent_dict.get(query), labels)
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
        plt.savefig("box_and_carriage_plot" + str(box_number + 1) + ".svg")

if __name__ == "__main__":
    main()
