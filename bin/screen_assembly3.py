#!/usr/bin/env python3

import argparse
import screen_assembly.results.results as res
import screen_assembly.run.analysis as ana
import matplotlib as mpl
import os
from glob import glob
from multiprocessing import Pool, TimeoutError
import shutil
import sys

def main ():

    '''
    Screen assemblies
    Can use conitgs (nuc), CDS or protein seqs
    If assemblies and queries are nucleotide sequence then blastn used.
    If 'assemblies' are accually proteins sequences and queries are protein sequence then blastp used.
    If assemblies are nucleotide sequence and query is protein then tblastn is used.
    '''
    #Args
    parser = argparse.ArgumentParser(
            description='Screens assemblies for given gene/genes/proteins/operons.')
    
    required = parser.add_argument_group(
            'Required',
            'Need query fastas and a database to look in')
    required.add_argument(
            "-q",
            "--query",
            type=str,
            help="Query sequences fasta. Must only have one . in name. ie sample_1.fa, NOT sample.1.fa")
    required.add_argument(
            "-i",
            "--input_folder",
            type=str,
            help="Folder with ONLY db seqs in it. Can be one file or many. Must be the same sequence type.")

    group = parser.add_argument_group(
            'Options',
            'Alignment required for SNP plots, this will increase run time')
    group.add_argument(
            "-t",
            "--threads",
            type=str,
            help="How many threads to give blast, muscle, raxml etc.",
            default = '4') 
    group.add_argument(
            "-a",
            "--aln",
            action='store_false',
            help="Make alignments with muscle. On by default.",
            default=True)
    group.add_argument(
            "-x",
            "--plots",
            action='store_true',
            help="Make gap and aa variation plots. Some versions of MAC OS are incompatable with these plots and cause crashes. Off by default.",
            default=False)
    group.add_argument(
            "-p",
            "--percent_identity",
            type=str,
            help="Percent_identity. Default 80.",
            default='80')
    group.add_argument(
            "-l",
            "--percent_length",
            type=str,
            help="Percent_length. Default 80.",
            default='80')
    group.add_argument(
            '-o',
            "--operon",
            action='store_true',
            default=False,
            help='Seq is operon or kmer, something where stops are expected')
    group.add_argument(
            "-k",
            "--keep_flagged",
            action='store_true',
            help="Include hits flagged with Ns, premature stops or contig breaks. Default off.",
            default=False)
    group.add_argument(
            "-y",
            "--ignore_warnings",
            action='store_true',
            help="Ignore warnings",
            default=False)
   
    phy = parser.add_argument_group(
            'Phylogeny',
            'Uses nuc SNP aln. 1k ultrafast bootstraps.')
    phy.add_argument(
            "-r",
            "--IQtree",
            action='store_true',
            help="Run IQtree, requires nuc alignments. Off by default.",
            default=False)

    box_wisker = parser.add_argument_group(
            'box_wisker parameters',
            'Box and wisker plot')
    box_wisker.add_argument(
            "-c",
            "--label_rotation",
            type=int,
            help="Labels on box plot. Default 90.",
            default=90)
    box_wisker.add_argument(
            "-s",
            "--style",
            action='store_false',
            help="Box plot style. Default classic",
            default=True)
    box_wisker.add_argument(
            "-g",
            "--print_count",
            action='store_false',
            help="Overlay total count on percent carriage bars. Default True",
            default=True)
    box_wisker.add_argument(
            "-n",
            "--number_samples_on_box_plot",
            type=int,
            help="number_samples_on_box_plot. Default 10",
            default = 10)
  
    args = parser.parse_args()
    if args.style:
        mpl.style.use('classic')
        
    #run

    #put stuff back if re-running
    query_seqs = list(ana.get_query_seqs(args))
    for query in query_seqs:
        if os.path.exists(query):
            for gene_file in glob(query+'/*'):
                try: shutil.move(gene_file, gene_file.split('/')[-1])
                except: pass
    for log_file in glob('logs/*'):
        try: shutil.move(log_file, log_file.split('/')[-1])
        except: pass

    #start pipe    
    assemblies = ana.cat(args)
    try: assert '' not in assemblies
    except: print ('issue with assembly name', list(sorted(assemblies)))
    print ('Starting Blast...')
    blast_type = ana.blast(args)
    print ('Snipping sequences from db...')
        #Get fastas
    with Pool(processes=int(args.threads)) as pool:
        tmp = [(args, query) for query in query_seqs]
        pool.map(ana.fasta, tmp)
    print ('Looking for hits at contig breaks...')
    with Pool(processes=int(args.threads)) as pool:
        tmp = [(args, query) for query in query_seqs]
        pool.map(ana.boot_hits_at_contig_breaks, tmp)
        #Process seqs
    print ('Processing seqs (running QC (if enabled)) and translating (if required)...')
    with Pool(processes=int(args.threads)) as pool:
        tmp = [(args, blast_type, query) for query in query_seqs]
        pool.map(ana.process_seqs, tmp)
    
    with Pool(processes=int(args.threads)) as pool:
        tmp = [(args, blast_type, query) for query in query_seqs]
        pool.map(ana.process_seqs2, tmp)     
    
    print ('Making an itol template...')
    res.itol(args, assemblies)
    print ('Making a CSV summary of results...')
    hits_per_query_dik2 = res.csv(args, assemblies)
    try: ana.box(args, assemblies, blast_type)
    except Exception as e: print(e, 'cant make box plot')
    if args.aln:
        print ('Making a Box plot')
        hits_per_query_dik = ana.sanity_check(args, blast_type)
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
                    print ('hits_per_query_dik.get(query)', hits_per_query_dik.get(query)
                            , 2, 'hits_per_query_dik2.get(query)', hits_per_query_dik2.get(query))
                    sys.exit(0)
        if args.plots:
            print ('Plotting...')
            if args.operon:
                with Pool(processes=int(args.threads)) as pool:
                    tmp = [(args, query, 'nuc') for query in query_seqs]
                    pool.map(ana.plot_vars, tmp)
            else:
                ana.variant_types(args, assemblies) 
                with Pool(processes=int(args.threads)) as pool:
                    tmp = [(args, query, 'aa') for query in query_seqs]
                    pool.map(ana.plot_vars, tmp)   
            #ana.var_pos_csv(args, ana.blast_type)
    if args.IQtree:
        print ('Running IQtree... please be patient!')
        ana.gene_tree(args)
    #rename
    with open('names_map.csv', 'r') as fin:
        name_map_dict = dict(x.strip().split(',') for x in fin.readlines())
    for any_file in glob('*'):
        if 'names_map' not in any_file:
            for i in reversed(range(0,111111)): #do big first so s1 doesn't replace s11 etc
                tig_ID = 'sample_' +str(i)
                ass_file_name = name_map_dict.get(tig_ID)
                if ass_file_name:
                    fin = open(any_file, "rt")
                    data = fin.read()
                    if tig_ID + '_' in data:
                        data = data.replace(tig_ID + '_', ass_file_name + '_')
                    else:
                        data = data.replace(tig_ID, ass_file_name)
                    fin.close()
                    fin = open(any_file, "wt")
                    fin.write(data)
                    fin.close()
    #Tidy
    for query in ana.get_query_seqs(args):
        if not os.path.exists(query):
            os.makedirs(query)
        for query_output in glob(query+'*'):
            try:
                shutil.move(query_output, query + '/' + query_output)
            except Exception as e:
                #print ('wtf', e)
                pass
    for tmp in glob('*tmp*'):
        os.remove(tmp)
    ana.versions(args)
    ana.boil_down(args)
    print ('Done!')

if __name__ == "__main__":
    main()
