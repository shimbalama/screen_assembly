#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import sys
import random
import matplotlib as mpl 
import screen_assembly.run.analysis as ana

def main ():

    '''
    Functions for plotting
    '''

    pass

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
    percent_dict, query_seqs = ana.parse_blast(args, assemblies, dict_key = 'query', dict_value = 'percent')

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

if __name__ == "__main__":
    main()


