#!/usr/bin/env python3

import screen_assembly.run.analysis as ana
import pandas as pd
import collections
import itertools

def main ():

    '''
    '''

    pass

def itol(args, assemblies):

    '''
    Makes an itol heatmap template
    '''

    name = args.query.strip().split('/')[-1][:-3]+'_itol.txt'
    fout = open(name,'w')
    #Group by percent similarity
    percent_dict, query_seqs = ana.parse_blast(args, assemblies, dict_key = 'assembly', dict_value = 'percent')
    df = pd.DataFrame.from_dict(percent_dict, orient='index')
    df = df.reindex(df.mean().sort_values(ascending=False).index, axis=1)#sort so most hits and highest homology first
    #Header
    fout.write('DATASET_HEATMAP\n')
    fout.write('SEPARATOR COMMA\n')
    fout.write('DATASET_LABEL,' + args.query[:-3] + '\n')
    fout.write('COLOR,#ff0000\n')
    FIELD_LABELS = ','.join(list(df.columns))
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



def csv(args, assemblies, binary = True):
    
    '''
    Makes to CSVs; a binary presence absence and a redundant total hits
    '''
    
    print ('Making csv ...')
    omit_set, omit_dik = ana.rejects(args)#todo make this work with specific multi hits
    hits_dict, query_seqs = ana.parse_blast(args, assemblies, dict_key = 'assembly', dict_value = 'query')
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
                            hits_per_query = ana.helper(csv, hits_per_query, query, fout, hits_dict, ass)
                        else:#omit flagged, use other hits if there are any
                            number_of_hits = len(hits_dict.get(ass).get(query))
                            if number_of_hits == 1:
                                length, percent, coords = hits_dict.get(ass).get(query)[0]
                                if coords in omit_set:
                                    fout.write('-'.join(list(itertools.chain.from_iterable(
                                        [omit_dik[ass][query][coord] for coord in omit_dik[ass][query]]))) + ',')
                                else:
                                    hits_per_query = ana.helper(csv, hits_per_query, query, fout, hits_dict, ass)
                            else:
                                not_omited = []
                                try:
                                    for tup in hits_dict.get(ass).get(query):
                                        length, percent, coords = tup
                                        if coords not in omit_set:
                                            not_omited.append(coords)
                                    if not_omited:#use the best one that passed qc
                                        hits_per_query = ana.helper(csv, hits_per_query, query, fout, hits_dict, ass, omit_set)
                                    else:#report all reasons for failing qc
                                        fout.write('-'.join(list(itertools.chain.from_iterable(
                                            [omit_dik[ass][query][coord] for coord in omit_dik[ass][query]]))) + ',')
                                except:
                                    print ('issue with', ass, query, hit_dict)        
                    else:
                        fout.write('0,')
                fout.write('\n')
    fout.close()
    
    return hits_per_query




if __name__ == "__main__":
    main()
