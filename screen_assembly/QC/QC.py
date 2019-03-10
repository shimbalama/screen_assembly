#!/usr/bin/env python3

import collections
import screen_assembly.run.analysis as ana

def main ():

    '''
    QC functions
    '''

    pass

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
    query_seqs = ana.get_query_seqs(args)
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
