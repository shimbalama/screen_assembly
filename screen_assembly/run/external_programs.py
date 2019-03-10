#!/usr/bin/env python3

import screen_assembly.run.analysis as ana
from subprocess import call, Popen

def main():
    pass

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
    seq_type_query = ana.get_seq_type(args.query)
    print ('fine -> ', seq_type_query)
    print ('testing db seqs')
    cat = args.input_folder + '/concatenated_assemblies/concatenated_assemblies.fasta'
    seq_type_db = ana.get_seq_type(cat)
    print ('fine -> ', seq_type_db)
    #db
    if seq_type_db == 'DNA':
        seq_type = 'nucl'
    else:
        seq_type = 'prot'
    try:
        print ('makeblastdb')
        call(['makeblastdb',
        '-in', cat,
        '-parse_seqids',
        '-dbtype', seq_type])#'-parse_seqids',
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
    call(cmd)

def versions(args):

    '''
    Save wrapped program version to text
    '''

    with open('versions.txt', 'w') as fout:
        fout.write(str(args)+'\n')
        call(['blastn', '-version'], stdout=fout)
        if args.aln:
            call(['muscle', '-version'], stdout=fout)
        #if args.raxml:
        #    call([args.raxml_executable, '-v'], stdout=fout)
        if args.IQtree:
            call(['iqtree', '-version'], stdout=fout)

def clustal(args, query, seq_type):
	
    if not os.path.exists(query + '_' + seq_type + '_non_redundant.aln'):	
        call(['clustalo',
            '-i', query + '_' + seq_type + '_non_redundant.fasta',
            '-o', query + '_' + seq_type + '_non_redundant.aln',
            '--threads', '1'])#seems faster than just giving clustal muli threads
    
    if args.plots:
        if seq_type == 'aa':
            gap_plot(args, query, seq_type)


if __name__ == "__main__":
        main()
