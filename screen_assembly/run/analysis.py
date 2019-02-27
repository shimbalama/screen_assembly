#!/usr/bin/env python3

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
    if os.path.isfile(args.input_folder + '/concatenated_assemblies/concatenated_assemblies.fasta'):
        make = False
    else:
        if not os.path.exists(args.input_folder + '/concatenated_assemblies'):
            os.makedirs(args.input_folder + '/concatenated_assemblies')
        make = True
        fout = open(args.input_folder + '/concatenated_assemblies/concatenated_assemblies.fasta', 'w')
    assemblies = set([])
    if make:
        for assembly in glob(args.input_folder + '/*'):
            ass_file_name = assembly.strip().split('/')[-1].split('.')
            if ass_file_name == ['concatenated_assemblies']:
                continue
            try:
                assert len(ass_file_name) == 2
            except:
                print ('Please ensure that your database file name only has one fullstop i.e., sample_1.fa NOT sample.1.fa')
                sys.exit(0)
            ass_file_name = ass_file_name[0].replace('#','_')#Hashes break stuff
            assemblies.add(ass_file_name)
            for i, record in enumerate(SeqIO.parse(assembly, 'fasta')):
                record.id = 'gnl|MYDB|'+ass_file_name + '_' + str(i)#can't handle all the variablity any other way
                SeqIO.write(record,fout,'fasta')
    if make:
        fout.close()
    #If pre cated
    print ('Getting assembly names...')
    with open('assembly_names.txt','w') as fout:
        for record in SeqIO.parse(args.input_folder + '/concatenated_assemblies/concatenated_assemblies.fasta', 'fasta'):
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
            #seq = str(record.seq)
            #seq_len = len(seq)
            SeqIO.write(record,fout,'fasta')
            '''
            if seq_len%3==0:
                SeqIO.write(record,fout,'fasta')
            else:
                while seq_len%3!=0:
                    seq_len -= 1
                seq = seq[:seq_len]
                try:assert len(seq)%3==0
                except: print (len(seq), len(seq)%3)
                fout.write('>'+record.id + ' ' + record.description +'\n')
                for i in range(0, seq_len, 60):
                    fout.write(seq[i:i+60] +'\n')
            '''

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
    except:
        print ('Not making unique seq fasta for',query + '_seqs_and_ref_' + seq_type)

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
