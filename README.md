# screen_assembly - last updated on 23rd March 2018

# Applologies for the formatting - please see the wiki for an easier to read version.

Screen a bacterial assemblies (contigs/CDS or proteins) for nucleotide or protein sequences.

Required dependencies:
Blast
Muscle
common_modules/lab_modules.py (from my Github)

Optional dependencies:
Raxml

Required arguments:
-q QUERY, --query QUERY
                        Query sequences fasta.
-i INPUT_FOLDER, --input_folder INPUT_FOLDER
                        Folder with ONLY db seqs in it. Can be one file or
                        many. Must be the same sequence type. Names of files and names of contigs with files will be used.
                        Do not use fullstops ('.') in file names, i.e., use sample_1.fa NOT sample.1.fa

Optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, 
                        How many threads to give blast, muscle, raxml etc.
  -r, --raxml           Run raxml, requires alignments. Off by default.
  -b BOOTSTRAPS, 
                        Number of bootstraps for raxml. Default 100.
  -m RAXML_MODEL, 
                        Model for raxml. Default is GTRGAMMA
  -z RAXML_EXECUTABLE, 
                        Default is raxmlHPC-PTHREADS-AVX2
  -d, --dNdS            Calculate dNdS. Off by default.
  -a, --aln             Make alignments with muscle. On by default.
  -x, --plots           Make plots. On by default.
  -pi PERCENT_IDENTITY, 
                        Percent_identity. Default 80.
  -pl PERCENT_LENGTH,
                        Percent_length. Default 80
  -o, --operon          Seq is operon or kmer, something where stops are
                        expected
                        
Where a parameter is set to on/off by default use that flag to toggle the parameter to the opposite on/off.

Output:

  Top level:

    blast_out.txt - raw blast output
    binary_hits.csv, total_hits.csv and hits_as_percentage.csv are the same thing except one is binary (presence absence), one is a count of hits and one reports the actual percentage of sequence identity of hits (biggest if more than one). Negative results are reported as 0, Ns (excluded due to high level of Ns), stops (excluded due to frame shift) and contigs (excluded due to hitting a contig break). The non zero negative hits are excluded from all plots but this DOES NOT mean they are not there.
    *_itol.txt - heatmap presence absence to drop into itol on tree (tree needs same samlpe names as used in contigs)
    sequence_variants.csv - the amount of different sequences for each gene screened.
    box_and_carriage_plot1.svg - high level summary of results. Shows carriage of genes and sequence variation between hits.
 
  Per gene folder:
    *fa - various fasta files of the gene snipped from contigs
    *aln - various alignments
    *svg - various plots to show the point of variation 
    
