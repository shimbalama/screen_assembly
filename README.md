# Screen assemblies

Pipeline that screens for presence of genes of interest (GOI) in bacterial assemblies. Generates multiple CSVs and plots that describe which genes are present and how variable their sequence is. Can use DNA or protein query sequences (GOIs) and DNA contigs/fastas or protein fastas as database (db) to search in. 

If you use this script in scientific publications, then please reference Davies et al., 2019., Nature Genetics., http://dx.doi.org/10.1038/s41588-019-0417-8

## Getting Started

You need one fasta file with all GOIs as the query and a folder with db contigs/fastas. Db files can only have one '.' in the name (i.e., sample_1.fa NOT sample.1.fa) 

### Prerequisites

#### Required

Python 3

Command line blast

ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

#### Optional

Clustal Omega

http://www.clustal.org/omega/

IQtree

http://www.iqtree.org/doc/Quickstart

### Installing

pip3 install --user screen_assembly

Make sure screen_assembly3.py is in you PATH

If tkinter is missing, do 'sudo apt-get install python3-tk' (on Ubuntu)

### Check for updates

pip3 install --user screen_assembly

## Running the tests

Once screen_assembly3.py is in your PATH type screen_assembly3.py -h . If you have all dependencies then the help menu will display. Otherwise read the erorr and install whichever dependency is missing.

## Running the program

Please see the WIKI

## Authors

* **Liam McIntyre** - https://github.com/shimbalama/

## License

This project is licensed under the MIT License - see the LICENSE https://github.com/shimbalama/screen_assembly/blob/master/LICENSE file for details

## Acknowledgments

* Mark Davies lab and Jake for testing

