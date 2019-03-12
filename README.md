# Screen assemblies

Pipeline that screens for presence of genes of interest (GOI) in bacterial assemblies. Generates multiple CSVs and plots that describe which genes are present and how variable their sequence is. Can use DNA or protein query sequences (GOIs) and DNA contigs/fastas or protein fastas as database (db) to search in. 

## Getting Started

You need one fasta file with all GOIs as the query and a folder with db contigs/fastas. Db files can only have one '.' in the name (i.e., sample_1.fa NOT sample.1.fa) 

### Prerequisites

#### Required

Python 3 and scypi/Biopython
Command line blast
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

#### Optional

Clustal Omega,
http://www.clustal.org/omega/
IQtree

### Installing

* Download the screen_assembly3.py script and place it in your PATH: 
  * git clone https://github.com/shimbalama/screen_assembly.git
  * Make sure its executable (chmod +x screen_assembly/screen_assembly3.py)
  * Export PATH="your_path:$PATH" (the command pwd will give you your PATH)
  * Best to permanently add it to you path by adding it to .bash_profile (mac) or .profile (unix)
* Download lab_modules.py:
  * git clone https://github.com/shimbalama/common_modules.git
  * Make sure its executable (chmod +x common_modules/lab_modules.py)
  * Export PYTHONPATH="your_path:$PYTHONPATH"
  * Best to permanently add it to you path by adding it to .bash_profile (mac) or .profile (unix)
* Place the common_modules folder next to the screen_assembly folder (as thats where it looks by default). OR use a text editor to set this line in screen_assembly3.py to point at the dir you put lab_modules.py in:  sys.path.append('../common_modules') becomes sys.path.append('your_path/common_modules')

### Check for updates

* git pull

## Running the tests

Once screen_assembly3.py is in your PATH type screen_assembly3.py -h . If you have all dependencies then the help menu will display. Otherwise read the erorr and install whichever dependency is missing.

## Running the program

Please see the WIKI

## Authors

* **Liam McIntyre** - *Initial work* - https://github.com/shimbalama/

## License

This project is licensed under the MIT License - see the LICENSE https://github.com/shimbalama/screen_assembly/blob/master/LICENSE file for details

## Acknowledgments

* Mark Davies lab and Jake for testing

