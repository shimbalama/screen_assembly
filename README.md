# Screen assemblies

Pipeline that screens for presence of genes of interest (GOI) in bacterial assemblies. Generates multiple CSVs and plots that descirbe which genes are present and how variabel their sequence is. Can you DNA or protein query sequences (GOIs) and DNA contigs/fastas or proein fastas as database (db) to search in. 

## Getting Started

You need one fasta file with all GOIs as the query and a folder with db contigs/fastas. Db files can only have one '.' in the (i.e., sample_1.fa NOT sample.1.fa) 

### Prerequisites

#### Required

Python 3 and scypi/Biopython
Command line blast

#### Optional

Clustal Omega
RAXML and or IQtree

### Installing

* Download the script and place it in your PATH: 
  * git clone https://github.com/shimbalama/screen_assembly.git
  * export PATH="your_path:$PATH"
* Make sure its executable (chmod +x).
* Download lab_modules.py git clone https://github.com/shimbalama/common_modules.git
* Use a text editor to set this line in screen_assembly3.py to point at the dir you put lab_modules.py in:  sys.path.append('/Users/lmcintyre/Dropbox/work/uniMelb/code/github/common_modules') becomes sys.path.append('your_path/common_modules')

## Running the tests

Once screen_assembly3.py is in your PATH type screen_assembly3.py -h . If you have all dependencies then the menu will display. Otherwise read the erorr and install whichever dependency is missing.


## Authors

* **Liam McIntyre** - *Initial work* - https://github.com/shimbalama/

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Mark Davies lab and Jake for testing

