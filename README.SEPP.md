------------------------------------
Summary
------------------------------------

SEPP stands for `SATe-enabled phylogenetic placement`, and so is a method for the following problem:

- Input: tree `T` and alignment `A` for a set of full-length gene sequences, and set `X` of fragmentary sequences for the same gene

- Output: placement of each fragment in `X` into the tree T, and alignment of each fragment in `X` to the alignment `A`.

SEPP operates by using a divide-and-conquer strategy adopted from SATe [Liu et. al., Science, 2009](http://www.sciencemag.org/content/324/5934/1561.abstract) to improve the alignment produced by running HMMER (code by Sean Eddy). It then places each fragment into the user-provided tree using pplacer (code by Erick Matsen). Our study shows that SEPP provides improved accuracy for quickly evolving genes as compared to other methods.

Developers: Tandy Warnow, Nam Nguyen, and Siavash Mirarab

Publication:
S. Mirarab, N. Nguyen, and T. Warnow, SEPP: SATe-enabled phylogenetic placement, Proceedings of the Pacific Symposium of Biocomputing 2012, pages 247-58 (http://www.ncbi.nlm.nih.gov/pubmed/22174280#).

### Note and Acknowledgment: 
- SEPP bundles the following two programs into its distribution:
  1- pplacer: http://matsen.fhcrc.org/pplacer/
  2- hmmer: http://hmmer.janelia.org/


-------------------------------------
Installation
-------------------------------------
This section details steps for installing and running SEPP. We have run SEPP on Linux and MAC. If you experience difficulty installing or running the software, please contact one of us (Tandy Warnow, Nam Nguyen, or Siavash Mirarab).

Requirements:
-------------------
Before installing the software you need to make sure the following programs are installed on your machine.

1. Python: Version > 2.6. 
2. Java: Version > 1.5

Installation Steps:
-------------------
SEPP is distributed as Python source code. Once you have the above required software installed, do the following. 

1. Obtain the latest SEPP distribution from git repository (using `git clone` or by simply downloading the Zip file). If you downloaded the zip file, uncompress the distribution file.
2. Go to the distribution directory
3. Install: run `sudo python setup.py install`. 
4. Configure: run `sudo python setup.py config`. 

The last step creates a ~/.sepp/ directory and put the default config file under ~/.sepp/main.config. Since this is specific to a user, each user that runs sepp needs to execute the last step. 

Common Problems:
-------------------
1. The last step by default requires root access to the system. If you do not have root access, invoke the setup script as follows: `python setup.py install --prefix=/some/path/on/your/system`, where `/some/path/on/your/system` is the path to a directory on your system to which you do have read and write access. If you use the `--prefix` option, you must ensure that the `lib/python2.x/site-packages` subdirectory (where `x` denotes the minor version number of your Python install) of the directory you specify following `--prefix=` is on Python's search path. To add a directory to Python's search path, modify your PYTHONPATH environment variable.

2. SEPP relies on pplacer and HMMER for alignment and placement steps. These tools are packaged with SEPP. If for some reason the packaged version of HMMER and pplacer do not run in your environment, you need to download and build those programs for your system (see below for links), and point SEPP to them. To point sepp to your installation of hmmer and pllacer modify ~/.sepp/main.config. 
   pplacer: http://matsen.fhcrc.org/pplacer/
   hmmer://hmmer.janelia.org/


---------------------------------------------
Running SEPP
---------------------------------------------
To run SEPP, invoke the `run_sepp.py` script from the `bin` sub-directory of the location in which you installed the Python packages. To see options for running the script, use the command:

`python <bin>/run_sepp.py -h`

The general command for running SEPP is:

`python <bin>/run_sepp.py -t <tree_file> -a <alignment_file> -f <fragment_file> -r <raxml_info_file> -A <alignment_set_size> -P <placement_set_size> `

SEPP can also be run using a configuration file. Sample configuration files and input files can be found under test/unittest/data/simulated/. Change to that directory to run SEPP on the sample files. To run using command options, run

`python <bin>/run_sepp.py -t test.tree -a test.fasta -f test.fas -r test.RAxML_info -A 250 -P 250`

and to run using a configuration file, run

`python <bin>/run_sepp.py -c sample.config`

The main output of SEPP is a .json file, created according to pplacer format. Please refer to pplacer website (currently http://matsen.github.com/pplacer/generated_rst/pplacer.html#json-format-specification) for more information on the format of the josn file. Also note that pplacer package provides a program called guppy that can read .json files and perform downstream steps such as visualization.

In addition to the .json file, SEPP outputs alignments of fragments to reference sets. There could be multiple alignment files created, each corresponding to a different placement subset. 

By setting SEPP_DEBUG environmental variable to `True`, you can instruct SEPP to output more information that can be helpful for debugging.  

---------------------------------------------
Bugs and Errors
---------------------------------------------
SEPP is under active research development at UTCS by the Warnow Lab (and especially with her PhD students Siavash Mirarab and Nam Nguyen). Please report any errors to Siavash Mirarab (smirarab@gmail.com) and Nam Nguyen (namphuon@cs.utexas.edu).

