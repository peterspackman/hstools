## Simple command line computational chemistry (SCLCC)
The aim of this project is to provide a simple interface
for a variety of computational techniques associated
with hirshfeld surface calculations in tonto.

This software has only been tested using python 2.7.

A list of contributors is available in contributors.txt
Dependencies:
* python-progressbar2
* fastcluster
* numpy
* scipy
* matplotlib
* docopt

Installing these dependencies is straightforward using pip or easy_install

once the dependencies are installed, the following should compile the c
extensions used in this program:

    cd pack
    make

(there may be warnings about pragmas if your gcc is outdated)
