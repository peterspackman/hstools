#!/usr/bin/python
import sys
import hist
import calc
import fileio as fio
import argparse
import os
import glob
import time


# Generate n histograms from a directory, then do something with the list of them
def batch_process(dirname, suffix='.cxs', resolution=10,do_calc=calc.get_correl_mat):
  
  files = os.path.join(dirname,'*'+suffix)
  histograms = []

  for f in glob.glob(files):
    #TIMING
    start_time = time.time()
    #END TIMING

    x,y = fio.readcxsfile(f)
   
    
    #TIMING
    print 'Reading {0} took {1} seconds'.format(f, time.time() - start_time)
    #END TIMING 

    histograms.append(hist.bin_data(x,y,resolution))

  print do_calc(histograms)




def cl_options():
  parser = argparse.ArgumentParser(description='Calculate histograms')
  
  parser.add_argument('-f','--file',metavar='fname', 
                      help='the file to be processed')
  parser.add_argument('-b','--bins',metavar='nbins', type=int,
                      help='the number of bins to use in the histogram')
  parser.add_argument('-d','--dir',
                      help='A directory of cxs files to be read')
  return parser.parse_args()

def main():
  opts = cl_options()
  
  batch_process(opts.dir)

 
if __name__ == '__main__':
  main()

