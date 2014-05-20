#!/usr/bin/python
import sys
import hist
import calc
import fileio as fio
import argparse



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
  print opts

 
if __name__ == '__main__':
  main()

