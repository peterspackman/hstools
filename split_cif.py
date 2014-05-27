import sys
import argparse
import os
import shutil
import glob

def split_cif(fname):  
  dir_name = os.path.dirname(fname)

  with open(fname) as f:
    lines_out = []
    outname = ''
    for line in f:
      
      if line.startswith('data_'):
        if(outname):
          with open(outname,'w') as of:
            for item in lines_out:
              of.write(item)
          lines_out = []

        outname = dir_name +'/'+ line.strip().replace('data_','')+'.cif'  
      lines_out.append(line)

  with open(outname,'w') as of:
   for item in lines_out:
     of.write(item)
  

  processed_dir = dir_name + '/processed_cifs'
  if not os.path.exists(processed_dir):
    os.makedirs(processed_dir)
  shutil.move(fname,processed_dir+'/'+os.path.basename(fname))

def needs_splitting(fname):
  count = 0
  with open(fname) as f:
    for line in f:
      if line.startswith('data_'):
        count += 1
  return count > 1

def batch_split(dname):
  files = os.path.join(dname,'*.cif')

  for f in glob.glob(files):
    if(needs_splitting(f)):
      split_cif(f)
    else:
      print "Passing on {0}, doesn't need splitting".format(f)

def cl_options():
  parser = argparse.ArgumentParser(description='Split cif files into constituent molecules')
  parser.add_argument('-f','--file',help='the file to be processed')
  parser.add_argument('-d','--dir',help='A directory of files to be processed')

  return parser.parse_args()

def main():
  opts = cl_options()

  if opts.dir:
    print 'Processing cif files in '+opts.dir
    batch_split(opts.dir)
    sys.exit(0)

  if opts.file:
    print 'Processing '+opts.file
    
    if(needs_splitting(opts.file)):
      split_cif(opts.file)
    else:
      print "Passing on {0}, doesn't need splitting".format(opts.file)
    
  sys.exit(0)

if __name__ == '__main__':
  main()
