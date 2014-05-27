import sys
import os
import shutil

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

if __name__ == '__main__':
  if(needs_splitting(sys.argv[1])):
    split_cif(sys.argv[1])
  else:
    print "Passing on {0}, doesn't need splitting".format(sys.argv[1])
