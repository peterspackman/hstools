#core imports
import re
import sys
import os
import glob
import time
# library imports
from joblib import Parallel, delayed
import progressbar
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.colors import from_levels_and_colors
# local imports
import hist
ticklabelpad = mpl.rcParams['xtick.major.pad']

def readcxsfile(fname):
  """ Hacky way to find the de_vals and di_vals  in a cxs file
      Ideally I'd write a file parser especially for cxs files 
      but this is not yet done
  """
  i = 0
  devals = []
  divals = []

  with open(fname) as f:
    
    content = f.readlines()

    for n in range( len(content) ):
      # I have a feeling this is really inefficient
      # match = re.search('begin d_(.) (\d+)', content[n])

      if content[n].startswith('begin d_e '):
        words = content[n].split()
        for i in range( int(words[2]) ):
          n = n + 1
          devals.append(float(content[n]))

      elif content[n].startswith('begin d_i '):
        words = content[n].split()
        for i in range( int(words[2]) ):
          n = n + 1
          divals.append(float(content[n]))
      
      if devals and divals:
        break
  # We have a problem. i.e. de_vals or di_vals will be empty
  if not devals or not divals:
    print 'Input file is likely missing necessary data from tonto'
    sys.exit(0)

  return devals, divals


def plotfile(x , y, fname = 'out.png', type='linear', nbins=10):
  """ Construct a histogram, plot it, then write the image as a png
  to a given filename.
  """

  # Not sure why i've bothered with linear and log as options
  if(type == 'linear'):
    H , xedges, yedges = hist.bin_data(x,y,bins=nbins)
  else:
    H , xedges, yedges = hist.bin_data_log(x,y,bins=nbins)

  H = np.rot90(H)
  H = np.flipud(H)

  Hmasked = np.ma.masked_where(H==0,H) # Mask 0 pixels

  fig = plt.figure()
  
  # provides a different way to colorize
  # cmap, norm = from_levels_and_colors([1,10,50,100,200],['gray','blue','green','black'])

  plt.pcolormesh(xedges,yedges,H,norm=mpl.colors.LogNorm())
  plt.xlim([0.5,2.5])
  plt.ylim([0.5,2.5])
  plt.xticks(np.arange(0.5,2.5,0.2))
  plt.yticks(np.arange(0.5,2.5,0.2))
  plt.grid()

  plt.annotate(r'$d_e$',fontsize=20,xy=(1,0),xytext=(5,-ticklabelpad),
               ha='left',va='top', xycoords='axes fraction',textcoords='offset points')
  plt.annotate(r'$d_i$',fontsize=20,xy=(0,1.02),xytext=(5,-ticklabelpad), 
               ha='right',va='bottom', xycoords='axes fraction',textcoords='offset points')
 # print 'Saving histogram plot to {0} bin size:{1}'.format(fname,nbins)
  fig.savefig(fname,bbox_inches='tight') 
  plt.clf()
  plt.close()

def process_file(fname,resolution=10,write_png=False):
  """ Read a file from fname, generate a histogram and potentially write
      the png of it to file
  """
  x,y = readcxsfile(fname)
  cname =  os.path.basename(os.path.splitext(fname)[0])
  h = hist.bin_data(x,y,resolution)


  if(write_png):
    outfile = os.path.splitext(fname)[0] + '{0}bins.png'.format(resolution)
    plotfile(x,y,fname=outfile,nbins=resolution) 
  
  return h, cname



def batch_process(dirname, suffix='.cxs', resolution=10,
                  write_png=False,threads=4):
  """Generate n histograms from a directory, returning a list of them
     and their corresponding substance names

     threads would be more accurately called 'processes'
  """

  files = glob.glob(os.path.join(dirname,'*'+suffix))

  histograms = []
  names = []
  
  start_time = time.time()  
  
  vals = Parallel(n_jobs=threads,verbose=3)(
               delayed(process_file)(f, resolution=resolution,write_png=write_png)
               for f in files)
  
  # unzip the output
  histograms,names = zip(*vals)
  
  print 'Reading {0} files took {1} seconds using {2} threads.'.format(len(files), time.time() - start_time,threads)
 
  return (histograms, names)


