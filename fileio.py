#core imports
import re
import sys
import os
import glob
import time
# library imports
import progressbar
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.colors import from_levels_and_colors
# local imports
import hist


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
  plt.axis("equal")
  print 'Saving histogram plot to {0} bin size:{1}'.format(fname,nbins)
  fig.savefig(fname,bbox_inches='tight')
  
  plt.clf()
  plt.close()

def process_file(fname,resolution=10,write_png=False):
  """ Read a file from fname, generate a histogram and potentially write
      the png of it to file
  """
  x,y = readcxsfile(fname)
  h = hist.bin_data(x,y,resolution)
  cname =  os.path.basename(os.path.splitext(fname)[0])

  if(write_png):
    outfile = os.path.splitext(fname)[0] + '{0}bins.png'.format(resolution)
    plotfile(x,y,fname=outfile,nbins=resolution) 
  
  return h, cname



def batch_process(dirname, suffix='.cxs', resolution=10,
                  write_png=False):
  """Generate n histograms from a directory, returning a list of them
  """

  files = glob.glob(os.path.join(dirname,'*'+suffix))

  histograms = []
  names = []
  widgets = ['Reading files:',' ', 
              progressbar.Bar(marker=unichr(0x2592)),
              ' ',progressbar.ETA()]
  pbar = progressbar.ProgressBar(widgets = widgets, maxval = len(files))
  pbar.start()
  start_time = time.time()  
  for i,f in enumerate(files):
    
    h, cname = process_file(f,resolution=resolution,write_png =write_png) 
    
    #Can add bounds = true to make the bins the same values throughout histograms
    histograms.append(h)
    names.append(cname)
    pbar.update(i)
  
  pbar.finish()
  print 'Reading {0} files took {1} seconds'.format(len(files), time.time() - start_time)
 
  return (histograms, names)


