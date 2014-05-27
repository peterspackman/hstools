import re
import sys
import hist
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt

# Hacky way to find the de_vals etc in a cxs file
def readcxsfile(fname):
  i = 0
  devals = []
  divals = []
  matches = 0

  with open(fname) as f:
    
    content = f.readlines()

    for n in range( len(content) ):
      match = re.search('begin d_(.) (\d+)', content[n])

      if match:
        matches += 1

        for i in range( int(match.group(2)) ):

          n = n + 1

          if(match.group(1) == 'e'):
            devals.append(float(content[n]))
          else:
            divals.append(float(content[n]))

  if(matches < 2):
    print 'Input file is likely missing necessary data from tonto'
    sys.exit(0)

  return devals, divals

# Plot the data to file

def plotfile(x , y, fname = 'out.png', type='linear', nbins=10):

  if(type == 'linear'):
    H , xedges, yedges = hist.bin_data(x,y,bins=nbins)
  else:
    H , xedges, yedges = hist.bin_data_log(x,y,bins=nbins)

  H = np.rot90(H)
  H = np.flipud(H)

  Hmasked = np.ma.masked_where(H==0,H) # Mask 0 pixels

  fig = plt.figure()

  plt.pcolormesh(xedges,yedges,H,norm=mpl.colors.LogNorm())
  plt.axis("equal")
  fig.savefig(fname,bbox_inches='tight')
