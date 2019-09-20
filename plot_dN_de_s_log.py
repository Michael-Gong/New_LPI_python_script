#%matplotlib inline
#import sdf
import matplotlib
import matplotlib as mpl
#mpl.style.use('https://raw.githubusercontent.com/Michael-Gong/DLA_project/master/style')
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os
from mpl_toolkits.mplot3d import Axes3D
import random
from mpl_toolkits import mplot3d
from matplotlib import rc
import matplotlib.transforms as mtransforms
import sys
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
font = {'family' : 'monospace',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 25,
       }

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')


def pxpy_to_energy(gamma, weight):
    binsize = 300
    en_grid = np.linspace(25,4975,15000)
    en_bin  = np.linspace(0,5000.0,15001)
    en_value = np.zeros_like(en_grid) 
    for i in range(binsize):
#        if i == binsize-1:
#            en_value[i] = sum(weight[en_bin[i]<=gamma])
#        else:
            en_value[i] = sum(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1]) ])
    return (en_grid, en_value)



if __name__ == "__main__":
  part_number = 50000
  #nsteps      = 101 #sum(1 for line in open(from_path+'x_0000.txt'))/part_number


  from_path = './p50000_no_T300/'
  nsteps      = int(sum(1 for line in open(from_path+'t_tot_s.txt'))/part_number)
  t0  = np.loadtxt(from_path+'t_tot_s.txt')/2/np.pi
  #x0  = np.loadtxt(from_path+'x_tot.txt')/2/np.pi
  #y0  = np.loadtxt(from_path+'y_tot.txt')/2/np.pi
  px0 = np.loadtxt(from_path+'px_tot_s.txt')
  py0 = np.loadtxt(from_path+'py_tot_s.txt')
  t0  = np.reshape(t0,(part_number,nsteps))
  #x0  = np.reshape(x0,(part_number,nsteps))
  #y0  = np.reshape(y0,(part_number,nsteps))
  px0 = np.reshape(px0,(part_number,nsteps))
  py0 = np.reshape(py0,(part_number,nsteps))
  gg0 = (px0**2+py0**2+1)**0.5
  ww0 = np.zeros_like(gg0)+1


  from_path = './p50000_rr_T300/'
  t1  = np.loadtxt(from_path+'t_tot_s.txt')/2/np.pi
  #x1  = np.loadtxt(from_path+'x_tot.txt')/2/np.pi
  #y1  = np.loadtxt(from_path+'y_tot.txt')/2/np.pi
  px1 = np.loadtxt(from_path+'px_tot_s.txt')
  py1 = np.loadtxt(from_path+'py_tot_s.txt')
  t1  = np.reshape(t1,(part_number,nsteps))
  #x1  = np.reshape(x1,(part_number,nsteps))
  #y1  = np.reshape(y1,(part_number,nsteps))
  px1 = np.reshape(px1,(part_number,nsteps))
  py1 = np.reshape(py1,(part_number,nsteps))
  gg1 = (px1**2+py1**2+1)**0.5
  ww1 = np.zeros_like(gg1)+1


  from_path = './p50000_qe_T300/'
  to_path   = './jpg_log/'
  t2  = np.loadtxt(from_path+'t_tot_s.txt')/2/np.pi
  #x2  = np.loadtxt(from_path+'x_tot.txt')/2/np.pi
  #y2  = np.loadtxt(from_path+'y_tot.txt')/2/np.pi
  px2 = np.loadtxt(from_path+'px_tot_s.txt')
  py2 = np.loadtxt(from_path+'py_tot_s.txt')
  t2  = np.reshape(t2,(part_number,nsteps))
  #x2  = np.reshape(x1,(part_number,nsteps))
  #y2  = np.reshape(y1,(part_number,nsteps))
  px2 = np.reshape(px2,(part_number,nsteps))
  py2 = np.reshape(py2,(part_number,nsteps))
  gg2 = (px2**2+py2**2+1)**0.5
  ww2 = np.zeros_like(gg2)+1

  for i in range(nsteps):
          width = 100
          kwargs = dict(histtype='stepfilled', alpha=0.3, normed=None, color='r', bins=100, range=(0,25000), weights=ww0[:,i],label='no')
          plt.hist(gg0[:,i], **kwargs)
          kwargs = dict(histtype='stepfilled', alpha=0.3, normed=None, color='b', bins=100, range=(0,25000), weights=ww1[:,i],label='ll')
          plt.hist(gg1[:,i], **kwargs)
          kwargs = dict(histtype='stepfilled', alpha=0.3, normed=None, color='g', bins=100, range=(0,25000), weights=ww2[:,i],label='qe')
          plt.hist(gg2[:,i], **kwargs)
          #### manifesting colorbar, changing label and axis properties ####
          plt.legend(loc='upper right')
          plt.xlabel('$\gamma$',fontdict=font)
          plt.ylabel('dN/dE [A.U.]',fontdict=font)
          plt.xticks(fontsize=20); plt.yticks(fontsize=20);
          plt.yscale('log')
          #plt.xlim(0,4000)
          plt.legend(loc='best',fontsize=20,framealpha=0.5)
          #plt.text(250,6e9,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
          plt.title('t='+str(round(t1[0,i],0))+' $T_0$',fontdict=font)
          fig = plt.gcf()
          fig.set_size_inches(10, 7.)
          fig.savefig(to_path+'dN_dE_comb_s_'+str(i).zfill(4)+'.png',format='png',dpi=80)
          plt.close("all")
          print('plotting '+str(i).zfill(4))
          
  
          
          
