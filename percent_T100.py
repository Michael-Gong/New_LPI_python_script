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

font_2 = {'family' : 'monospace',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 20,
       }


font_size = 25
font_size2 = 20
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

upper = matplotlib.cm.jet(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolor_jet = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

upper = matplotlib.cm.rainbow(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolor_rainbow = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

def pxpy_to_energy(gamma, weight):
    binsize = 200
    en_grid = np.linspace(50,19950,200)
    en_bin  = np.linspace(0,20000.0,201)
    en_value = np.zeros_like(en_grid) 
    for i in range(binsize):
#        if i == binsize-1:
#            en_value[i] = sum(weight[en_bin[i]<=gamma])
#        else:
            en_value[i] = sum(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1]) ])
    return (en_grid, en_value)



if __name__ == "__main__":
  part_number = 50000
  from_path = './p50000_rr_T150/'
  nsteps      = int(sum(1 for line in open(from_path+'t_tot_s.txt'))/part_number)

  from_path = './p50000_qe_T150/'
  t2  = np.loadtxt(from_path+'t_tot_s.txt')/2/np.pi
  #x2  = np.loadtxt(from_path+'x_tot_s.txt')/2/np.pi
  #y2  = np.loadtxt(from_path+'y_tot_s.txt')/2/np.pi
  px2 = np.loadtxt(from_path+'px_tot_s.txt')
  py2 = np.loadtxt(from_path+'py_tot_s.txt')
  t2  = np.reshape(t2,(part_number,nsteps))
  #x2  = np.reshape(x1,(part_number,nsteps))
  #y2  = np.reshape(y1,(part_number,nsteps))
  px2 = np.reshape(px2,(part_number,nsteps))
  py2 = np.reshape(py2,(part_number,nsteps))
  gg2 = (px2**2+py2**2+1)**0.5
  ww2 = np.zeros_like(gg2)+1

  from_path = './p50000_rr_T150/'
  t1  = np.loadtxt(from_path+'t_tot_s.txt')/2/np.pi
  #x1  = np.loadtxt(from_path+'x_tot_s.txt')/2/np.pi
  #y1  = np.loadtxt(from_path+'y_tot_s.txt')/2/np.pi
  px1 = np.loadtxt(from_path+'px_tot_s.txt')
  py1 = np.loadtxt(from_path+'py_tot_s.txt')
  t1  = np.reshape(t1,(part_number,nsteps))
  #x1  = np.reshape(x1,(part_number,nsteps))
  #y1  = np.reshape(y1,(part_number,nsteps))
  px1 = np.reshape(px1,(part_number,nsteps))
  py1 = np.reshape(py1,(part_number,nsteps))
  gg1 = (px1**2+py1**2+1)**0.5
  ww1 = np.zeros_like(gg1)+1

  axis_time = np.linspace(0,nsteps*3.333,nsteps)
 
  per_001_qe = np.zeros_like(axis_time)   
  per_002_qe = np.zeros_like(axis_time)   
  per_005_qe = np.zeros_like(axis_time)   
  per_010_qe = np.zeros_like(axis_time)   
  per_020_qe = np.zeros_like(axis_time)   
  per_050_qe = np.zeros_like(axis_time)   
  per_100_qe = np.zeros_like(axis_time)   
 
  per_001_rr = np.zeros_like(axis_time)   
  per_002_rr = np.zeros_like(axis_time)   
  per_005_rr = np.zeros_like(axis_time)   
  per_010_rr = np.zeros_like(axis_time)   
  per_020_rr = np.zeros_like(axis_time)   
  per_050_rr = np.zeros_like(axis_time)   
  per_100_rr = np.zeros_like(axis_time)   

  for i in range(nsteps):
      sort_en    = np.sort(gg2[:,i])
      per_001_qe[i] = np.sum(sort_en[int(part_number*0.99):])/np.size(sort_en[int(part_number*0.99):])  
      per_002_qe[i] = np.sum(sort_en[int(part_number*0.98):])/np.size(sort_en[int(part_number*0.98):])  
      per_005_qe[i] = np.sum(sort_en[int(part_number*0.95):])/np.size(sort_en[int(part_number*0.95):]) 
      per_010_qe[i] = np.sum(sort_en[int(part_number*0.90):])/np.size(sort_en[int(part_number*0.90):])
      per_020_qe[i] = np.sum(sort_en[int(part_number*0.80):])/np.size(sort_en[int(part_number*0.80):])
      per_050_qe[i] = np.sum(sort_en[int(part_number*0.50):])/np.size(sort_en[int(part_number*0.50):])
      per_100_qe[i] = np.sum(sort_en[int(part_number*0.00):])/np.size(sort_en[int(part_number*0.00):])  

      sort_en    = np.sort(gg1[:,i])
      per_001_rr[i] = np.sum(sort_en[int(part_number*0.99):])/np.size(sort_en[int(part_number*0.99):])  
      per_002_rr[i] = np.sum(sort_en[int(part_number*0.98):])/np.size(sort_en[int(part_number*0.98):])  
      per_005_rr[i] = np.sum(sort_en[int(part_number*0.95):])/np.size(sort_en[int(part_number*0.95):]) 
      per_010_rr[i] = np.sum(sort_en[int(part_number*0.90):])/np.size(sort_en[int(part_number*0.90):])
      per_020_rr[i] = np.sum(sort_en[int(part_number*0.80):])/np.size(sort_en[int(part_number*0.80):])
      per_050_rr[i] = np.sum(sort_en[int(part_number*0.50):])/np.size(sort_en[int(part_number*0.50):])
      per_100_rr[i] = np.sum(sort_en[int(part_number*0.00):])/np.size(sort_en[int(part_number*0.00):])
 
      print('finish sort: '+str(round(100.*i/nsteps,1))+'%')

  plt.subplot(1,2,1)        
  plt.plot(axis_time, per_001_qe*0.51e-3, '-', color='tomato',         linewidth=3,label='QED 1%')
  plt.plot(axis_time, per_010_qe*0.51e-3, '-', color='mediumseagreen',      linewidth=3,label='QED 10%')
  plt.plot(axis_time, per_100_qe*0.51e-3, '-', color='royalblue', linewidth=3,label='QED 100%')

  plt.plot(axis_time, per_001_rr*0.51e-3, ':', color='tomato',         linewidth=3,label='Classic 1%')
  plt.plot(axis_time, per_010_rr*0.51e-3, ':', color='mediumseagreen',      linewidth=3,label='Classic 10%')
  plt.plot(axis_time, per_100_rr*0.51e-3, ':', color='royalblue', linewidth=3,label='Classic 100%')
#### manifesting colorbar, changing label and axis properties ####
  plt.xlabel('t [fs]',fontdict=font)
  plt.ylabel(r'$\overline{\epsilon}_e$'+' [GeV]',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks(fontsize=font_size);
#  plt.title('dN/dE for QED RR', fontsize=font_size)
  plt.xlim(0,500)
  plt.ylim(0,5.1)
#  plt.yscale('log')
  plt.legend(loc='best',fontsize=font_size2)
#  plt.ylim(0,500)
  plt.grid()
#plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)
#  plt.text(3.,2,'t = 500 fs',fontdict=font)

  plt.subplot(1,2,2)
  axis_x = ['100','50','20','10','5','2','1'] 
  axis_y_qe = np.array([per_100_qe[-12], per_050_qe[-12], per_020_qe[-12], per_010_qe[-12], per_005_qe[-12], per_002_qe[-12], per_001_qe[-12]])*0.51e-3
  axis_y_rr = np.array([per_100_rr[-12], per_050_rr[-12], per_020_rr[-12], per_010_rr[-12], per_005_rr[-12], per_002_rr[-12], per_001_rr[-12]])*0.51e-3
  plt.scatter(axis_x, axis_y_qe, c='tomato',marker='o',s=250, label='QED', edgecolors='black', linewidth='3',alpha=1,zorder=2)
  plt.scatter(axis_x, axis_y_rr, c='mediumseagreen',marker='o',s=250, label='Classic', edgecolors='black', linewidth='3',alpha=1,zorder=2)
#### manifesting colorbar, changing label and axis properties ####
  plt.xlabel('$\kappa$ [%]',fontdict=font)
  plt.ylabel(r'$\overline{\epsilon}_e$'+' [GeV]',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks(fontsize=font_size);
#  plt.title('dN/dE for QED RR', fontsize=font_size)
#  plt.xlim(0,500)
#  plt.yscale('log')
  plt.ylim(0,5.1)
  plt.legend(loc='best',fontsize=font_size2)
#  plt.ylim(0,500)
  plt.grid()
#plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)
#  plt.text(3.,2,'t = 500 fs',fontdict=font)




  plt.subplots_adjust(top=0.98, bottom=0.14, left=0.1, right=0.98, hspace=0.20, wspace=0.23)

  fig = plt.gcf()
  fig.set_size_inches(14, 6)
#fig.set_size_inches(5, 4.5)
  fig.savefig('./wrap_percent.png',format='png',dpi=160)
  plt.close("all")

  
          
          
