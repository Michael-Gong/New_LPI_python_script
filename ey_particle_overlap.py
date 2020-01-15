#!/public/home/users/bio001/tools/python-2.7.11/bin/python
import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
import matplotlib.colors as mcolors 
import scipy.ndimage as ndimage
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
import multiprocessing as mp
  
if __name__ == "__main__":
  print ('This is main of module "test2d.py"')
  ######## Constant defined here ########
  pi        =     3.1415926535897932384626
  q0        =     1.602176565e-19 # C
  m0        =     9.10938291e-31  # kg
  v0        =     2.99792458e8    # m/s^2
  kb        =     1.3806488e-23   # J/K
  mu0       =     4.0e-7*np.pi       # N/A^2
  epsilon0  =     8.8541878176203899e-12 # F/m
  h_planck  =     6.62606957e-34  # J s
  wavelength=     0.8e-6
  frequency =     v0*2*pi/wavelength
  
  exunit    =     m0*v0*frequency/q0
  bxunit    =     m0*frequency/q0
  denunit    =     frequency**2*epsilon0*m0/q0**2
  jalf      =     4*np.pi*epsilon0*m0*v0**3/q0/wavelength**2
  print('electric field unit: '+str(exunit))
  print('magnetic field unit: '+str(bxunit))
  print('density unit nc: '+str(denunit))
  
  font = {'family' : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 25,  
          }  
  font2 = {'family' : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 8,  
          } 
  space_1 = 1
  space_2 = 1
  font_size = 25
  marker_size=0.05
##below is for generating mid transparent colorbar
  c_red = matplotlib.colors.colorConverter.to_rgba('red')
  c_blue= matplotlib.colors.colorConverter.to_rgba('blue')
  c_yellow= matplotlib.colors.colorConverter.to_rgba('yellow')
  c_cyan= matplotlib.colors.colorConverter.to_rgba('cyan')
  c_white_trans = matplotlib.colors.colorConverter.to_rgba('white',alpha = 0.0)
  cmap_rb = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_red,c_white_trans,c_blue],128) 
  cmap_br = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_blue,c_white_trans,c_red],128)
  cmap_ryb= matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_red,c_yellow,c_blue],128)
  cmap_yc = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_yellow,c_white_trans,c_cyan],128) 

  c_brown = matplotlib.colors.colorConverter.to_rgba('gold')
  c_green = matplotlib.colors.colorConverter.to_rgba('springgreen')
  cmap_bg = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_brown,c_white_trans,c_white_trans,c_green],128)
  c_black = matplotlib.colors.colorConverter.to_rgba('black')
  c_white = matplotlib.colors.colorConverter.to_rgba('white')
  cmap_bw = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_white_trans,c_white,c_black],128)
   
##end for transparent colorbar##
 
##below is for norm colorbar
  class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y)) 
##end for norm colorbar####

  color_list = ['blue','limegreen','red'] 
  
def processplot(n): 
    #from_path = './new_single_a30/'
    #to_path   = './new_single_a30_fig/'
    from_path = './PW_w020/'
    to_path   = './PW_w020_fig/'
    ######### Script code drawing figure ################
    data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']
    x  = data['Grid/Grid_mid'].data[0]/1.0e-6
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    X, Y = np.meshgrid(x, y) 
    
    fig=plt.subplot(1,1,1)
    name = 'ey'
    ex = data['Electric Field/'+str.capitalize(name)].data/exunit
    eee = 30.0
    levels = np.linspace(-eee, eee, 41)
    ex.T[ex.T < -eee]=-eee
    ex.T[ex.T >  eee]= eee
    plt.contourf(X, Y, ex.T, levels=levels, cmap=cmap_yc)
    #### manifesting colorbar, changing label and axis properties ####
    cbar=plt.colorbar(pad=0.01,ticks=np.linspace(-eee, eee, 5))
    cbar.ax.tick_params(labelsize=font2['size']) 
    #cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=font2['size'])        
    cbar.set_label('$E_y$ [$m_ec\omega_0/|e|$]',fontdict=font2)        

    if 'Grid/Particles/electron' in data.keys():
          ion_x = data['Grid/Particles/electron'].data[0]/1e-6
          ion_y = data['Grid/Particles/electron'].data[1]/1e-6
#          ion_px = data['Particles/Px/electron'].data/(1000*m0*v0)
          plt.scatter(ion_x[::space_2], ion_y[::space_2], c=color_list[0], s=marker_size, marker='.',alpha=0.8,label='electron',zorder=1,lw=0)
    if 'Grid/Particles/carbon' in data.keys():
          ion_x = data['Grid/Particles/carbon'].data[0]/1e-6
          ion_y = data['Grid/Particles/carbon'].data[1]/1e-6
          plt.scatter(ion_x[::space_2], ion_y[::space_2], c=color_list[1], s=marker_size, marker='.',alpha=0.8,label='carbon',zorder=2,lw=0)
    if 'Grid/Particles/proton' in data.keys():
          ion_x = data['Grid/Particles/proton'].data[0]/1e-6
          ion_y = data['Grid/Particles/proton'].data[1]/1e-6
#          ion_px = data['Particles/Px/proton'].data/(1836*m0*v0)
          plt.scatter(ion_x[::space_2], ion_y[::space_2], c=color_list[2], s=marker_size, marker='.',alpha=0.8,label='proton',zorder=3,lw=0)


    plt.text(2.,10,str(round(time/1.0e-15,0))+' fs',fontdict=font2,color='k')
    #ax.text(21.,1.75,'t = 70 fs',fontdict=font)
    
    plt.ylim(-12,12)
    plt.xlim(-4,12)
    plt.xlabel('$x$ [$\mu m$]',fontdict=font)
    plt.ylabel('$y$ [$\mu m$]',fontdict=font)
    plt.xticks(fontsize=font_size); 
    plt.yticks(fontsize=font_size);
    plt.legend(loc='best',fontsize=font_size,framealpha=0.5)
    plt.subplots_adjust(left=0.12, bottom=0.14, right=0.98, top=0.97, wspace=0.011, hspace=0.051)

    fig = plt.gcf()
    fig.set_size_inches(10, 6.)
    fig.savefig(to_path+'ey_part_'+str(n).zfill(4)+'.png',format='png',dpi=160)
    plt.close("all")
    print(to_path+'ey_part_'+str(n).zfill(4)+'.png')
    return 0

if __name__ == '__main__':
  start   =  1 # start time
  stop    =  19  # end time
  step    =  1  # the interval or step
    
  inputs = range(start,stop+step,step)
  pool = mp.Pool(processes=1)
  results = pool.map(processplot,inputs)
  print(results)
