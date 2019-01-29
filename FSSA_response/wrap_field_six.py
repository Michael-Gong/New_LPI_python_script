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
  wavelength=     1.0e-6
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
          'size'   : 12,  
          } 

  font_size = 25 
##below is for generating mid transparent colorbar
  c_red = matplotlib.colors.colorConverter.to_rgba('red')
  c_blue= matplotlib.colors.colorConverter.to_rgba('blue')
  c_white_trans = matplotlib.colors.colorConverter.to_rgba('white',alpha = 0.0)
  cmap_rb = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_red,c_white_trans,c_blue],128) 
  cmap_br = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_blue,c_white_trans,c_red],128)

  c_brown = matplotlib.colors.colorConverter.to_rgba('saddlebrown')
  c_green = matplotlib.colors.colorConverter.to_rgba('seagreen')
  cmap_bg = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_brown,c_white_trans,c_green],128)
  
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

 
  
if __name__ == '__main__':
  ######### Parameter you should set ###########
  #start   =  210  # start time
  #stop    =  210  # end time
  #step    =  1  # the interval or step

  
#  youwant = ['electron_x_px','electron_density','electron_en','electron_theta_en','ey'] #,'electron_ekbar']
  #youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
  #youwant Derived electron_density,electron_ekbar...
  #youwant dist_fn electron_x_px, electron_py_pz, electron_theta_en...
  #if (os.path.isdir('jpg') == False):
  #  os.mkdir('jpg')
    #from_path = './uniform_a190_n30/'  
    n1   = 4
    dim  =  '2d'
    from_path = './one_double/' # './one_20/'
    to_path   = from_path
    ######### Script code drawing figure ################
    #for n in range(start,stop+step,step):
    #### header data ####
    data = sdf.read(from_path+'q'+str(n1).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']
    x  = data['Grid/Grid_mid'].data[0]/1.0e-6
    print('ok')
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    X, Y = np.meshgrid(x, y) 
    
    den = (data['Derived/Number_Density/Electron'].data+data['Derived/Number_Density/E_1'].data)/denunit
    if dim == '3d':
        n3d = len(den[0,0,:])
        den = (den[:,:,n3d//2-1]+den[:,:,n3d//2])/2.0
    eee = 20.0
    levels = np.linspace(0, eee, 101)
    den.T[den.T > eee]=eee 

    den_a = (data['Derived/Number_Density_averaged/Electron'].data+data['Derived/Number_Density_averaged/E_1'].data)/denunit
    if dim == '3d':
        n3d = len(den_a[0,0,:])
        den_a = (den_a[:,:,n3d//2-1]+den_a[:,:,n3d//2])/2.0
    eee = 20.0
    levels = np.linspace(0, eee, 101)
    den_a.T[den_a.T > eee]=eee 

    #gs = gridspec.GridSpec(2, 2, width_ratios=[6, 1], height_ratios=[1, 3])

    
    #ax=plt.subplot(3,2,1)
    ax=plt.subplot2grid((3, 3), (0, 0))
    #axin1 = inset_axes(ax, width='15%', height='5%', loc='upper left',pad=0.2)
    #axin2 = inset_axes(ax, width='15%', height='5%', loc='lower left',pad=0.2)
    #axin1 = inset_axes(ax,width="5%",height="45%",loc='upper right', bbox_to_anchor=(1.05, 0., 1, 1), bbox_transform=ax.transAxes, borderpad=0,)
    #axin2 = inset_axes(ax,width="5%",height="45%",loc='lower right', bbox_to_anchor=(1.05, 0., 1, 1), bbox_transform=ax.transAxes, borderpad=0,)
    #### manifesting colorbar, changing label and axis properties ####
    image1=ax.contourf(X, Y, den.T, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap='pink_r')
    ax.contour(X, Y, den_a.T, levels=[0.9], colors='limegreen', linewidths=3.5,origin='lower')
#    cbar=plt.colorbar(image1,pad=0.01,ticks=np.linspace(0.0, eee, 5),orientation="vertical")
#    cbar.set_label('$n_e$ [$n_c$]', fontdict=font2)
#    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=font2['size'])        

    name = 'ey'
    data = sdf.read(from_path+'e_fields'+str(n1).zfill(4)+".sdf",dict=True)
    ex = data['Electric Field/'+str.capitalize(name)].data/exunit
    if dim == '3d':
        ex = (ex[:,:,n3d//2-1]+ex[:,:,n3d//2])/2.0
    eee = 50.0
    levels = np.linspace(-eee, eee, 51)
    ex.T[ex.T < -eee]=-eee
    ex.T[ex.T >  eee]= eee
    image2=ax.contourf(X, Y, ex.T, levels=levels, cmap=cmap_br)
    #### manifesting colorbar, changing label and axis properties ####
#    cbar=plt.colorbar(image2,pad=0.01,ticks=np.linspace(-eee, eee, 5),orientation="vertical")
#    cbar.set_label(r'$E_y\ [m_ec\omega_0/e]$',fontdict=font2)        
#    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=font2['size'])        
#    plt.title('Ey_den_e at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
    ax.text(2.,4.5,'t = '+str(round(time/1.0e-15,0))+' fs',fontdict=font,color='white')
    #ax.text(21.,1.75,'t = 70 fs',fontdict=font)
    ax.set_ylim(-6.5,6.5)
    ax.set_xlim(-5,15)
#    ax.set_xlabel('X [$\lambda$]',fontdict=font)
    ax.set_ylabel('Y [$\mu m$]',fontdict=font)
    ax.tick_params(axis='y',labelsize=25) 
    ax.tick_params(axis='x',labelsize=0.0) 
#    ax.tick_params(axis='both',labelsize=25) 
    ax.set_xticks(ticks=[0,10])
    ax.set_yticks(ticks=[-5,0,5])


    #ax=plt.subplot(3,2,3)
    ax=plt.subplot2grid((3, 3), (1, 0))
    name = 'Bz_averaged'
    data = sdf.read(from_path+'b_fields'+str(n1).zfill(4)+".sdf",dict=True)
    ex = data['Magnetic Field/'+str.capitalize(name)].data/bxunit
    if dim == '3d':
        ex = (ex[:,:,n3d//2-1]+ex[:,:,n3d//2])/2.0
    eee = 10.0
    levels = np.linspace(-eee, eee, 51)
    ex.T[ex.T < -eee]=-eee
    ex.T[ex.T >  eee]= eee
    image2=ax.contourf(X, Y, ex.T, levels=levels, cmap=cmap_bg)#'BrBG')
    #### manifesting colorbar, changing label and axis properties ####
#    cbar=plt.colorbar(image2,pad=0.01,ticks=np.linspace(-eee, eee, 5),orientation="vertical")
#    cbar.set_label(r'$B_z\ [m_ec\omega_0/e]$',fontdict=font2)        
#    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=font2['size'])        
#    plt.title('Bz_averaged at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
    #ax.text(21.,1.75,'t = '+str(round(time/1.0e-15,0))+' fs',fontdict=font)
    #ax.text(21.,1.75,'t = 70 fs',fontdict=font)
#    ax.contour(X, Y, den_a.T, levels=[0.9], colors='limegreen', linewidths=2,origin='lower')
    image1=ax.contour(X, Y, den_a.T, levels=[19.5], colors='k', linewidths=2,origin='lower')
    
    ax.set_ylim(-6.5,6.5)
    ax.set_xlim(-5,15)
#    ax.set_xlabel('X [$\mu m$]',fontdict=font)
    ax.set_ylabel('Y [$\mu m$]',fontdict=font)
    ax.tick_params(axis='y',labelsize=25) 
    ax.tick_params(axis='x',labelsize=0) 
    ax.set_xticks(ticks=[0,10])
    ax.set_yticks(ticks=[-5,0,5])


    #ax=plt.subplot(3,2,5)
    ax=plt.subplot2grid((3, 3), (2, 0))
    name = 'Ey_averaged'
    data = sdf.read(from_path+'e_fields'+str(n1).zfill(4)+".sdf",dict=True)
    ex = data['Electric Field/'+str.capitalize(name)].data/exunit
    if dim == '3d':
        ex = (ex[:,:,n3d//2-1]+ex[:,:,n3d//2])/2.0
    eee = 10.0
    levels = np.linspace(-eee, eee, 51)
    ex.T[ex.T < -eee]=-eee
    ex.T[ex.T >  eee]= eee
    image2=ax.contourf(X, Y, ex.T, levels=levels, cmap=cmap_bg)#'BrBG')
    #### manifesting colorbar, changing label and axis properties ####
#    cbar=plt.colorbar(image2,pad=0.01,ticks=np.linspace(-eee, eee, 5),orientation="vertical")
#    cbar.set_label(r'$E_y\ [m_ec\omega_0/e]$',fontdict=font2)        
#    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=font2['size'])        
#    plt.title('Ey_averaged at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
    #ax.text(21.,1.75,'t = '+str(round(time/1.0e-15,0))+' fs',fontdict=font)
    #ax.text(21.,1.75,'t = 70 fs',fontdict=font)
#    ax.contour(X, Y, den_a.T, levels=[0.9], colors='limegreen', linewidths=2,origin='lower')
    image1=ax.contour(X, Y, den_a.T, levels=[19.5], colors='k', linewidths=2,origin='lower')
    
    ax.set_ylim(-6.5,6.5)
    ax.set_xlim(-5,15)
    ax.set_xlabel('X [$\mu m$]',fontdict=font)
    ax.set_ylabel('Y [$\mu m$]',fontdict=font)
    ax.tick_params(axis='y',labelsize=25) 
    ax.tick_params(axis='x',labelsize=25) 
    ax.set_xticks(ticks=[0,10])
    ax.set_yticks(ticks=[-5,0,5])



    n2   = 13
    data = sdf.read(from_path+'q'+str(n2).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']
    x  = data['Grid/Grid_mid'].data[0]/1.0e-6
    print('ok')
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    X, Y = np.meshgrid(x, y) 
    
    den = (data['Derived/Number_Density/Electron'].data+data['Derived/Number_Density/E_1'].data)/denunit
    if dim == '3d':
        n3d = len(den[0,0,:])
        den = (den[:,:,n3d//2-1]+den[:,:,n3d//2])/2.0
    eee = 20.0
    levels = np.linspace(0, eee, 101)
    den.T[den.T > eee]=eee 

    den_a = (data['Derived/Number_Density_averaged/Electron'].data+data['Derived/Number_Density_averaged/E_1'].data)/denunit
    if dim == '3d':
        n3d = len(den_a[0,0,:])
        den_a = (den_a[:,:,n3d//2-1]+den_a[:,:,n3d//2])/2.0
    eee = 20.0
    #levels = np.linspace(0, eee, 101)
    den_a.T[Y > 4]=20.0

    #gs = gridspec.GridSpec(2, 2, width_ratios=[6, 1], height_ratios=[1, 3])

    
    #ax=plt.subplot(3,2,2)
    ax=plt.subplot2grid((3, 3), (0, 1),colspan=2)
    #axin1 = inset_axes(ax, width='15%', height='5%', loc='upper left',pad=0.2)
    #axin2 = inset_axes(ax, width='15%', height='5%', loc='lower left',pad=0.2)
    #axin1 = inset_axes(ax,width="5%",height="45%",loc='upper right', bbox_to_anchor=(1.05, 0., 1, 1), bbox_transform=ax.transAxes, borderpad=0,)
    #axin2 = inset_axes(ax,width="5%",height="45%",loc='lower right', bbox_to_anchor=(1.05, 0., 1, 1), bbox_transform=ax.transAxes, borderpad=0,)
    #### manifesting colorbar, changing label and axis properties ####
    image1=ax.contourf(X, Y, den.T, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap='pink_r')
    #ax.contour(X, Y, den_a.T, levels=[0.9], colors='limegreen', linewidths=2,origin='lower')
#    cbar=plt.colorbar(image1,pad=0.01,ticks=np.linspace(0.0, eee, 5),orientation="vertical")
#    cbar.set_label('$n_e$ [$n_c$]', fontdict=font2)
#    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=font2['size'])        

    name = 'ey'
    data = sdf.read(from_path+'e_fields'+str(n2).zfill(4)+".sdf",dict=True)
    ex = data['Electric Field/'+str.capitalize(name)].data/exunit
    if dim == '3d':
        ex = (ex[:,:,n3d//2-1]+ex[:,:,n3d//2])/2.0
    eee = 50.0
    levels = np.linspace(-eee, eee, 51)
    ex.T[ex.T < -eee]=-eee
    ex.T[ex.T >  eee]= eee
    image2=ax.contourf(X, Y, ex.T, levels=levels, cmap=cmap_br)
    #### manifesting colorbar, changing label and axis properties ####
#    cbar=plt.colorbar(image2,pad=0.01,ticks=np.linspace(-eee, eee, 5),orientation="vertical")
#    cbar.set_label(r'$E_y\ [m_ec\omega_0/e]$',fontdict=font2)        
#    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=font2['size'])        
#    plt.title('Ey_den_e at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
    #ax.text(21.,1.75,'t = 70 fs',fontdict=font)
    ax.text(42.,4.5,'t = '+str(round(time/1.0e-15,0))+' fs',fontdict=font,color='white')
    ax.set_ylim(-6.5,6.5)
    ax.set_xlim(15,55)
#    ax.set_xlabel('X [$\lambda$]',fontdict=font)
#    ax.set_ylabel('Y [$\mu m$]',fontdict=font)
    ax.tick_params(axis='y',labelsize=0.0) 
    ax.tick_params(axis='x',labelsize=0.0) 
#    ax.tick_params(axis='both',labelsize=25) 
    ax.set_xticks(ticks=[20,30,40,50])
    ax.set_yticks(ticks=[-5,0,5])


    #ax=plt.subplot(3,2,4)
    ax=plt.subplot2grid((3, 3), (1, 1),colspan=2)
    name = 'Bz_averaged'
    data = sdf.read(from_path+'b_fields'+str(n2).zfill(4)+".sdf",dict=True)
    ex = data['Magnetic Field/'+str.capitalize(name)].data/bxunit
    if dim == '3d':
        ex = (ex[:,:,n3d//2-1]+ex[:,:,n3d//2])/2.0
    eee = 10.0
    levels = np.linspace(-eee, eee, 51)
    ex.T[ex.T < -eee]=-eee
    ex.T[ex.T >  eee]= eee
    image2=ax.contourf(X, Y, ex.T, levels=levels, cmap=cmap_bg)#'BrBG')
    #### manifesting colorbar, changing label and axis properties ####
#    cbar=plt.colorbar(image2,pad=0.01,ticks=np.linspace(-eee, eee, 5),orientation="vertical")
#    cbar.set_label(r'$B_z\ [m_ec\omega_0/e]$',fontdict=font2)        
#    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=font2['size'])        
#    plt.title('Bz_averaged at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
    #ax.text(21.,1.75,'t = '+str(round(time/1.0e-15,0))+' fs',fontdict=font)
    #ax.text(21.,1.75,'t = 70 fs',fontdict=font)
#    ax.contour(X, Y, den_a.T, levels=[0.9], colors='limegreen', linewidths=2,origin='lower')
    image1=ax.contour(X, Y, den_a.T, levels=[19.5], colors='k', linewidths=2,origin='lower')
    
    ax.set_ylim(-6.5,6.5)
    ax.set_xlim(15,55)
#    ax.set_xlabel('X [$\mu m$]',fontdict=font)
#    ax.set_ylabel('Y [$\mu m$]',fontdict=font)
    ax.tick_params(axis='y',labelsize=0) 
    ax.tick_params(axis='x',labelsize=0) 
    ax.set_xticks(ticks=[20,30,40,50])
    ax.set_yticks(ticks=[-5,0,5])


    #ax=plt.subplot(3,2,6)
    ax=plt.subplot2grid((3, 3), (2, 1),colspan=2)
    name = 'Ey_averaged'
    data = sdf.read(from_path+'e_fields'+str(n2).zfill(4)+".sdf",dict=True)
    ex = data['Electric Field/'+str.capitalize(name)].data/exunit
    if dim == '3d':
        ex = (ex[:,:,n3d//2-1]+ex[:,:,n3d//2])/2.0
    eee = 10.0
    levels = np.linspace(-eee, eee, 51)
    ex.T[ex.T < -eee]=-eee
    ex.T[ex.T >  eee]= eee
    image2=ax.contourf(X, Y, ex.T, levels=levels, cmap=cmap_bg)#'BrBG')
    #### manifesting colorbar, changing label and axis properties ####
#    cbar=plt.colorbar(image2,pad=0.01,ticks=np.linspace(-eee, eee, 5),orientation="vertical")
#    cbar.set_label(r'$E_y\ [m_ec\omega_0/e]$',fontdict=font2)        
#    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=font2['size'])        
#    plt.title('Ey_averaged at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
    #ax.text(21.,1.75,'t = '+str(round(time/1.0e-15,0))+' fs',fontdict=font)
    #ax.text(21.,1.75,'t = 70 fs',fontdict=font)
#    ax.contour(X, Y, den_a.T, levels=[0.9], colors='limegreen', linewidths=2,origin='lower')
    image1=ax.contour(X, Y, den_a.T, levels=[19.5], colors='k', linewidths=2,origin='lower')
    
    ax.set_ylim(-6.5,6.5)
    ax.set_xlim(15,55)
    ax.set_xlabel('X [$\mu m$]',fontdict=font)
#    ax.set_ylabel('Y [$\mu m$]',fontdict=font)
    ax.tick_params(axis='y',labelsize=0) 
    ax.tick_params(axis='x',labelsize=25) 
    ax.set_xticks(ticks=[20,30,40,50])
    ax.set_yticks(ticks=[-5,0,5])






    plt.subplots_adjust(left=0.08, bottom=0.14, right=0.99, top=0.95, wspace=0.055, hspace=0.025)
    #ax.set_xticklabels(xticklabels,fontdict=font)
    #ax.set_yticklabels(yticklabels,fontdict=font)

#    plt.subplot(gs[0])
#    plt.scatter(ppp_x[abs(ppp_y)<=3.2],ppp_px[abs(ppp_y)<=3.2],s=10,c='green',edgecolors='None',alpha=0.5)
#    plt.xlim(0,30)
#    plt.ylabel('p$_x$ [m$_e$c]', fontdict=font)
#    plt.xticks([])
#    plt.yticks(fontsize=20)
#  
#    plt.subplot(gs[3])
#    plt.scatter(ppp_py[abs(ppp_y)<=3.2],ppp_y[abs(ppp_y)<=3.2],s=10,c='green',edgecolors='None',alpha=0.5)
#    plt.ylim(-6.5,6.5)
#    plt.xlabel('p$_y$ [m$_e$c]', fontdict=font)
#    plt.yticks([])
#    plt.xticks(fontsize=20)
#
#    
#    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.011, hspace=0.051)
#

    #fig=plt.subplot(gs[1])
    #ax1 = fig.add_axes([0.05, 0.85, 0.9, 0.10])
    #ax2 = fig.add_axes([0.05, 0.35, 0.9, 0.10])
    #cmap = mpl.cm.rainbow
    #norm = mpl.colors.Normalize(vmin=0.0, vmax=50)
    #cb1 = mpl.colorbar.ColorbarBase(ax1, cmap='Greys',
    #                            norm=norm,
    #                            orientation='horizontal',ticks=np.linspace(0.00, 50, 6))
    #cb1.set_label('n$_e$ [n$_c$]')
    #cmap = mpl.colors.ListedColormap(['r', 'g', 'b', 'c'])
    #cmap.set_over('0.25')
    #cmap.set_under('0.75')

    #cmap = mpl.cm.BrBG
    #Bz = 22.5
    #norm = mpl.colors.Normalize(vmin=-abs(Bz), vmax=abs(Bz))
    #cb2 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap_br,
    #                            norm=norm,
    #                            orientation='horizontal',ticks=np.linspace(-abs(Bz), abs(Bz), 5),alpha=0.7)
    #cb2.set_label(r'E$_y$ [m$_e\omega$/e]')
    #cmap = mpl.colors.ListedColormap(['r', 'g', 'b', 'c'])
    #cmap.set_over('0.25')
    #cmap.set_under('0.75')
    fig = plt.gcf()
    fig.set_size_inches(14, 14.)
    fig.savefig(to_path+'wrap_field.png',format='png',dpi=160)
    plt.close("all")

