multiplier = 1
textcolor = '#000000'
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica Neue'],'size': 15*multiplier})
rc('text', usetex=True)

from matplotlib import font_manager
ticks_font = font_manager.FontProperties(family='Bitstream Vera Sans')

import matplotlib.pyplot as plt
import numpy as np


nullfmt = plt.NullFormatter()
nullloc = plt.NullLocator()

tick_vals = np.array( (0,25,50,75,100),dtype=int )

skip = 4

p   = np.load('p'+str(skip)+'.npy')
pc  = np.load('pc'+str(skip)+'.npy')
pz  = np.load('pz'+str(skip)+'.npy')
pcz = np.load('pcz'+str(skip)+'.npy')
xx1 = np.load('x'+str(skip)+'_1.npy')
yy1 = np.load('y'+str(skip)+'_1.npy')
xx2 = np.load('x'+str(skip)+'_2.npy')
yy2 = np.load('y'+str(skip)+'_2.npy')

yfrac = yy2-np.floor(yy2)
mat = yfrac>0.5

fig = plt.figure(1,figsize=(10*multiplier,10*multiplier))

d2 = 5./16 # side length of main figures
d1 = d2/5. # side length of trace figures
d0 = d1*2. # border

mainframe = [0]*4
xframe    = [0]*4
yframe    = [0]*4

# Indexing: 0 = Top left
#           1 = Top right
#           2 = Bottom left
#           3 = Bottom right

mainframe[0] = [d0,d0+2*d1+d2   ,d2,d2]
xframe[0] =    [d0,d0+d1+d2     ,d2,d1]
yframe[0] =    [d0+d2,d0+2*d1+d2   ,d1,d2]

mainframe[1] = [d0+2*d1+d2   ,d0+2*d1+d2,d2,d2]
xframe[1]    = [d0+2*d1+d2   ,d0+d1+d2,d2,d1]
yframe[1]    = [d0+d1+d2,d0+2*d1+d2,d1,d2]

mainframe[2] = [d0,d0,d2,d2]
xframe[2] =    [d0,d0+d2   ,d2,d1]
yframe[2] =    [d0+d2,d0,d1,d2]

mainframe[3] = [d0+2*d1+d2,d0,d2,d2]
xframe[3] =    [d0+2*d1+d2,d0+d2,d2,d1]
yframe[3] =    [d0+d1+d2,d0,d1,d2]

on_top   = [ 1,1,-1,-1]
on_right = [-1,1,-1, 1]

axmain = []
axx = []
axy = []
axd = []

kx = np.linspace(0.01,3,500)
ky = np.linspace(0.01,3,500)
kkx,kky = np.meshgrid(kx,ky)

for i,p in enumerate([p,pc,pz,pcz]):
    clim = [-0.25,0.25]
    if i in (1,2,3): # new runs
        xx = xx1
        yy = yy1
    else: # old runs
        xx = xx2
        yy = yy2
    # Main plot
    axmain.append(plt.axes(mainframe[i]))
    #if i>0:
        #plt.imshow(mat.T,extent=[xx.min(),xx.max(),yy.min(),yy.max()],interpolation='nearest',origin='lower',cmap='bone',alpha=0.5)
    p_scaled = (np.abs(p)**(1./3))*np.sign(p)
    #plt.pcolormesh(on_right[i]*xx,on_top[i]*yy,p_scaled,cmap='RdBu_r')
    p_scaled = p_scaled.T
    if on_right[i] == -1:
        p_scaled = np.fliplr(p_scaled)
        xmin = -xx.max(); xmax = -xx.min()
    else:
        xmin = xx.min(); xmax = xx.max()
    if on_top[i] == -1:
        p_scaled = np.flipud(p_scaled)
        ymin = -yy.max(); ymax = -yy.min()
    else:
        ymin = yy.min(); ymax = yy.max()
    plt.imshow(p_scaled, cmap='RdBu_r',extent=[xmin,xmax,ymin,ymax],interpolation='nearest',origin='lower',alpha=1.0)
    plt.clim(clim)
    plt.axis('equal'); plt.axis('tight')

    if on_top[i] == 1:
        axmain[-1].xaxis.set_major_locator(plt.NullLocator())
    else:
        axmain[-1].set_xticks(on_right[i]*tick_vals)
        axmain[-1].set_xticklabels(axmain[-1].get_xticks(),ticks_font,color=textcolor)

    if on_right[i] == -1:  # Tick labels on left
        axmain[-1].set_yticks(on_top[i]*tick_vals)
        axmain[-1].set_yticklabels(axmain[-1].get_yticks(),ticks_font,color=textcolor)
    else:
        axmain[-1].yaxis.set_major_locator(plt.NullLocator())

    # x-slice
    axx.append(plt.axes(xframe[i], frameon=False))
    axx[-1].grid(color='w', linewidth=2*multiplier, linestyle='solid',alpha=0.4)
    axx[-1].plot(on_right[i]*xx[:,0],p[:,0],'-r',lw=2*multiplier,alpha=0.8)
    axx[-1].xaxis.set_major_locator(nullloc)
    axx[-1].yaxis.set_major_locator(nullloc)
    if on_right[i] == -1:
        axx[-1].set_xlim((-100,0)); 
    else:
        axx[-1].set_xlim((0,100)); 

    # Homogenized solution x-slice
    if i == 1:
        import scipy.io
        D = scipy.io.loadmat('c-disp.mat')
    elif i == 2:
        D = scipy.io.loadmat('z-disp.mat')
    elif i == 3:
        D = scipy.io.loadmat('cz-disp.mat')
    if i>0:
        phom = D['U']
        if i==1:
            phom = phom/5.  # Accidentally ran this one with wrong initial amplitude
        mx = D['mx']
        x  = np.squeeze(D['x'])[mx/2:]-100
        xslice = np.squeeze(phom[mx/2,mx/2:])
        axx[i].plot(on_right[i]*x,xslice,'--k',lw=2*multiplier,alpha=1.0)

    # y-slice
    axy.append(plt.axes(yframe[i],frameon=False))
    axy[-1].plot(on_right[i]*p[0,:],on_top[i]*yy[0,:],'-r',lw=2*multiplier,alpha=0.8)
    axy[-1].xaxis.set_major_locator(nullloc)
    axy[-1].yaxis.set_major_locator(nullloc)
    if on_top[i] == -1:
        axy[-1].set_ylim((-100,0)); 
    else:
        axy[-1].set_ylim((0,100)); 

    # Homogenized solution x-slice
    if i > 0:
        my = D['my']
        y  = np.squeeze(D['y'])[my/2:]-100
        yslice = np.squeeze(phom[my/2:,my/2])
        axy[i].plot(on_right[i]*yslice,on_top[i]*y,'--k',lw=2*multiplier)



# Top left
axmain[0].text(0.02,0.98,'Homogeneous\n medium',transform=axmain[0].transAxes,verticalalignment='top')
axx[0].set_ylim((-0.05,0.05))
axy[0].set_xlim((-0.05,0.05))

# Top right
axmain[1].text(0.98,0.98,'Sound speed\n mismatched',horizontalalignment='right',transform=axmain[1].transAxes,verticalalignment='top')
axx[1].set_ylim((-0.05,0.05))
axy[1].set_xlim((-0.05,0.05))

#Bottom left
axmain[2].text(0.02,0.02,'Impedance\n mismatched',transform=axmain[2].transAxes)
axx[2].set_ylim((-0.07,0.07))
axy[2].set_xlim((-0.04,0.04))
axx[2].xaxis.set_label_coords(1, -0.5)

# Bottom right
axmain[3].text(0.98,0.02,'Both\n mismatched',horizontalalignment='right',transform=axmain[3].transAxes)
axx[3].set_ylim((-0.05,0.05))
axy[3].set_xlim((-0.04,0.04))


#Inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

axins = zoomed_inset_axes(axmain[2],6,loc=1)

yfrac = yy-np.floor(yy)
mat = yfrac>0.5

axins.imshow(mat[:80/skip,:320/skip].T,vmin=-0.5, vmax=1.5, extent=[-40,-35,-40,-35],interpolation='nearest',origin='lower',cmap='binary',alpha=0.8)
axins.xaxis.set_major_locator(nullloc)
axins.yaxis.set_major_locator(nullloc)
mark_inset(axmain[2],axins,2,4,fc='none',ec='0.5')

#plt.suptitle('Acoustic wave propagation in layered periodic media',fontsize=25*multiplier,color=textcolor)
#plt.figtext(0.5,0.05,'Units scaled to medium period', fontsize=15*multiplier, horizontalalignment='center',color=textcolor)
plt.figtext(0.5,0.075,'$x$', fontsize=25*multiplier, horizontalalignment='center',color=textcolor)
plt.figtext(0.05,0.5,'$y$', fontsize=25*multiplier, verticalalignment='center',color=textcolor)

import matplotlib
xl,yl = np.array([[-20., -20.], [-140., 100.]])
vert = matplotlib.lines.Line2D(xl, yl, lw=1., color='k', alpha=0.8)
axmain[1].add_line(vert)
vert.set_clip_on(False)
xl,yl = np.array([[-140., 100.], [-20., -20.]])
horiz = matplotlib.lines.Line2D(xl, yl, lw=1., color='k', alpha=0.8)
axmain[1].add_line(horiz)
horiz.set_clip_on(False)

plt.savefig('_dispersion_io.pdf',dpi=500)
plt.close()
