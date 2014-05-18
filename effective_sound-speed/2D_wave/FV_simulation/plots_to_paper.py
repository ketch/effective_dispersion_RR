from clawpack.petclaw.solution import Solution
#from petclaw.io.petsc import read_petsc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
import matplotlib.cm as cm
from matplotlib import rc
#rc('text', usetex=True)
import numpy as np
import os
import scipy.io

def plot_p(frame):
    sol=Solution(frame,file_format='petsc',read_aux=False,path='./_output/_p/',file_prefix='claw_p')
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers
    mx=len(x); my=len(y)
    
    mp=sol.state.num_eqn    
    yy,xx = np.meshgrid(y,x)

    p=sol.state.q[0,:,:]
    fig = pl.figure(figsize=(8, 3.5))
    #pl.title("t= "+str(sol.state.t),fontsize=20)
    pl.xticks(size=20); pl.yticks(size=20)
    pl.xlabel('x',fontsize=20); pl.ylabel('y',fontsize=20)
    #pl.pcolormesh(xx,yy,p_subxy,cmap=cm.OrRd)
    pl.pcolormesh(xx,yy,p,cmap='RdBu_r')
    pl.autoscale(tight=True)
    cb = pl.colorbar(ticks=[0.5,1,1.5,2]);
    
    #pl.clim(ticks=[0.5,1,1.5,2])
    imaxes = pl.gca(); pl.axes(cb.ax)
    pl.yticks(fontsize=20); pl.axes(imaxes)
    #pl.xticks(fontsize=20); pl.axes(imaxes)
    #pl.axis('equal')
    pl.axis('tight')
    fig.tight_layout()
    pl.savefig('./_plots_to_paper/sound-speed_FV_t'+str(frame)+'_pcolor.png')
    pl.close()

def plot_p_leading_order(frame):
    mat = scipy.io.loadmat('sound-speed_2D-wave.mat')
    T=5; nt=T/0.5
    pp=mat['U'][nt,:,:]
    xx=mat['xx']
    yy=mat['yy']

    fig=pl.figure(figsize=(8, 3.5))
    #pl.title("t= "+str(sol.state.t),fontsize=20)
    pl.xticks(size=20); pl.yticks(size=20)
    pl.xlabel('x',fontsize=20); pl.ylabel('y',fontsize=20)
    #pl.pcolormesh(xx,yy,p_subxy,cmap=cm.OrRd)
    pl.pcolormesh(xx,yy,pp,cmap='RdBu_r')
    pl.autoscale(tight=True)
    cb = pl.colorbar(ticks=[0.5,1,1.5,2]);
    
    #pl.clim(ticks=[0.5,1,1.5,2])
    imaxes = pl.gca(); pl.axes(cb.ax)
    pl.yticks(fontsize=20); pl.axes(imaxes)
    #pl.xticks(fontsize=20); pl.axes(imaxes)
    #pl.axis('equal')
    pl.axis('tight')
    fig.tight_layout()
    pl.savefig('./_plots_to_paper/sound-speed_LO_t'+str(frame)+'_pcolor.png')
    pl.close()

if __name__== "__main__":
    if not os.path.exists('./_plots_to_paper'): os.mkdir('_plots_to_paper')

    print('**********************')
    print('**********************')
    print('Plotting p ...')
    plot_p(frame=5)
    plot_p_leading_order(frame=5)
