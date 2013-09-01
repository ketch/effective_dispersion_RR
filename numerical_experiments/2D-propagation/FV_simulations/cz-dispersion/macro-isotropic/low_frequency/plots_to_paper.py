from clawpack.petclaw.solution import Solution
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
import matplotlib.cm as cm
from matplotlib import rc
import numpy as np
import os

def plot_p(frame,slices_xlimits=None,init_cond=False):
    sol=Solution(frame,file_format='petsc',read_aux=False,path='./_output/_p/',file_prefix='claw_p')
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers
    mx=len(x); my=len(y)
    
    mp=sol.state.num_eqn    
    yy,xx = np.meshgrid(y,x)

    if frame < 10:
        str_frame = "00"+str(frame)
    elif frame < 100:
        str_frame = "0"+str(frame)
    else:
        str_frame = str(frame)

    p=sol.state.q[0,:,:]
    pl.pcolormesh(xx,yy,p,cmap=cm.OrRd)
    pl.title("t= "+str(sol.state.t),fontsize=20)
    pl.xlabel('x',fontsize=20); pl.ylabel('y',fontsize=20)
    pl.xticks(size=20); pl.yticks(size=20)
    cb = pl.colorbar();
    #pl.clim(colorbar_min,colorbar_max);
    imaxes = pl.gca(); pl.axes(cb.ax)
    pl.yticks(fontsize=20); pl.axes(imaxes)
    pl.axis('equal')
    pl.axis([np.min(x),np.max(x),np.min(y),np.max(y)])
    pl.savefig('./_plots_to_paper/cz-dispersion_macro-iso_low-freq_pcolor_frame'+str_frame+'.png')
    pl.close()

    # write slices
    fx=open('cz-dispersion_macro-iso_low-freq_layered_xslice.txt','w')
    fx.writelines(str(xc)+" "+str(p[i,0])+"\n" for i,xc in enumerate(x))
    fx.close()

    fy=open('cz-dispersion_macro-iso_low-freq_layered_yslice.txt','w')
    fy.writelines(str(yc)+" "+str(p[0,i])+"\n" for i,yc in enumerate(y))
    fy.close()

    #pl.figure()
    #pl.plot(x,p[:,0],'b')
    #pl.plot(y,p[0,:],'r')
    #pl.savefig('./_plots_to_paper/cz-dispersion_macro-iso_low-freq_frame'+str_frame+'_slices.png')
    #pl.close()

if __name__== "__main__":
    if not os.path.exists('./_plots_to_paper'): os.mkdir('_plots_to_paper')
    
    print('**********************')
    print('**********************')
    print('Plotting p ...')

    plot_p(65)

    #frames=[0,10,20,30,40,50,60,70,80,90,100]
    #for i in xrange(0,101):
    #for i in frames:
    #    plot_p(frame=i)
    #    print i,
