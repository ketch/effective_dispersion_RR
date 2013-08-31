from clawpack.petclaw.solution import Solution
#from petclaw.io.petsc import read_petsc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
from matplotlib import rc
#rc('text', usetex=True)
import numpy as np
import os

def write_slices(frame,file_prefix,path,name):
    sol=Solution(frame,file_format='petsc',path=path,read_aux=False,file_prefix=file_prefix)
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers; my=len(y)
    q=sol.state.q
    
    f=open(name+'.txt','w')
    #f.writelines(str(xc)+" "+str(q[0,i,my/4])+" "+str(q[0,i,3*my/4])+"\n" for i,xc in enumerate(x))
    f.writelines(str(xc)+" "+str(sum(q[0,i,:])/my)+"\n" for i,xc in enumerate(x))
    f.close()

if __name__== "__main__":

    write_slices(970,'claw_p','./_output/_p/','stress')
