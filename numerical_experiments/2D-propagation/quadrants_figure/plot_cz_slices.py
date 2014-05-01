import matplotlib.pyplot as plt
import numpy as np
import scipy.io

skip = 4

#Load FV data
pcz = np.load('pcz'+str(skip)+'.npy')
xx1 = np.load('x'+str(skip)+'_1.npy')
yy1 = np.load('y'+str(skip)+'_1.npy')
FV_xslice = pcz[:,0]

# Load pseudospectral data
D = scipy.io.loadmat('cz-disp.mat')
phom = D['U']
mx = D['mx']
x  = np.squeeze(D['x'])[mx/2:]-100
xslice = np.squeeze(phom[mx/2,mx/2:])

my = D['my']
y  = np.squeeze(D['y'])[my/2:]-100
yslice = np.squeeze(phom[my/2:,my/2])

plt.plot(xx1[:,0], pcz[:,0],'b',linewidth=2)
plt.plot(yy1[0,:], pcz[0,:],'r',linewidth=2)
plt.plot(x,xslice,'--k',linewidth=2)
plt.plot(y,yslice,'--k',linewidth=2)

plt.legend(('Variable-coefficient y=0','Variable-coefficient x=0','Homgenized y=0','Homogenized x=0'),loc='best')
