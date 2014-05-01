import numpy as np
import matplotlib.pyplot as plt

p = np.load('pc4.npy')
q = np.load('qc4.npy')
x = np.load('x4_1.npy')
y = np.load('y4_1.npy')

a = 1200
c = 1520
xskip = 2
b = 17
xx = np.squeeze(x[a:c+1,0])
yy = np.squeeze(y[0,:b+1])

plt.pcolormesh(x[a:c+1,:b+1],y[a:c+1,:b+1],p[a:c,:b],cmap='RdBu_r')
plt.clim(-0.025,0.025)
plt.streamplot(xx,yy,np.squeeze(q[1,a:c+1,:b+1]).T,np.squeeze(q[2,a:c+1,:b+1]).T,density=2,linewidth=2,color='k')
plt.ylim(0,1)
plt.show()
