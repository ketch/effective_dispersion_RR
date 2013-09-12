from clawpack import petclaw
from psystem import setaux
import matplotlib.pyplot as plt
import numpy as np

sol = petclaw.Solution(100,file_format='petsc')
xc = sol.state.grid.x.centers
yc = sol.state.grid.y.centers
q = sol.q
aux = setaux(xc,yc)
u = q[1,...]/aux[0,...]
v = q[2,...]/aux[0,...]

X, Y = np.meshgrid(xc,yc,indexing='ij')

plt.figure()
plt.pcolor(X,Y,q[0,...],cmap='RdBu')
plt.hold(True)
x_skip = 8
y_skip = 4
plt.quiver(X[::x_skip,::y_skip],Y[::x_skip,::y_skip],u[::x_skip,::y_skip],10*v[::x_skip,::y_skip],pivot='center',units='width',linewidth=1.5,headaxislength=5,scale=100)
plt.axis([18,32,0,1])
plt.xlabel('x'); plt.ylabel('y')
plt.hold(False)
plt.show()
