import matplotlib.pyplot as plt
from clawpack.petclaw.solution import Solution
import numpy as np

skip = 4

#c-dispersive
frame = 90
sol = Solution(frame,file_format='petsc',read_aux=False,path='./c-disp/',file_prefix='claw_p')
p=sol.state.q[0,:,:]

np.save('pc'+str(skip)+'.npy',p[::skip,::skip])

x=sol.state.grid.x.centers; y=sol.state.grid.y.centers
yy,xx = np.meshgrid(y,x)
np.save('x'+str(skip)+'_1.npy',xx[::skip,::skip])
np.save('y'+str(skip)+'_1.npy',yy[::skip,::skip])


#cZ-dispersive
frame = 90
sol = Solution(frame,file_format='petsc',read_aux=False,path='./cz-disp/',file_prefix='claw_p')
p=sol.state.q[0,:,:]

np.save('pcz'+str(skip)+'.npy',p[::skip,::skip])


#Z-dispersive
frame = 90
sol = Solution(frame,file_format='petsc',read_aux=False,path='./z-disp/',file_prefix='claw_p')
p=sol.state.q[0,:,:]

np.save('pz'+str(skip)+'.npy',p[::skip,::skip])

#"Homogeneous"
frame = 90
sol = Solution(frame,file_format='petsc',read_aux=False,path='./homog/',file_prefix='claw_p')
p=sol.state.q[0,:,:]

np.save('p'+str(skip)+'.npy',p[::skip,::skip])

# Lower resolution for homogeneous simulation
x=sol.state.grid.x.centers; y=sol.state.grid.y.centers
yy,xx = np.meshgrid(y,x)
np.save('x'+str(skip)+'_2.npy',xx[::skip,::skip])
np.save('y'+str(skip)+'_2.npy',yy[::skip,::skip])

