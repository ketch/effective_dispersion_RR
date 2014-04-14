#!/usr/bin/env bash
#
# @ job_name            = test_pyclaw
# @ job_type            = bluegene
# @ output              = ./$(job_name)_$(jobid).out
# @ error               = ./$(job_name)_$(jobid).err
# @ environment         = COPY_ALL; 
# @ wall_clock_limit    = 4:00:00,4:00:00
# @ notification        = always
# @ bg_size             = 256
# @ account_no          = k47

# @ queue


projdir=/home/ketch/test_build_64
builddir=${projdir}/opt/share
pythondir=${builddir}/python/2.7.2/bgp
ldpath=${pythondir}/lib
numpy_path=${builddir}/numpy/1.6.2/bgp/lib/python
nose_path=${builddir}/nose/1.1.2/bgp/lib/python
clawpack_path=${builddir}/clawpack/dev/bgp/lib/python
petsc4py_path=${builddir}/petsc4py/1.2/bgp/lib/python
bgp_python_path=${numpy_path}:${nose_path}:${clawpack_path}:${petsc4py_path}
sandbox=${projdir}/sandbox
ldpath=${pythondir}/lib
logdir=${builddir}/logs

mpirun -env LD_LIBRARY_PATH=${ldpath} -env PYTHONPATH=${bgp_python_path} \
    -mode VN -exp_env HOME -n 1024 ${pythondir}/bin/python \
    psystem_lin.py Nx=64 Ny=64 \
    | tee ${logdir}/psystem.log
