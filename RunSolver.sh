#!/bin/sh
# This script is used to run the compiled ./Solve for the different parameters as required in
# the report. The data is then saved and read via a matlab file from which the plots are generated.

make

mpiexec -np 9 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=3 --Py=3 --dt=0.0005 --T=1.0 --Re=100
mpiexec -np 9 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=3 --Py=3 --dt=0.0005 --T=1.0 --Re=400
mpiexec -np 9 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=3 --Py=3 --dt=0.0001 --T=5.0 --Re=1000
mpiexec -np 9 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=3 --Py=3 --dt=0.0001 --T=10.0 --Re=3200
mpiexec -np 9 ./Solve --Lx=1 --Ly=2 --Nx=161 --Ny=161 --Px=3 --Py=3 --dt=0.0005 --T=1.0 --Re=100
mpiexec -np 9 ./Solve --Lx=2 --Ly=1 --Nx=161 --Ny=161 --Px=3 --Py=3 --dt=0.0005 --T=1.0 --Re=100

"/apps/matlab/2018b/bin/matlab" -nodisplay -nosplash -nodesktop -r "run('/home/jsc4017/Desktop/hpc_coursework/GeneratePlots.m');exit;"
