import sys, petsc4py
petsc4py.init(sys.argv)

import numpy as np
from petsc4py import PETSc
from fdtd_da_pml import em_da_q3  as em_q3
from fdtd_da_pml import em_da_q12 as em_q12

nx = 32
ny = 32
num_aux = 3
num_pml = 8

stype  = PETSc.DA.StencilType.STAR
swidth = 1
da = PETSc.DA().create([nx,ny], dof=1,
                       stencil_type=stype,
                       stencil_width=swidth)

( xi, xf), ( yi, yf) = da.getRanges()
(gxi,gxf), (gyi,gyf) = da.getGhostRanges()

Q1 = da.createGlobalVec()
Q2 = da.createGlobalVec()
Q3 = da.createGlobalVec()
Q1loc = da.createLocalVec()
Q2loc = da.createLocalVec()
Q3loc = da.createLocalVec()

q1 = Q1loc.getArray().reshape([gxf-gxi,gyf-gyi], order='F')
q2 = Q2loc.getArray().reshape([gxf-gxi,gyf-gyi], order='F')
q3 = Q3loc.getArray().reshape([gxf-gxi,gyf-gyi], order='F')

s1  = np.zeros([gxf-gxi,gyf-gyi], order='F')
s2  = np.zeros([gxf-gxi,gyf-gyi], order='F')
s3  = np.zeros([gxf-gxi,gyf-gyi], order='F')
s4  = np.zeros([gxf-gxi,gyf-gyi], order='F')
aux = np.ones ([num_aux,gxf-gxi,gyf-gyi], order='F')
pml = np.ones ([num_pml,gxf-gxi,gyf-gyi], order='F')

xi  += 1
yi  += 1
gxi += 1
gyi += 1

def bc(Q1,Q2,Q3):
    q1 = da.getVecArray(Q1)
    q2 = da.getVecArray(Q2)
    q3 = da.getVecArray(Q3)
    if xi == 1: 
        q2[ 0,:] = 1.0
    if yi == 1:
        q1[:, 0] = 0.0
    if yf == ny:
        #q1[:,-1] = 0.0
        #q2[:,-1] = 0.0
        #q3[:,-1] = 0.0
        pass
    if xf == nx:
        #q1[-1,:] = 0.0
        #q2[-1,:] = 0.0
        #q3[-1,:] = 0.0
        pass

io   = PETSc.Viewer.BINARY()
draw = PETSc.Viewer.DRAW()
root = da.comm.getRank() == 0
for t in range(1,100):
    #if root: print "timestep:", t

    if t == 1: bc(Q1,Q2,Q3)

    da.globalToLocal(Q1, Q1loc)
    da.globalToLocal(Q2, Q2loc)
    em_q3(aux,pml,s3,s4,q1,q2,q3,xi,xf,yi,yf,gxi,gxf,gyi,gyf)
    da.localToGlobal(Q3loc, Q3)

    da.globalToLocal(Q3, Q3loc)
    em_q12(aux,pml,s1,s2,q1,q2,q3,xi,xf,yi,yf,gxi,gxf,gyi,gyf)
    da.localToGlobal(Q1loc, Q1)
    da.localToGlobal(Q2loc, Q2)

    #bc(Q1,Q2,Q3)

    #Q3.view(draw)
    #Q3.view(io)
    
sys.exit()
