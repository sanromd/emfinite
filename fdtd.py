#!/usr/bin/env python
# encoding: utf-8
# ddd
# sanity check for submodule working and pushing to...
import sys, petsc4py
petsc4py.init(sys.argv)

import numpy as np

# space
xx = 1e-6
yy = 500e-6

ddx = 1#5e-9
ddy = 1#5e-9

nx = 35
ny = 35

# -------- GLOBAL SCALAR DEFINITIONS -----------------------------
# ======== all definitions are in m,s,g unit system.
n_frames = 30
# ....... dimensions .............................................
x_lower = 0.0e-6
x_upper = 100e-6                    # lenght [m]
y_lower = 0.0e-6
y_upper = 1.0e-6                    # notice that for multilayer this is value will be over-written
# ........ material properties ...................................

# vacuum
eo = 8.854187817e-12            # vacuum permittivity   - [F/m]
mo = 4e-7*np.pi                 # vacuum peremeability  - [V.s/A.m]
co = 1/np.sqrt(eo*mo)           # vacuum speed of light - [m/s]
zo = np.sqrt(eo/mo)
# material
mat_shape = 'homogeneous'           # material definition: homogeneous, interface, rip (moving perturbation), multilayered
num_aux = 3
# background refractive index
bkg_er = 1.5
bkg_mr = 1.5
bkg_n  = np.sqrt(bkg_er*bkg_mr)
bkg_e  = eo*bkg_er
bkg_m  = mo*bkg_mr

# if interface declare position
x_change = x_upper/2

# set moving refractive index parameters
rip_vx_e    = 0.0*co    # replace here the value of x
rip_vx_m    = rip_vx_e
rip_vy_e    = 0.0*co
rip_vy_m    = rip_vy_e

rip_xoff_e  = 10e-6
rip_xoff_m  = rip_xoff_e
rip_yoff_e  = rip_xoff_e
rip_yoff_m  = rip_xoff_e

rip_xsig_e  = 10.0e-6
rip_xsig_m  = rip_xsig_e
rip_ysig_e  = .9*y_upper/2
rip_ysig_m  = rip_ysig_e
s_x_e       = rip_xsig_e**2
s_x_m       = rip_xsig_m**2
s_y_e       = rip_ysig_e**2
s_y_m       = rip_ysig_m**2

prip        = 0.1
deltan      = prip*(bkg_n) # assumes epsilon = mu
d_e         = deltan #*(2.0*1.5+deltan)
d_m         = deltan #*(2.0*1.5+deltan)

# set multilayer parameters

# multilayered definition
n_layers = 2
layers = np.zeros([n_layers,7]) # _layer:  eps mu N t chi2e chi2m chi3e chi3m
layers[0,0] = 1.5
layers[0,1] = 1.5
layers[0,2] = 10
layers[0,3] = 15e-9
layers[1,0] = 2.5
layers[1,1] = 2.5
layers[1,2] = layers[0,2] - 1
layers[1,3] = 50e-9
N_layers = 5
if mat_shape=='multilayer':
    y_upper = N_layers*np.sum(layers[:,3])+layers[0,3]
    tlp = np.sum(layers[:,3])
    mlp = np.floor(tlp/1e-9)

# ........ Boundary settings settings .................

bc_x_lower = 'scattering'
bc_x_upper = 'none'
bc_y_lower = 'scattering'
bc_y_upper = 'none'

aux_bc_x_lower = 'scattering'
aux_bc_x_upper = 'pml'
aux_bc_y_upper = 'metallic'
aux_bc_y_upper = 'metallic'

# parameters needed for pml calculation
num_pml = 8
norder = 3
Ro = 1.0e-6

# ........ excitation - initial conditoons .......................

ex_type  = 'plane'
alambda  = 1e-6             # wavelength
ex_t_sig = 1.0*alambda          # width in time (pulse only)
ex_x_sig = 1.0*alambda          # width in the x-direction (pulse)
ex_y_sig = y_upper-y_lower
ex_toff  = 0.0                  # offset in time
ex_xoff  = 0.0                  # offset in the x-direction
ex_yoff  = y_upper/2            # offset in the y-direction
omega    = 2.0*np.pi*co/alambda # frequency
k        = 2.0*np.pi/alambda
amp_Ex   = 0.
amp_Ey   = 1.
amp_Hz   = 1.


# ........ pre-calculations for wave propagation .................

v_r = 1./bkg_n
v = co*v_r
ex_vx = v
ex_vy = 0.0
ex_kx = k
ex_ky = 0.0

# Grid - mesh settings
mx = np.floor(40*(x_upper-x_lower)/alambda)
if mat_shape=='multilayer':
    my = np.floor((y_upper-y_lower)/1e-9)
else:
    my = np.floor(20*(y_upper-y_lower)/alambda)

ddx = (x_upper-x_lower)/mx
ddy = (y_upper-y_lower)/my
ddt = dt = 1#0.90/(co*np.sqrt(1.0/(ddx**2)+1.0/(ddy**2)))
max_steps = 250000
t_final = (x_upper-x_lower)/v

dxdt = 1
dydt = 1

# -------- GLOBAL FUNCTION DEFINITIONS --------------

def etar(num_aux,xi,xf,yi,yf,ddx,ddy):
    """
    eta = etar(num_aux,xi,xf,yi,yf,ddx,ddy)

    Sets the auxiliary arrays for permittivity and permeability.

    Implemented mappings

    ..gaussian1dx:  stationary and moving gaussian shape for eps and mu
    ..homogeneous:  homogeneous refractive index in eps and mu
    ..interface:    simple interface (jump) acroos the 2d domain
    ..interfacex:   simple interface (jump) 1D in x-direction
    ..interfacey:   ibid in y-direction
    ..multilayer:   2D multilayers in x or y direction.


    y,x are the point coordinates of the grid.


    on output aux holds:

                               EM equivalent

         idim   curvilinear  |   TE      TM
         0:     eta_1        |   mu1     eps1
         1:     eta_2        |   mu2     eps2
         2:     eta_3        |   eps3    mu3

    """

    X = np.linspace(xi,xf,xf-xi)*ddx
    Y = np.linspace(yi,yf,yf-yi)*ddy
    y,x = np.meshgrid(Y,X)
    eta = np.empty( [3,len(X),len(Y)], order='F')

    if mat_shape=='gaussian1dx':
        u_x_e = x - rip_vx_e*t - rip_xoff_e
        u_x_m = x - rip_vx_m*t - rip_xoff_m
        u_y_e = y - rip_vy_e*t - rip_yoff_e
        u_y_m = y - rip_vy_m*t - rip_yoff_m

        u_e = (u_x_e/rip_xsig_e)**2 + (u_y_e/rip_ysig_e)**2
        u_m = (u_x_m/rip_xsig_m)**2 + (u_y_m/rip_ysig_m)**2

        eta[0,:,:] = d_e*np.exp(-u_e) + bkg_er
        eta[1,:,:] = d_e*np.exp(-u_e) + bkg_er
        eta[2,:,:] = d_m*np.exp(-u_m) + bkg_mr
    elif mat_shape=='homogeneous':
        eta[0,:,:] = bkg_er
        eta[1,:,:] = bkg_er
        eta[2,:,:] = bkg_mr
    elif mat_shape=='interfacex':
        eta[0,:,:] = 1*(x<x_change) + 4*(x>=x_change)
        eta[1,:,:] = 1*(x<x_change) + 4*(x>=x_change)
        eta[2,:,:] = 1*(x<x_change) + 4*(x>=x_change)
    elif mat_shape=='interfacey':
        yy = y_upper-y_lower
        eta[0,:,:] = 1*(y<yy/2) + 4*(x>=yy/2)
        eta[1,:,:] = 1*(y<yy/2) + 4*(x>=yy/2)
        eta[2,:,:] = 1*(y<yy/2) + 4*(x>=yy/2)
    elif mat_shape=='multilayer':
        for n in range(0,N_layers):
            yi = n*tlp
            for m in range(0,n_layers):
                if m==0:
                    eta[0,:,:] = layers[m,0]*(yi<y)*(y<=yi+layers[m,3])
                    eta[1,:,:] = layers[m,1]*(yi<y)*(y<=yi+layers[m,3])
                    eta[2,:,:] = layers[m,1]*(yi<y)*(y<=yi+layers[m,3])
                else:
                    eta[0,:,:] = layers[m,0]*(yi+layers[m-1,3]<y)*(y<=yi+layers[m,3])
                    eta[1,:,:] = layers[m,1]*(yi+layers[m-1,3]<y)*(y<=yi+layers[m,3])
                    eta[2,:,:] = layers[m,1]*(yi+layers[m-1,3]<y)*(y<=yi+layers[m,3])


        eta[0,:,:] = layers[0,0]*(N_layers*tlp<y)*(y<=N_layers*tlp+layers[0,3])
        eta[1,:,:] = layers[0,1]*(N_layers*tlp<y)*(y<=N_layers*tlp+layers[0,3])
        eta[2,:,:] = layers[0,1]*(N_layers*tlp<y)*(y<=N_layers*tlp+layers[0,3])

    return eta

def qinit(Q1,Q2,Q3,da):
    """
    Set initial conditions for q in the grid, but not on the boundary lines
    """
    nx, ny = da.getSizes()
    (xi, xf), (yi, yf) = da.getRanges()
    q1 = da.getVecArray(Q1)
    q2 = da.getVecArray(Q2)
    q3 = da.getVecArray(Q3)
    q1[:,:] = 0.0
    q2[:,:] = 0.0
    q3[:,:] = 0.0

def qbc(Q1,Q2,Q3,da):
    """
    Set the boundary conditions for q. Implemented conditions are:

    ..metallic:     set all qs to 0
    ..scattering:   scattering boundary conditions (line), set by function bc_scattering (user controlled)
    ..periodic:     periodic boundary conditions            (not implemented)
    ..neumann:      neumann boundary conditions             (not implemented)
    ..rounded:      end boundary becomes beginning boundary (not implemented)
    ..none:         pass, does not set any value
    """
    nx, ny = da.getSizes()
    (xi, xf), (yi, yf) = da.getRanges()

    q1  = da.getVecArray(Q1)
    q2  = da.getVecArray(Q2)
    q3  = da.getVecArray(Q3)

    if xi == 0:
        if bc_x_lower == 'metallic':
            q1[0,:] = 0.0
            q2[0,:] = 0.0
            q3[0,:] = 0.0
        elif bc_x_lower == 'scattering':
            bc_scattering(Q1,Q2,Q3,da,0,0)
        elif bc_x_lower == 'none':
            pass
    if xi == nx-1:
        if bc_x_upper == 'metallic':
            q1[-1,:] = 0.0
            q2[-1,:] = 0.0
            q3[-1,:] = 0.0
        elif bc_x_upper == 'scattering':
            bc_scattering(Q1,Q2,Q3,da,0,1)
        elif bc_x_upper == 'none':
            pass
    if yi == 0:
        if bc_y_lower == 'metallic':
            q1[:,0] = 0.0
            q2[:,0] = 0.0
            q3[:,0] = 0.0
        elif bc_y_lower == 'scattering':
            bc_scattering(Q1,Q2,Q3,da,1,0)
        elif bc_y_lower == 'none':
            pass
    if yf == ny-1:
        if bc_y_upper == 'metallic':
            q1[:,-1] = 0.0
            q2[:,-1] = 0.0
            q3[:,-1] = 0.0
        elif bc_y_upper == 'scattering':
            bc_scattering(Q1,Q2,Q3,da,1,1)
        elif bc_y_upper == 'none':
            pass

def bc_scattering(Q1,Q2,Q3,da,axis,side):
    """
    Boundary scattering conditions (source)
    """
    nx, ny = da.getSizes()
    (xi, xf), (yi, yf) = da.getRanges()

    q1  = da.getVecArray(Q1)
    q2  = da.getVecArray(Q2)
    q3  = da.getVecArray(Q3)

    if (axis,side) == (0,0) and xi == 0:
        q2[ 0,:] = 1.0
    if (axis,side) == (1,0) and yi == 0:
        q1[:, 0] = 0.0
    if (axis,side) == (0,1) and xf == nx-1:
        #q1[-1,:] = 0.0
        #q2[-1,:] = 0.0
        #q3[-1,:] = 0.0
        pass
    if (axis,side) == (1,1) and yf == ny-1:
        #q1[:,-1] = 0.0
        #q2[:,-1] = 0.0
        #q3[:,-1] = 0.0
        pass

def aux_bc_pml(pml,pml_type,xi,xf,yi,yf,nx,ny):
    """
    Set  PML on the auxiliary boundary conditions.
    """
    from build_pml import build_pml
    ddx,ddy,dt,co = 1,1,1,1
    build_pml(pml,pml_type,ddx,ddy,dt,norder,Ro,co,xi+1,xf,yi+1,yf,nx,ny)

def aux_bc():
    """
    user controlled boundary auxiliary conditions
    """
    pass


# create DA and allocate global and local variables
from petsc4py import PETSc
from fdtd_da_pml import fdtd_2d

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

s1  = np.zeros([xf-xi,yf-yi], order='F')
s2  = np.zeros([xf-xi,yf-yi], order='F')
s3  = np.zeros([xf-xi,yf-yi], order='F')
s4  = np.zeros([xf-xi,yf-yi], order='F')

#aux = etar(num_aux,gxi,gxf,gyi,gyf,ddx,ddy)
aux = np.ones ([num_aux,xf-xi,yf-yi], order='F')

#da_pml = PETSc.DA().create([mx,my], dof=num_pml,
#                       stencil_type=stype,
#                       stencil_width=swidth)
#PML = da_pml.createGlobalVec()
#PMLloc = da_pml.createLocalVec()
#pml = PMLloc.getArray().reshape([num_pml,gxf-gxi,gyf-gyi], order='F')
pml = np.ones ([num_pml,xf-xi,yf-yi], order='F')
pml_axis = 0
pml_side = 1
pml_type = pml_axis*2+pml_side
aux_bc_pml(pml,pml_type,xi,xf,yi,yf,nx,ny)
#from matplotlib import pylab
#for i in range(8):
#    pylab.figure()
#    pylab.contourf(pml[i,:,:].copy())
#    pylab.colorbar()
#pylab.show()

io   = PETSc.Viewer.BINARY()
draw = PETSc.Viewer.DRAW()
root = da.comm.getRank() == 0

for t in range(1,100):
    if t == 1:
        qinit(Q1,Q2,Q3,da)
        qbc(Q1,Q2,Q3,da)

    da.globalToLocal(Q1, Q1loc)
    da.globalToLocal(Q2, Q2loc)
    fdtd_2d(aux,pml,dxdt,dydt,s1,s2,s3,s4,q1,q2,q3,xi+1,xf,yi+1,yf,gxi+1,gxf,gyi+1,gyf,0,1)
    da.localToGlobal(Q3loc, Q3)

    da.globalToLocal(Q3, Q3loc)
    fdtd_2d(aux,pml,dxdt,dydt,s1,s2,s3,s4,q1,q2,q3,xi+1,xf,yi+1,yf,gxi+1,gxf,gyi+1,gyf,1,1)
    da.localToGlobal(Q1loc, Q1)
    da.localToGlobal(Q2loc, Q2)

    #bc(Q1,Q2,Q3)
    Q3.view(draw)
    Q3.view(io)

sys.exit()
