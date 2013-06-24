#!/usr/bin/env python
# encoding: utf-8
# ddd
# sanity check for submodule working and pushing to...
import sys, petsc4py
petsc4py.init(sys.argv)

import numpy as np

# -------- GLOBAL SCALAR DEFINITIONS -----------------------------

# ======== all definitions are in m,s,g unit system.
n_frames = 30
# ....... dimensions .............................................
x_lower = 0.0
x_upper = 10.0e-6                    # lenght [m]
y_lower = 0.0
y_upper = 10.0e-6                   # notice that for multilayer this is value will be over-written
# ........ material properties ...................................

# vacuum
vac      = np.ones([3])
eo = 8.854187817e-12            # vacuum permittivity   - [F/m]
mo = 4e-7*np.pi                 # vacuum peremeability  - [V.s/A.m]
co = 1.0/np.sqrt(eo*mo)           # vacuum speed of light - [m/s]
zo = np.sqrt(mo/eo)
vac[0] = eo
vac[1] = eo
vac[2] = mo 
# material
mat_shape = 'multilayer'           # material definition: homogeneous, interface, rip (moving perturbation), multilayered
# parameters needed for pml calculation
num_pml = 8
norder = 3
Ro = 1.0e-6

# initialize background refractive index and etas
eta      = np.ones([3])
bkg_n    = np.ones([2])

bkg_n[0] = np.sqrt(eta[0]*eta[2])
bkg_n[1] = np.sqrt(eta[1]*eta[2])


# if interface declare position
x_change = (x_upper-x_lower)/2.0
y_change = (y_upper-y_lower)/2.0

# set moving refractive index or gaussian2D parameters
rip_velocity = np.zeros([2,3])
rip_offset   = np.zeros([2,3])
rip_sigma    = np.zeros([2,3])
delta_eta    = np.zeros([3])

rip_offset[0,:].fill((x_upper-x_lower)/2.0)
rip_offset[1,:].fill((y_upper-y_lower)/2.0)
rip_sigma[0,:].fill((x_upper-x_lower)/25.0)
rip_sigma[1,:].fill((y_upper-y_lower)/25.0)
rip_sigma.fill(10e-6)
rip_sigma.fill(10e-6)
rip_sigma2 = rip_sigma**2

delta_eta = np.zeros([3])
delta_eta = 0.1*eta

# set multilayer parameters
num_materials = 2
num_layers = 10

layers = np.zeros([num_materials,9]) # _layer:  eps mu N t chi2e chi2m chi3e chi3m

layers[0,0:3] = 1.0
layers[1,0:3] = 2.5
layers[0,3] = 10
layers[1,3] = layers[0,3] - 1
layers[0,4] = 15.0e-9
layers[1,4] = 50.0e-9

if mat_shape=='multilayer':
    y_upper = num_layers*np.sum(layers[:,4])+layers[0,4]
    tlp = np.sum(layers[:,4])
    mlp = np.floor(tlp/1e-9)

# set non-linear parameters of the material
chi2 = chi3 = np.zeros( [3], order='F')

# ........ excitation - initial conditoons .......................

# pre-allocate arrays
ex_sigma     = np.ones([3])    # x,y,t
ex_offset    = np.zeros([3])
ex_amplitude = np.ones([3])
ex_kvector   = np.zeros([2])

# fill arrays and set respective values
ex_type   = 'plane'
ex_lambda = 1e-6
ex_sigma[0:1] = 1.0*ex_lambda
ex_sigma[2]   = (y_upper-y_lower)/2.0
ex_offset[2]  = (y_upper-y_lower)/2.0

# post calculations
omega    = 2.0*np.pi*co/ex_lambda   # frequency
k        = 2.0*np.pi/ex_lambda      # k vector magnitude
ex_kvector[0] = k                   # propagation along the x-direction

# ........ Boundary settings settings .................

bc_x_lower = 'scattering'
bc_x_upper = 'none'
bc_y_lower = 'scattering'
bc_y_upper = 'none'

aux_bc_x_lower = 'scattering'
aux_bc_x_upper = 'pml'
aux_bc_y_upper = 'metallic'
aux_bc_y_upper = 'metallic'



# ........ pre-calculations for wave propagation .................

v = co/bkg_n.min()

# Grid - mesh settings
nx = np.floor(20*(x_upper-x_lower)/ex_lambda)
if mat_shape=='multilayer':
    ny = np.floor((y_upper-y_lower)/1e-9)
else:
    ny = np.floor(20*(y_upper-y_lower)/ex_lambda)

ddx = (x_upper-x_lower)/nx
ddy = (y_upper-y_lower)/ny
ddt = dt=0.50/(co*np.sqrt(1.0/(ddx**2)+1.0/(ddy**2)))
dt = ddt
t_final = (x_upper-x_lower)/v
max_steps = np.floor(t_final/ddt)+1
print t_final,ddx,ddy,ddt,max_steps

dxdt = dt/ddx
dydt = dt/ddy


# -------- GLOBAL FUNCTION DEFINITIONS --------------

def etar(da,ddx,ddy,t=0):
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
    t is the time coordinate

    on output aux holds:

                               EM equivalent

         idim   curvilinear  |   TE      TM
         0:     eta_1        |   mu1     eps1
         1:     eta_2        |   mu2     eps2
         2:     eta_3        |   eps3    mu3

    """
    nx, ny = da.getSizes()
    (xi, xf), (yi, yf) = da.getRanges()
    X = np.linspace(xi,xf,xf-xi)*ddx
    Y = np.linspace(yi,yf,yf-yi)*ddy
    print X
    y,x = np.meshgrid(Y,X)
    eta_out = np.ones( [3,len(X),len(Y)], order='F')

    if mat_shape=='gaussian1dx':
        u_x_eta1 = x - rip_velocity[0,0]*t - rip_offset[0,0]
        u_x_eta2 = x - rip_velocity[0,1]*t - rip_offset[0,1]
        u_x_eta3 = x - rip_velocity[0,2]*t - rip_offset[0,2]
        u_y_eta1 = y - rip_velocity[1,0]*t - rip_offset[1,0]
        u_y_eta2 = y - rip_velocity[1,1]*t - rip_offset[1,1]
        u_y_eta3 = y - rip_velocity[1,2]*t - rip_offset[1,2]

        u_eta1 = (u_x_eta1/rip_sigma[0,0])**2 + (u_y_eta1/rip_sigma[1,0])**2
        u_eta2 = (u_x_eta2/rip_sigma[0,1])**2 + (u_y_eta2/rip_sigma[1,1])**2
        u_eta3 = (u_x_eta3/rip_sigma[0,2])**2 + (u_y_eta3/rip_sigma[1,2])**2

        eta_out[0,:,:] = delta_eta[0]*np.exp(-u_eta1) + eta[0]
        eta_out[1,:,:] = delta_eta[1]*np.exp(-u_eta2) + eta[1]
        eta_out[2,:,:] = delta_eta[2]*np.exp(-u_eta3) + eta[2]
    elif mat_shape=='homogeneous':
        eta_out[0,:,:] = eta[0]
        eta_out[1,:,:] = eta[1]
        eta_out[2,:,:] = eta[2]
    elif mat_shape=='interfacex':
        eta_out[0,:,:] = 1*(x<x_change) + 4*(x>=x_change)
        eta_out[1,:,:] = 1*(x<x_change) + 4*(x>=x_change)
        eta_out[2,:,:] = 1*(x<x_change) + 4*(x>=x_change)
    elif mat_shape=='interfacey':
        eta_out[0,:,:] = 1*(y<y_change/2) + 4*(x>=y_change/2)
        eta_out[1,:,:] = 1*(y<y_change/2) + 4*(x>=y_change/2)
        eta_out[2,:,:] = 1*(y<y_change/2) + 4*(x>=y_change/2)
    elif mat_shape=='multilayer':
        yi = np.arange(0,num_layers+1)*tlp
        for i in range(0,num_layers+1):
            for m in range(0,num_materials):
                if m==0:
                    eta_out[0,:,:] += layers[0,0]*(y>=yi[i])*(y<=(yi[i]+layers[0,4]))
                    eta_out[1,:,:] += layers[0,1]*(y>=yi[i])*(y<=(yi[i]+layers[0,4]))
                    eta_out[2,:,:] += layers[0,2]*(y>=yi[i])*(y<=(yi[i]+layers[0,4]))
                else:
                    eta_out[0,:,:] += layers[m,0]*(y>=(yi[i]+layers[m-1,4]))*(y<=(yi[i]+layers[m,4]))
                    eta_out[1,:,:] += layers[m,1]*(y>=(yi[i]+layers[m-1,4]))*(y<=(yi[i]+layers[m,4]))
                    eta_out[2,:,:] += layers[m,2]*(y>=(yi[i]+layers[m-1,4]))*(y<=(yi[i]+layers[m,4]))


        eta_out[0,:,:] = layers[0,0]*(num_layers*tlp<y)*(y<=num_layers*tlp+layers[0,4])
        eta_out[1,:,:] = layers[0,1]*(num_layers*tlp<y)*(y<=num_layers*tlp+layers[0,4])
        eta_out[2,:,:] = layers[0,2]*(num_layers*tlp<y)*(y<=num_layers*tlp+layers[0,4]) 

    eta_out[0,:,:] = eta_out[0,:,:]*vac[0]
    eta_out[1,:,:] = eta_out[1,:,:]*vac[1]
    eta_out[2,:,:] = eta_out[2,:,:]*vac[2]

    return eta_out

def qinit(Q1,Q2,Q3,da):
    """
    Set initial conditions for q in the grid, but not on the boundary lines
    """
    nx, ny = da.getSizes()
    (xi, xf), (yi, yf) = da.getRanges()

    q1 = da.getVecArray(Q1)
    q2 = da.getVecArray(Q2)
    q3 = da.getVecArray(Q3)

    if ex_type=='off':
        X = np.linspace(xi,xf,xf-xi)*ddx
        Y = np.linspace(yi,yf,yf-yi)*ddy
        y,x = np.meshgrid(Y,X)
        dd1 = x_upper-x_lower
        dd2 = y_upper-y_lower
        sdd = ex_lambda
        r2 = (x-dd1/2.0)**2 + (y-dd2/2.0)**2
        q1[:,:] = 0.0
        q2[:,:] = zo*np.exp(-r2/(sdd**2))
        q3[:,:] = 1.0*np.exp(-r2/(sdd**2))
    else:
        q1[:,:] = 0.0
        q2[:,:] = 0.0
        q3[:,:] = 0.0

def qbc(Q1,Q2,Q3,da,t):
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
            bc_scattering(Q1,Q2,Q3,da,t,0,0)
        elif bc_x_lower == 'none':
            pass
    if xi == nx-1:
        if bc_x_upper == 'metallic':
            q1[-1,:] = 0.0
            q2[-1,:] = 0.0
            q3[-1,:] = 0.0
        elif bc_x_upper == 'scattering':
            bc_scattering(Q1,Q2,Q3,da,t,0,1)
        elif bc_x_upper == 'none':
            pass
    if yi == 0:
        if bc_y_lower == 'metallic':
            q1[:,0] = 0.0
            q2[:,0] = 0.0
            q3[:,0] = 0.0
        elif bc_y_lower == 'scattering':
            bc_scattering(Q1,Q2,Q3,da,t,1,0)
        elif bc_y_lower == 'none':
            pass
    if yf == ny-1:
        if bc_y_upper == 'metallic':
            q1[:,-1] = 0.0
            q2[:,-1] = 0.0
            q3[:,-1] = 0.0
        elif bc_y_upper == 'scattering':
            bc_scattering(Q1,Q2,Q3,da,t,1,1)
        elif bc_y_upper == 'none':
            pass

def bc_scattering(Q1,Q2,Q3,da,t,axis,side):
    """
    Boundary scattering conditions (source)
    """
    nx, ny = da.getSizes()
    (xi, xf), (yi, yf) = da.getRanges()
    X = np.linspace(xi,xf,xf-xi)*ddx
    Y = np.linspace(yi,yf,yf-yi)*ddy
    y,x = np.meshgrid(Y,X)
    pulseshape = np.zeros( [len(X),len(Y)], order='F')
    harmonic = np.zeros( [len(X),len(Y)], order='F')
    if ex_type=='plane':
        pulseshape = 1.0
        harmonic = np.sin(omega*t)
    elif ex_type=='gauss-beam':
        pulseshape = np.exp(-(y - ex_offset[2])**2/ex_sigma[2]**2)
        harmonic = np.sin(omega*t)
    elif ex_type=='gauss_pulse':
        pulseshape = np.exp(-(ex_vx*(t-ex_offset[3]))**2/ex_sigma[3]**2 - (y - ex_offset[2] - ex_vy*(t-ex_offset[3]))**2/ex_sigma[2]**2)
        harmonic = np.sin(omega*t)
    elif ex_type=='plane_pulse':
        pulseshape = np.exp(-(ex_vx*(t-ex_offset[3]))**2/ex_sigma[3]**2)
        harmonic = np.sin(omega*t)
    elif ex_type=='simple_pulse2D':
        pulseshape = np.exp(-(ex_vx*(t-ex_offset[3]))**2/ex_sigma[3]**2 - (y - ex_offset[2] - ex_vy*(t-ex_offset[3]))**2/ex_sigma[2]**2)
        harmonic = 1.0
    elif ex_type=='simple_pulse2D_x':
        pulseshape = np.exp(-(ex_vx*(t-ex_offset[3]))**2/ex_sigma[3]**2)
        harmonic = 1.0
    elif ex_type=='off':
        pulseshape = 0.0
        harmonic = 0.0
    elif ex_type=='jump':
        pulseshape = 1.0
        harmonic = 1.0
    
    q1  = da.getVecArray(Q1)
    q2  = da.getVecArray(Q2)
    q3  = da.getVecArray(Q3)

    if (axis,side) == (0,0) and xi == 0:
        q2[0 ,:] = zo*ex_amplitude[1]*pulseshape*harmonic
        q3[0 ,:] = ex_amplitude[1]*pulseshape*harmonic
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
    #    ddx,ddy,dt,co = 1,1,1,1
    build_pml(pml,pml_type,ddx,ddy,dt,norder,Ro,co,xi+1,xf,yi+1,yf,nx,ny)

def aux_bc():
    """
    user controlled boundary auxiliary conditions
    """
    pass

def write(Q1,Q2,Q3,filename):
    io = PETSc.Viewer().createBinary(filename,mode="w")
    Q1.view(io)
    Q2.view(io)
    Q3.view(io)
    io.destroy()

def gauges(da):
    from probe import Probe
    entries=[[ 0, 0],
         [34,10],
         [34,34],
         [10,34],
         [0,15],
         ]
    prb = Probe(da, entries)
    return prb

# -------- MAIN PROGRAM --------------

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
#aux = np.ones ([num_aux,xf-xi,yf-yi], order='F')
aux = etar(da,ddx,ddy)

pml = np.ones ([num_pml,xf-xi,yf-yi], order='F')
pml_axis = 0
pml_side = 1
pml_type = pml_axis*2+pml_side
aux_bc_pml(pml,pml_type,xi,xf,yi,yf,nx,ny)
from matplotlib import pylab
#for i in range(8):
pylab.figure()
pylab.contourf(aux[1,:,:].copy())
pylab.colorbar()
pylab.show()

draw = PETSc.Viewer.DRAW()
for n in range(0,int(max_steps)):
    if n == 0:
        t = n*dt
        qinit(Q1,Q2,Q3,da)
        qbc(Q1,Q2,Q3,da,t)
    t = n*dt
    bc_scattering(Q1,Q2,Q3,da,t,0,0)
    bc_scattering(Q1,Q2,Q3,da,t,1,0)
    da.globalToLocal(Q1, Q1loc)
    da.globalToLocal(Q2, Q2loc)
    fdtd_2d(aux,pml,dxdt,dydt,s1,s2,s3,s4,q1,q2,q3,xi+1,xf,yi+1,yf,gxi+1,gxf,gyi+1,gyf,0,0)
    da.localToGlobal(Q3loc, Q3)

    da.globalToLocal(Q3, Q3loc)
    fdtd_2d(aux,pml,dxdt,dydt,s1,s2,s3,s4,q1,q2,q3,xi+1,xf,yi+1,yf,gxi+1,gxf,gyi+1,gyf,1,0)
    da.localToGlobal(Q1loc, Q1)
    da.localToGlobal(Q2loc, Q2)

    #bc(Q1,Q2,Q3)
    Q3.view(draw)
    #write(Q1,Q3,Q3,"./_test/step%d.dat" % t)
    prb = gauges(da)
    prb.probe('Q1', Q1)
    prb.probe('Q2', Q2)

#prb.save("probe.dat")
#from pprint import pprint
#pprint(prb.cache)
#print prb.cache['Q2'][0,15]

sys.exit()
