#!/usr/bin/env python
# encoding: utf-8
# ddd
# sanity check for submodule working and pushing to...
import sys, petsc4py, os
petsc4py.init(sys.argv)
import pickle as pkl
import numpy as np
from mpi4py import MPI

# -------- GLOBAL SCALAR DEFINITIONS -----------------------------

n_frames    = 300
save_outdir = '/media/noor-labs/fdtd-test/_spheres'
save_name   = 'moving'
liveview    = False
write_q     = True
write_aux   = True
gauge       = False
write_gauge = False
mode        = 'TM'
debug_eta   = True
debug_auxbc = False
vaccum_ones = False
before_step = False
# ======== all definitions are in m,s,g unit system.
# ....... dimensions .............................................
x_lower = -50.0e-6
x_upper = 100.0e-6                    # lenght [m]
y_lower = -30.0e-6
y_upper = 30.0e-6                   # notice that for multilayer this is value will be over-written
mid_x = (x_upper-x_lower)/2.0
mid_y = (y_upper-y_lower)/2.0
# ........ material properties ...................................

# vacuum
vac     = np.ones([3])
if vaccum_ones:
    eo  = 1.0
    mo  = 1.0
else:
    eo  = 8.854187817e-12            # vacuum permittivity   - [F/m]
    mo  = 4e-7*np.pi                 # vacuum peremeability  - [V.s/A.m]

co      = 1.0/np.sqrt(eo*mo)           # vacuum speed of light - [m/s]
zo      = np.sqrt(mo/eo)

# material
mat_shape       = 'custom'           # material definition: homogeneous, interface, rip (moving perturbation), multilayered
mat_nonliner    = False
mat_dispersion  = False
# initialize material properties and fill with default values (this should become class material)
if mat_shape == 'interface' or 'interfacex' or 'interfacey':
    eta             = np.ones([2,3])
    eta[0,:]        = 1.0
    eta[1,:]        = 2.0
    if mat_shape == 'interfacex':
        mat_change  = (x_upper-x_lower)/2.0
    else:
        mat_change  = (y_upper-y_lower)/2.0

if mat_shape == 'gaussian_x' or 'gaussian_y' or 'gaussian':
    print 'gaussian'
    eta             = np.ones([3])
    delta_eta       = np.zeros([3])
    eta_velocity    = np.zeros([2,3])
    eta_offset      = np.zeros([2,3])
    eta_sigma       = np.zeros([2,3])
    # once the material class is created the settings below should be created as defaults
    eta             = eta*1.5
    delta_eta       = 0.1*eta
    eta_offset[0,:].fill(10e-6)
    eta_offset[1,:].fill(mid_y)
    eta_velocity[0,:].fill(0.61*co)
    eta_velocity[1,:].fill(0.0)
    print eta_velocity
    # eta_sigma[0,:].fill(5.0e-6)#(x_upper-x_lower)/25.0)
    # eta_sigma[1,:].fill(5.0e-6)#(y_upper-y_lower)/25.0)
    # eta_sigma.fill(5.0e-6)
    eta_sigma.fill(5.0e-6)

if mat_shape =='multilayer':
    num_materials   = 2
    num_layers      = 10
    layers          = np.zeros([num_materials,9]) # _layer:  eps mu N t chi2e chi2m chi3e chi3m
    layers[0,0:3]   = 1.0
    layers[1,0:3]   = 2.0
    layers[0,3]     = 10
    layers[1,3]     = layers[0,3] - 1
    layers[0,4]     = 15.0e-9
    layers[1,4]     = 50.0e-9

if mat_shape =='homogeneous':
    eta             = np.ones([3])*1.5

if mat_shape =='custom':
    num_materials   = 2
    num_layers      = 50
    layers          = np.zeros([num_materials,9]) # _layer:  eps mu N t chi2e chi2m chi3e chi3m
    layers[0,0:3]   = 1.0
    layers[1,0:3]   = 2.0
    layers[0,3]     = 10
    layers[1,3]     = layers[0,3] - 1
    layers[0,4]     = 10.0e-9
    layers[1,4]     = 50.0e-9

if mat_nonliner:
    chi2 = chi3 = np.zeros( [3], order='F')
if mat_dispersion:
    num_poles = 2


# ........ excitation - initial conditoons .......................
# pre-allocate arrays
ex_sigma        = np.ones([3])    # x,y,t
ex_offset       = np.zeros([3])
ex_amplitude    = np.ones([3])
ex_kvector      = np.zeros([2])

# fill arrays and set respective values
ex_type         = 'plane_pulse'
ex_lambda       = 1.0e-6
ex_sigma[0:1]   = 1.0*ex_lambda
ex_sigma[2]     = (y_upper-y_lower)/2.0
ex_offset[2]    = (y_upper-y_lower)/2.0

# ........ Boundary settings settings .................
bc_lower = ['scattering', 'scattering'] # left and bottom boundary contidions
bc_upper = ['none', 'none'] # right and top

aux_bc_lower = ['none', 'none'] # left and bottom boundary contidions
aux_bc_upper = ['pml', 'none'] # right and top

# ........ pre-calculations for wave propagation .................

# Wave propagation calculations
omega    = 2.0*np.pi*co/ex_lambda   # frequency
k        = 2.0*np.pi/ex_lambda      # k vector magnitude
ex_kvector[0] = k                   # propagation along the x-direction

# get the minimum speed in the medium
v = np.zeros([2])
if mat_shape == 'gaussian_x':
    v[0] = co/(eta.max()+delta_eta.max())
    v[1] = co/(eta.min())
elif mat_shape == 'multilayer':
    v[0] = co/(layers[:,0:3].max())
    v[1] = co/(layers[:,0:3].min())
elif mat_shape == 'interface' or 'interfacex' or 'interfacey':
    v[0] = co/(eta.max())
    v[1] = co/(eta.min())
elif mat_shape == 'homogeneous':
    v[0] = v[1] = co/(eta.max())
elif mat_shape == 'custom':
    v[0] = co/1.5#/(layers[:,0:3].max())
    v[1] = co/1.5#/(layers[:,0:3].min())


# Grid - mesh settings
nx = np.floor(50*(x_upper-x_lower)/ex_lambda)

if mat_shape=='multilayer':
    y_upper = num_layers*np.sum(layers[:,4])+layers[0,4]
    tlp = np.sum(layers[:,4])
    mlp = np.floor(tlp/1e-9)
    ny = np.floor((y_upper-y_lower)/2.5e-9)
elif mat_shape=='customale':
    y_upper = num_layers*np.sum(layers[:,4])+layers[0,4]
    tlp = np.sum(layers[:,4])
    mlp = np.floor(tlp/1.0e-9)
    ny = np.floor((y_upper-y_lower)/1.0e-9)
else:
    ny = np.floor(50*(y_upper-y_lower)/ex_lambda)

ddx = (x_upper-x_lower)/nx
ddy = (y_upper-y_lower)/ny
ddt = dt=0.20/(co*np.sqrt(1.0/(ddx**2)+1.0/(ddy**2)))



dt = ddt
t_final = (x_upper-x_lower)/v.min()
max_steps = np.floor(t_final/ddt)+1
n_write = np.floor(max_steps/n_frames)
if MPI.COMM_WORLD.rank == 0:    
    print 'speed velocity is ', v.min()
    print 'total distance in x-direction: ',x_upper-x_lower
    print 'ddx,ddy  ', ddx, ddy
    print 'nx, ny ', nx,ny
    print 't_final ', t_final,
    print 'max_steps ', max_steps

dxdt = dt/ddx
dydt = dt/ddy

if mode == 'TM':
    vac[0]  = eo
    vac[1]  = eo
    vac[2]  = mo
else:
    print 'TE mode not implemented ---  self-destructing (!)'
    1/0


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
    X = np.linspace(xi,xf,xf-xi)*ddx + x_lower
    Y = np.linspace(yi,yf,yf-yi)*ddy + y_lower
    y,x = np.meshgrid(Y,X)
    eta_out = np.zeros( [3,len(X),len(Y)], order='F')
    eta_temp = eta_out.copy()

    if mat_shape=='gaussian_x':
        u_x_eta1 = x - eta_velocity[0,0]*t - eta_offset[0,0]
        u_x_eta2 = x - eta_velocity[0,1]*t - eta_offset[0,1]
        u_x_eta3 = x - eta_velocity[0,2]*t - eta_offset[0,2]

        u_eta1 = (u_x_eta1/eta_sigma[0,0])**2
        u_eta2 = (u_x_eta2/eta_sigma[0,1])**2
        u_eta3 = (u_x_eta3/eta_sigma[0,2])**2

        eta_out[0,:,:] = delta_eta[0]*np.exp(-u_eta1) + eta[0]
        eta_out[1,:,:] = delta_eta[1]*np.exp(-u_eta2) + eta[1]
        eta_out[2,:,:] = delta_eta[2]*np.exp(-u_eta3) + eta[2]
    elif mat_shape=='gaussian_y':
        u_y_eta1 = y - eta_velocity[1,0]*t - eta_offset[1,0]
        u_y_eta2 = y - eta_velocity[1,1]*t - eta_offset[1,1]
        u_y_eta3 = y - eta_velocity[1,2]*t - eta_offset[1,2]

        u_eta1 = (u_y_eta1/eta_sigma[1,0])**2
        u_eta2 = (u_y_eta2/eta_sigma[1,1])**2
        u_eta3 = (u_y_eta3/eta_sigma[1,2])**2

        eta_out[0,:,:] = delta_eta[0]*np.exp(-u_eta1) + eta[0]
        eta_out[1,:,:] = delta_eta[1]*np.exp(-u_eta2) + eta[1]
        eta_out[2,:,:] = delta_eta[2]*np.exp(-u_eta3) + eta[2]
    elif mat_shape=='gaussian':
        u_x_eta1 = x - eta_velocity[0,0]*t - eta_offset[0,0]
        u_x_eta2 = x - eta_velocity[0,1]*t - eta_offset[0,1]
        u_x_eta3 = x - eta_velocity[0,2]*t - eta_offset[0,2]
        u_y_eta1 = y - eta_velocity[1,0]*t - eta_offset[1,0]
        u_y_eta2 = y - eta_velocity[1,1]*t - eta_offset[1,1]
        u_y_eta3 = y - eta_velocity[1,2]*t - eta_offset[1,2]

        u_eta1 = (u_x_eta1/eta_sigma[0,0])**2 + (u_y_eta1/eta_sigma[1,0])**2
        u_eta2 = (u_x_eta2/eta_sigma[0,1])**2 + (u_y_eta2/eta_sigma[1,1])**2
        u_eta3 = (u_x_eta3/eta_sigma[0,2])**2 + (u_y_eta3/eta_sigma[1,2])**2

        eta_out[0,:,:] = delta_eta[0]*np.exp(-u_eta1) + eta[0]
        eta_out[1,:,:] = delta_eta[1]*np.exp(-u_eta2) + eta[1]
        eta_out[2,:,:] = delta_eta[2]*np.exp(-u_eta3) + eta[2]
    elif mat_shape=='homogeneous':
        eta_out[0,:,:] = eta[0]
        eta_out[1,:,:] = eta[1]
        eta_out[2,:,:] = eta[2]
    elif mat_shape=='interfacex':
        print eta
        eta_out[0,:,:] = eta[0,0]*(x<mat_change) + eta[1,0]*(x>=mat_change)
        eta_out[1,:,:] = eta[0,1]*(x<mat_change) + eta[1,1]*(x>=mat_change)
        eta_out[2,:,:] = eta[0,2]*(x<mat_change) + eta[1,2]*(x>=mat_change)
    elif mat_shape=='interfacey':
        eta_out[0,:,:] = eta[0,0]*(y<mat_change) + eta[1,0]*(y>=mat_change)
        eta_out[1,:,:] = eta[0,1]*(y<mat_change) + eta[1,1]*(y>=mat_change)
        eta_out[2,:,:] = eta[0,2]*(y<mat_change) + eta[1,2]*(y>=mat_change)
    elif mat_shape=='multilayer':
        yi = np.arange(0,num_layers+1)*tlp
        for m in range(0,num_materials):
            for i in range(0,num_layers+1):
                if m==0:
                    eta_temp[0,:,:] = eta_out[0,:,:] + layers[0,0]*(y>=yi[i])*(y<=(yi[i]+layers[0,4]))
                    eta_temp[1,:,:] = eta_out[1,:,:] + layers[0,1]*(y>=yi[i])*(y<=(yi[i]+layers[0,4]))
                    eta_temp[2,:,:] = eta_out[2,:,:] + layers[0,2]*(y>=yi[i])*(y<=(yi[i]+layers[0,4]))
                else:
                    eta_temp[0,:,:] = eta_out[0,:,:] + layers[m,0]*(y>(yi[i]+layers[m-1,4]))*(y<(yi[i]+layers[m-1,4]+layers[m,4]))
                    eta_temp[1,:,:] = eta_out[1,:,:] + layers[m,1]*(y>(yi[i]+layers[m-1,4]))*(y<(yi[i]+layers[m-1,4]+layers[m,4]))
                    eta_temp[2,:,:] = eta_out[2,:,:] + layers[m,2]*(y>(yi[i]+layers[m-1,4]))*(y<(yi[i]+layers[m-1,4]+layers[m,4]))

                eta_out = eta_temp.copy()
    elif mat_shape=='custom':
        r       = np.asarray([1e-6,2e-6,3e-6,5e-6,8e-6,13e-6,21e-6]); 
        # theta   = np.asarray([0.0,np.pi/8.0,np.pi/4.0,3.0*np.pi/8.0,np.pi/2.0])
        # phi     = np.asarray([0.0,np.pi/8.0,np.pi/4.0,3.0*np.pi/8.0,np.pi/2.0,5.0*np.pi/8.0,3.0*np.pi/4.0,7.0*np.pi/8.0,np.pi,9.0*np.pi/8.0,5.0*np.pi/4.0,11.0*np.pi/8.0,3.0*np.pi/2.0,13.0*np.pi/8.0,7.0*np.pi/4.0,15.0*np.pi/8.0])
        phi     = np.asarray([0.0,np.pi/8.0,np.pi/4.0,3.0*np.pi/8.0,np.pi/2.0,5.0*np.pi/8.0,3.0*np.pi/4.0,7.0*np.pi/8.0,np.pi])
        r_sph   = 0.5e-6

        for j in range(0,r.size):
                for k in range(0,phi.size):
                    x_off = r[j]*np.cos(phi[k]-np.pi/2);
                    y_off = r[j]*np.sin(phi[k]-np.pi/2);
                    eta_out = eta_out  + 1*(np.sqrt((x-x_off)**2+(y-y_off)**2)<r_sph)
        eta_out = eta_out + 1
        # # this material is based on the multilayer but has a modification to it to include the substrate and finite length rods
        # yi = np.arange(0,num_layers+1)*tlp
        # for m in range(0,num_materials):
        #     for i in range(0,num_layers+1):
        #         if m==0:
        #             eta_temp[0,:,:] = eta_out[0,:,:] + layers[0,0]*(y>=yi[i])*(y<=(yi[i]+layers[0,4]))*(x>=mid_x)*(x<=(mid_x+10e-6))
        #             eta_temp[1,:,:] = eta_out[1,:,:] + layers[0,1]*(y>=yi[i])*(y<=(yi[i]+layers[0,4]))*(x>=mid_x)*(x<=(mid_x+10e-6))
        #             eta_temp[2,:,:] = eta_out[2,:,:] + layers[0,2]*(y>=yi[i])*(y<=(yi[i]+layers[0,4]))*(x>=mid_x)*(x<=(mid_x+10e-6))
        #         else:
        #             eta_temp[0,:,:] = eta_out[0,:,:] + layers[m,0]*(y>(yi[i]+layers[m-1,4]))*(y<(yi[i]+layers[m-1,4]+layers[m,4]))*(x>=mid_x)*(x<=(mid_x+10e-6))
        #             eta_temp[1,:,:] = eta_out[1,:,:] + layers[m,1]*(y>(yi[i]+layers[m-1,4]))*(y<(yi[i]+layers[m-1,4]+layers[m,4]))*(x>=mid_x)*(x<=(mid_x+10e-6))
        #             eta_temp[2,:,:] = eta_out[2,:,:] + layers[m,2]*(y>(yi[i]+layers[m-1,4]))*(y<(yi[i]+layers[m-1,4]+layers[m,4]))*(x>=mid_x)*(x<=(mid_x+10e-6))

        #         eta_out = eta_temp.copy()
        
        # eta_out = eta_out + 1.4*(x<mid_x)
        # eta_out[eta_out==0.0] = 1.0


        # eta_out[0,:,:] += layers[0,0]*(num_layers*tlp<y)*(y<=num_layers*tlp+layers[0,4])
        # eta_out[1,:,:] += layers[0,1]*(num_layers*tlp<y)*(y<=num_layers*tlp+layers[0,4])
        # eta_out[2,:,:] += layers[0,2]*(num_layers*tlp<y)*(y<=num_layers*tlp+layers[0,4]) 

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
        X = np.linspace(xi,xf,xf-xi)*ddx + x_lower
        Y = np.linspace(yi,yf,yf-yi)*ddy + y_lower
        y,x = np.meshgrid(Y,X)
        dd1 = -20e-6#x_upper-x_lower
        dd2 = y_upper-y_lower
        sdd = 1.0*ex_lambda
        r2 = (x-dd1/2.0)**2 #+ (y-dd2/2.0)**2
        q1[:,:] = 0.0
        q2[:,:] = 10.0*np.exp(-r2/(sdd**2))
        q3[:,:] = 0.0*np.exp(-r2/(sdd**2))
    else:
        q1[:,:] = 0.0
        q2[:,:] = 0.0
        q3[:,:] = 0.0

def qbc(Q1,Q2,Q3,da,t):
    """
    Set the boundary conditions for q. Implemented conditions are:

    .. metallic:     set all qs to 0
    .. scattering:   scattering boundary conditions (line), set by function bc_scattering (user controlled)
    .. pml
    .. periodic:     periodic boundary conditions            (not implemented)
    .. neumann:      neumann boundary conditions             (not implemented)
    .. rounded:      end boundary becomes beginning boundary (not implemented)
    .. none:         pass, does not set any value
    """
    nx, ny = da.getSizes()
    (xi, xf), (yi, yf) = da.getRanges()

    q1  = da.getVecArray(Q1)
    q2  = da.getVecArray(Q2)
    q3  = da.getVecArray(Q3)

    if xi == 0:
        if bc_lower[0] == 'metallic':
            q1[0,:] = 0.0
            q2[0,:] = 0.0
            q3[0,:] = 0.0
        elif bc_lower[0] == 'scattering':
            bc_scattering(Q1,Q2,Q3,da,t,0,0)
        elif bc_lower[0] == 'custom':
            bc_custom()           
        elif bc_lower[0] == 'none':
            pass
    if xi == nx-1:
        if bc_upper[0] == 'metallic':
            q1[-1,:] = 0.0
            q2[-1,:] = 0.0
            q3[-1,:] = 0.0
        elif bc_upper[0] == 'scattering':
            bc_scattering(Q1,Q2,Q3,da,t,0,1)
        elif bc_upper[0] == 'custom':
            bc_custom()
        elif bc_upper[0] == 'none':
            pass
    if yi == 0:
        if bc_lower[1] == 'metallic':
            q1[:,0] = 0.0
            q2[:,0] = 0.0
            q3[:,0] = 0.0
        elif bc_lower[1] == 'scattering':
            bc_scattering(Q1,Q2,Q3,da,t,1,0)
        elif bc_lower[1] == 'custom':
            bc_custom()
        elif bc_lower[1] == 'none':
            pass
    if yf == ny-1:
        if bc_upper[1] == 'metallic':
            q1[:,-1] = 0.0
            q2[:,-1] = 0.0
            q3[:,-1] = 0.0
        elif bc_upper[1] == 'scattering':
            bc_scattering(Q1,Q2,Q3,da,t,1,1)
        elif bc_upper[1] == 'custom':
            bc_custom()
        elif bc_upper[1] == 'none':
            pass

def auxbc(da):
    num_pml = 8
    nx, ny = da.getSizes()
    (xi, xf), (yi, yf) = da.getRanges()
    temp_aux_bc = np.ones ([num_pml,xf-xi,yf-yi], order='F')

    if aux_bc_lower[0] == 'pml':
        pml_axis = 0
        pml_side = 0
        pml_type = pml_axis*2+pml_side
        aux_bc_pml(temp_aux_bc,pml_type,xi,xf,yi,yf,nx,ny)
    elif aux_bc_lower[0] == 'custom':
        aux_bc_custom(temp_aux_bc)
    else:
        pass

    if aux_bc_lower[1] == 'pml':
        pml_axis = 1
        pml_side = 0
        pml_type = pml_axis*2+pml_side
        aux_bc_pml(temp_aux_bc,pml_type,xi,xf,yi,yf,nx,ny)
    elif aux_bc_lower[1] == 'custom':
        aux_bc_custom(temp_aux_bc)
    else:
        pass

    if aux_bc_upper[0] == 'pml':
        pml_axis = 0
        pml_side = 1
        pml_type = pml_axis*2+pml_side
        aux_bc_pml(temp_aux_bc,pml_type,xi,xf,yi,yf,nx,ny)
    elif aux_bc_upper[0] == 'custom':
        aux_bc_custom(temp_aux_bc)
    else:
        pass

    if aux_bc_upper[1] == 'pml':
        pml_axis = 1
        pml_side = 1
        pml_type = pml_axis*2+pml_side
        aux_bc_pml(temp_aux_bc,pml_type,xi,xf,yi,yf,nx,ny)
    elif aux_bc_upper[1] == 'custom':
        aux_bc_custom(temp_aux_bc)
    else:
        pass

    return temp_aux_bc

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
        pulseshape = np.exp(-((t-100.0*dt))**2/(6.6e-15)**2)
        harmonic = 1.0#np.sin(omega*t)
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
    # parameters needed for pml calculation
    norder  = 3
    Ro      = 1.0e-6
    #    ddx,ddy,dt,co = 1,1,1,1
    build_pml(pml,pml_type,ddx,ddy,dt,norder,Ro,co,xi+1,xf,yi+1,yf,nx,ny)

def bc_custom():
    """
    user controlled boundary auxiliary conditions
    """
    pass

def aux_bc_custom(temp_aux_bc):
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
if mat_dispersion:
    from fdtd_da_pml import fdtddispersion2d as fdtd_2d
    from fdtd_da_pml import calcdispersion2d
else:
    from fdtd_da_pml import fdtd2d as fdtd_2d

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

if mat_dispersion:
    p1    = np.zeros([num_poles,3,xf-xi,yf-yi], order='F')
    p2    = np.zeros([num_poles,3,xf-xi,yf-yi], order='F')
    p3    = np.zeros([num_poles,3,xf-xi,yf-yi], order='F')
    psum  = np.zeros([3,xf-xi,yf-yi], order='F')
 
    c1    = np.zeros([num_poles,xf-xi,yf-yi], order='F')
    c2    = np.zeros([num_poles,xf-xi,yf-yi], order='F')
    c3    = np.zeros([num_poles,xf-xi,yf-yi], order='F')

aux = etar(da,ddx,ddy)
aux_bc = auxbc(da)

try:
    os.makedirs(save_outdir)
except: OSError("directory already exist")

if debug_eta:
    if MPI.COMM_WORLD.rank == 0:
        print 'debug eta'
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as pylab
    pylab.figure()
    pylab.imshow(aux[1,:,:].transpose())
    pylab.colorbar()
    pylab.draw()
    if MPI.COMM_WORLD.rank == 0:
        print 'saving'
    save_fig_name = 'aux_'+save_name+str(MPI.COMM_WORLD.rank)+'.png'
    pylab.savefig(os.path.join(save_outdir,save_fig_name))


if debug_auxbc:
    from matplotlib import pylab
    pylab.figure()
    pylab.pcolor(aux_bc[0,:,:].copy())
    pylab.colorbar()
    pylab.show()

if liveview:
    draw = PETSc.Viewer.DRAW()



# create a temporary dictionary with the parameters simulation
if MPI.COMM_WORLD.rank == 0:
    if mat_shape=='gaussian':
        params  = { 'outdir':save_outdir,
                    'nx': nx,
                    'ny': ny,
                    'dt': dt,
                    'dx': ddx,
                    'dy': ddy,
                    'num_steps': max_steps,
                    't_final': t_final,
                    'dimensions': [x_lower,x_upper,y_lower,y_upper],
                    'shape': mat_shape,
                    'ex_type': ex_type,
                    'lambda': ex_lambda,
                    'eta': eta,
                    'eta_velocity':eta_velocity,
                    'bc_lower': bc_lower,
                    'bc_upper': bc_upper,
                    'aux_bc_lower': aux_bc_lower,
                    'aux_bc_upper': aux_bc_upper
                  }
    else:
        params  = { 'outdir':save_outdir,
                    'nx': nx,
                    'ny': ny,
                    'dt': dt,
                    'dx': ddx,
                    'dy': ddy,
                    'num_steps': max_steps,
                    't_final': t_final,
                    'dimensions':[x_lower,x_upper,y_lower,y_upper],
                    'shape': mat_shape,
                    'ex_type': ex_type,
                    'lambda': ex_lambda,
                    'eta' : eta,
                    'bc_lower': bc_lower,
                    'bc_upper': bc_upper,
                    'aux_bc_lower': aux_bc_lower,
                    'aux_bc_upper': aux_bc_upper
                  }

    pkl_out = open(os.path.join(save_outdir,save_name+'.pkl'), 'wb')
    pkl.dump(params, pkl_out)
    pkl_out.close()

ki = int(0)
for n in range(0,int(max_steps)):
    if n == 0:
        t = n*dt
        qinit(Q1,Q2,Q3,da)
        qbc(Q1,Q2,Q3,da,t)
    t = n*dt

    if before_step:
        aux = etar(da,ddx,ddy,t)

    bc_scattering(Q1,Q2,Q3,da,t,0,0)
    bc_scattering(Q1,Q2,Q3,da,t,1,0)
    if mat_dispersion:
        da.globalToLocal(Q1, Q1loc)
        da.globalToLocal(Q2, Q2loc)
        fdtd_2d(aux,aux_bc,dxdt,dydt,s1,s2,s3,s4,q1,q2,q3,psum,xi+1,xf,yi+1,yf,gxi+1,gxf,gyi+1,gyf,0,1)
        da.localToGlobal(Q3loc, Q3)

        da.globalToLocal(Q3, Q3loc)
        fdtd_2d(aux,aux_bc,dxdt,dydt,s1,s2,s3,s4,q1,q2,q3,psum,xi+1,xf,yi+1,yf,gxi+1,gxf,gyi+1,gyf,1,1)
        da.localToGlobal(Q1loc, Q1)
        da.localToGlobal(Q2loc, Q2)

        CalcDispersion2D(q1,q2,q3,c1,c2,c3,p1,p2,p3,psum,xi+1,xf,yi+1,yf,gxi+1,gxf,gyi+1,gyf,0,1)
        CalcDispersion2D(q1,q2,q3,c1,c2,c3,p1,p2,p3,psum,xi+1,xf,yi+1,yf,gxi+1,gxf,gyi+1,gyf,1,1)

    else:
        da.globalToLocal(Q1, Q1loc)
        da.globalToLocal(Q2, Q2loc)
        fdtd_2d(aux,aux_bc,dxdt,dydt,s1,s2,s3,s4,q1,q2,q3,xi+1,xf,yi+1,yf,gxi+1,gxf,gyi+1,gyf,0,1)
        da.localToGlobal(Q3loc, Q3)

        da.globalToLocal(Q3, Q3loc)
        fdtd_2d(aux,aux_bc,dxdt,dydt,s1,s2,s3,s4,q1,q2,q3,xi+1,xf,yi+1,yf,gxi+1,gxf,gyi+1,gyf,1,1)
        da.localToGlobal(Q1loc, Q1)
        da.localToGlobal(Q2loc, Q2)


    if liveview:
        Q3.view(draw)
    if np.mod(n,n_write)==0:
        if write_q:
            save_q_name = 'step'+str(ki).zfill(len(str(int(n_frames)))+2)+'.bin'
            write(Q1,Q3,Q3,os.path.join(save_outdir,save_q_name))
            ki = ki+1
            if MPI.COMM_WORLD.rank == 0:
                print 100*n/max_steps
    if gauge:
        prb = gauges(da)
        prb.probe('Q2', Q1)
        prb.probe('Q3', Q2)

if write_gauge and gauge:
    prb.save(save_outdir+"probe.dat")

sys.exit()
