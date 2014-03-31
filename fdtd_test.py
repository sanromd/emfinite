#!/usr/bin/env python
# encoding: utf-8
# ddd
# sanity check for submodule working and pushing to...
import sys, petsc4py, os
petsc4py.init(sys.argv)
import pickle as pkl
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

# -------- GLOBAL SCALAR DEFINITIONS -----------------------------

n_frames        = 10
save_outdir     = './homework/test'
save_name       = 'sphere'
liveview        = False
write_q         = True
write_aux       = True
gauge           = False
write_gauge     = False
mode            = 'TM'
debug_eta       = True
debug_auxbc     = False
vaccum_ones     = False
before_step     = True
mat_shape       = 'homogeneous'           # material definition: homogeneous, interface, rip (moving perturbation), multilayered
mat_nonliner    = False
mat_dispersion  = False
info            = True
cfl_desired     = 0.9
# ======== all definitions are in m,s,g unit system.

# ....... dimensions .............................................
x_lower = 0.0
x_upper = 20.0e-6                    # lenght [m]
y_lower = 0.0
y_upper = 20.0e-6                   # notice that for multilayer this is value will be over-written
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


eta = np.ones([3])*1.0


# ........ excitation - initial conditoons .......................
# pre-allocate arrays
ex_sigma        = np.ones([3])    # x,y,t
ex_offset       = np.zeros([3])
ex_amplitude    = np.ones([3])
ex_kvector      = np.zeros([2])

# fill arrays and set respective values
ex_type         = 'plane'
ex_lambda       = 1.0e-6
ex_sigma[0:1]   = 1.0*ex_lambda
ex_sigma[2]     = (y_upper-y_lower)/2.0
ex_offset[2]    = (y_upper-y_lower)/2.0

# ........ Boundary settings settings .................

bc_lower = ['scattering', 'scattering'] # left and bottom boundary contidions
bc_upper = ['none', 'none'] # right and bootom

aux_bc_lower = ['none', 'none'] # left and bottom boundary contidions
aux_bc_upper = ['pml', 'none'] # right and bootom

# ........ pre-calculations for wave propagation .................

# Wave propagation calculations
omega    = 2.0*np.pi*co/ex_lambda   # frequency
k        = 2.0*np.pi/ex_lambda      # k vector magnitude
ex_kvector[0] = k                   # propagation along the x-direction

# get the minimum speed in the medium
v = np.zeros([2])
v[0] = v[1] = co/(eta.max())

if mode == 'TM':
    vac[0] = eo
    vac[1] = eo
    vac[2] = mo
elif mode == 'TE':
    vac[0] = mo
    vac[1] = mo
    vac[2] = eo


# -------- GLOBAL FUNCTION DEFINITIONS --------------

def etar(da,t=0):
    """
    eta = etar(num_aux,xi,xf,yi,yf,ddx,ddy)

    Sets the auxiliary arrays for permittivity and permeability.

    Implemented mappings

    ..gaussian:     stationary and moving gaussian shape for eps and mu
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
    y,x = np.meshgrid(Y,X)
    eta_out = np.zeros( [3,len(X),len(Y)], order='F')
    eta_temp = eta_out.copy()

    eta_out[0,:,:] = eta[0]
    eta_out[1,:,:] = eta[1]
    eta_out[2,:,:] = eta[2]
   
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
        sdd = 5.0*ex_lambda
        r2 = (x-dd1/2.0)**2 #+ (y-dd2/2.0)**2
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
        pulseshape = np.exp(-((t-200.0*dt))**2/(50.0*dt)**2)
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

def pec_sphere(Q1,Q2,Q3,da):
    """
    Create sphere of PEC
    """
    nx, ny = da.getSizes()
    (xi, xf), (yi, yf) = da.getRanges()

    q1 = da.getVecArray(Q1)
    q2 = da.getVecArray(Q2)
    q3 = da.getVecArray(Q3)


    X = np.linspace(xi*ddx,xf*ddx,xf-xi)
    Y = np.linspace(yi*ddy,yf*ddy,yf-yi)
    y,x = np.meshgrid(Y,X)
    r2 = (x - mid_x)**2 + (y - mid_y)**2
    r = np.sqrt(r2)
    q1[np.where(r<=3e-6)[0]+xi,np.where(r<=3e-6)[1]+yi] = 0.0
    q2[np.where(r<=3e-6)[0]+xi,np.where(r<=3e-6)[1]+yi] = 0.0 #zo*np.exp(-r2/(sdd**2))
    q3[np.where(r<=3e-6)[0]+xi,np.where(r<=3e-6)[1]+yi] = 0.0 #1.0*np.exp(-r2/(sdd**2))

# -------- MAIN PROGRAM --------------

error = np.zeros([3,10,9])
error[:,:,0] = 2**np.arange(4,14,1)
print error
timeVec1 = PETSc.Vec().createWithArray([0])
timeVec2 = PETSc.Vec().createWithArray([0])
timeVec3 = PETSc.Vec().createWithArray([0])
tic1 = MPI.Wtime()
for test in np.arange(4,14,1):
    nx = ny = 2**test
    ddx = (x_upper-x_lower)/nx
    ddy = (y_upper-y_lower)/ny
    ddt = dt = cfl_desired/(co*np.sqrt(1.0/(ddx**2)+1.0/(ddy**2)))
    dxdt = dt/ddx
    dydt = dt/ddy

    t_final   = (x_upper-x_lower)/v.min()
    max_steps = np.floor(t_final/ddt)+1
    n_write   = np.floor(max_steps/n_frames)

    if MPI.COMM_WORLD.rank==0:
        print 'wave velocity is          ', v.min()
        print 'x-propagation distance    ', x_upper-x_lower
        print 'grid size, (ddx, ddy)     ', ddx, ddy
        print 'grid points (nx, ny)      ', nx,ny
        print 'total propagation time    ', t_final,
        print 'number of time steps      ', max_steps
    # create DA and allocate global and local variables

    if mat_dispersion:
        from fdtd_da_pml import fdtddispersion2d as fdtd_2d
        from fdtd_da_pml import calcdispersion2d
    else:
        from fdtd_da_pml import fdtd2d as fdtd_2d

    tic2 = MPI.Wtime()

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
    else:
        pass

    aux = etar(da)
    aux_bc = auxbc(da)

    try:
        os.makedirs(save_outdir+'_'+str(test))
    except: OSError("directory already exist")

    if debug_eta:
        if MPI.COMM_WORLD.rank==0:
            print 'printing debug eta'
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as pylab
            pylab.figure()
            pylab.imshow(aux[1,:,:].transpose())
            pylab.colorbar()
            pylab.draw()
            print 'saving'
            pylab.savefig(os.path.join(save_outdir+'_'+str(test),'aux_'+save_name+'.png'))

    if debug_auxbc:
        from matplotlib import pylab
        pylab.figure()
        pylab.pcolor(aux_bc[0,:,:].copy())
        pylab.colorbar()
        pylab.draw()
        pylab.savefig(os.path.join(save_outdir+'_'+str(test),'aux_bc_'+save_name+'.png'))

    if liveview:
        draw = PETSc.Viewer.DRAW()



    # create a temporary dictionary with the parameters simulation
    if MPI.COMM_WORLD.rank==0:
        params  = { 'outdir':save_outdir+'_'+str(test),
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

        pkl_out = open(os.path.join(save_outdir+'_'+str(test),save_name+'.pkl'), 'wb')
        pkl.dump(params, pkl_out)
        pkl_out.close()

    tic3 = MPI.Wtime()
    ki = int(0)
    for n in range(0,int(max_steps)):
        if n == 0:
            t = n*dt
            qinit(Q1,Q2,Q3,da)
            qbc(Q1,Q2,Q3,da,t)

        t = n*dt

        if before_step:
            pec_sphere(Q1,Q2,Q3,da)
            #aux = etar(da,ddx,ddy,t)


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
                write(Q1,Q3,Q3,os.path.join(save_outdir+'_'+str(test),save_q_name))
                ki = ki+1

            if info:
                if MPI.COMM_WORLD.rank == 0:
                    print 'percentage of simulation completed: ', str(np.floor(100.0*n/max_steps))[0:3], ' at time  ', str(t)[0:10]

        if gauge:
            prb = gauges(da)
            prb.probe('Q2', Q2)
            prb.probe('Q3', Q3)
    
    toc = MPI.Wtime()
    timeVec2.array = toc - tic2
    timeVec3.array = toc - tic3
    duration2 = timeVec2.max()[1]
    duration3 = timeVec3.max()[1]
    if MPI.COMM_WORLD.rank==0:
        print 'time of computation was: ', str(duration3),'\n'
        print '===================================================== \n'

    if test>4:
        q1r = q1[::2,::2] + q1[1::2,::2] + q1[::2,1::2] + q1[1::2,1::2]
        q1r = q1r/4.0
        q2r = q2[::2,::2] + q2[1::2,::2] + q2[::2,1::2] + q2[1::2,1::2]
        q2r = q2r/4.0
        q3r = q3[::2,::2] + q3[1::2,::2] + q3[::2,1::2] + q3[1::2,1::2]
        q3r = q3r/4.0
        error[:,test-4,1] = ddx
        error[:,test-4,2] = ddy
        error[:,test-4,8] = duration3
        error[:,test-4,6] = dt
        error[:,test-4,7] = max_steps
        error[0,test-4,3] = np.linalg.norm(q1_old - q1r,1)
        error[0,test-4,4] = np.linalg.norm(q1r,1)
        error[0,test-4,5] = ddx*ddy*error[0,test-4,3]/error[0,test-4,4]

        error[1,test-4,3] = np.linalg.norm(q2_old - q2r,1)
        error[1,test-4,4] = np.linalg.norm(q2r,1)
        error[1,test-4,5] = ddx*ddy*error[0,test-4,3]/error[0,test-4,4]

        error[2,test-4,3] = np.linalg.norm(q3_old - q3r,1)
        error[2,test-4,4] = np.linalg.norm(q3r,1)
        error[2,test-4,5] = ddx*ddy*error[0,test-4,3]/error[0,test-4,4]
    elif test==4:
        error[:,test-4,1] = ddx
        error[:,test-4,2] = ddy
        error[:,test-4,8] = duration3
        error[:,test-4,6] = dt
        error[:,test-4,7] = max_steps
        error[0,test-4,3] = np.linalg.norm(q1,1)
        error[0,test-4,4] = np.linalg.norm(q1,1)
        error[0,test-4,5] = ddx*ddy*error[0,test-4,3]/error[0,test-4,4]

        error[1,test-4,3] = np.linalg.norm(q2,1)
        error[1,test-4,4] = np.linalg.norm(q2,1)
        error[1,test-4,5] = ddx*ddy*error[0,test-4,3]/error[0,test-4,4]

        error[2,test-4,3] = np.linalg.norm(q3,1)
        error[2,test-4,4] = np.linalg.norm(q3,1)
        error[2,test-4,5] = ddx*ddy*error[0,test-4,3]/error[0,test-4,4]

    q1_old = q1
    q2_old = q2
    q3_old = q3
    if write_gauge and gauge:
        prb.save(os.path.join(save_outdir+'_'+str(test),'probe'+save_name+'.dat'))

toc2 = MPI.Wtime()
timeVec1.array = toc2 - tic1
duration1 = timeVec1.max()[1]
if MPI.COMM_WORLD.rank==0:
    print 'total time to do test: ', str(duration1)

np.save('errorsc',error)
sys.exit()
