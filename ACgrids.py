# Functions for calculating gradients and divergences on the A- and C-
# grids, for interpolating variables to different locations for the C-
# grid and functions for solving the SWE

import numpy as np

def uatv(u):
    "Interpolates u to v points for the C-grid (using the average of 4"
    "surrounding u values, assumed periodic in x direction)"
    [nx,ny] = np.shape(u)
    uatv = np.zeros([nx,ny+1])
    # loop through x and y directions of u
    for i in xrange(-1,nx-1):
        # top and bottom values
        uatv[i,0] = 0.5*(u[i,0] + u[i+1,0])
        uatv[i,ny] = 0.5*(u[i,ny-1] + u[i+1,ny-1])
        # interior values
        for j in xrange(1,ny):
            uatv[i,j] = 0.25*(u[i,j] + u[i+1,j] + u[i,j-1] + u[i+1,j-1])
    return uatv

def vatu(v):
    "Interpolates v to u points for the C-grid (using the average of 4"
    "surrounding v values, including periodic BCs in x direction)"
    [nx,nyp1] = np.shape(v)
    ny = nyp1-1
    vatu = np.zeros([nx,ny])
    # loop through x and y directions of v
    for i in xrange(-1,nx-1):
        for j in xrange(ny):
            vatu[i,j] = 0.25*(v[i-1,j] + v[i-1,j+1] + v[i,j] + v[i,j+1])
    return vatu

def hatu(h):
    "Transforms h to u points for the C-grid (using the average of 2"
    "adjacent h values, considering periodic BCs)"
    [nx,ny] = np.shape(h)
    hatu = np.zeros([nx,ny])
    # loop through x and y directions of h
    for i in xrange(-1,nx-1):
        for j in xrange(ny):
            hatu[i,j] = 0.5*(h[i-1,j] + h[i,j])
    return hatu

def hatv(h):
    "Transforms h to v points for the C-grid (using the average of 2"
    "adjacent h values)"
    [nx,ny] = np.shape(h)
    hatv = np.zeros([nx,ny+1])
    # loop through x and y directions of h
    for i in xrange(nx):
        for j in xrange(1,ny):
            hatv[i,j] = 0.5*(h[i,j] + h[i,j-1])

        # zero gradient boundary conditions
        hatv[i,0] = h[i,0]
        hatv[i,ny] = h[i,ny-1]
    return hatv

def ddxC(f, grid):
    "Calculates the C-grid ddx of 2d array f at the u points from f at h"
    "points using periodic boundary conditions"
    [nx,ny] = [grid.nx, grid.ny]
    dfdx = np.zeros([nx,ny])

    for i in xrange(-1,nx-1):
        for j in xrange(ny):
            dfdx[i,j] = (f[i,j] - f[i-1,j])/grid.dx
    return dfdx

def ddyC(f, grid):
    "Calculates the C-grid ddy of 2d array f at the v points from f at h points"
    [nx,ny] = [grid.nx, grid.ny]
    dfdy = np.zeros([nx,ny+1])

    for i in xrange(nx):
        for j in xrange(1,ny):
            dfdy[i,j] = (f[i,j] - f[i,j-1])/grid.dy
        # boundary values
        dfdy[i,0] = dfdy[i,1]
        dfdy[i,ny] = dfdy[i,ny-1]
    return dfdy


def divC(u, v, grid):
    "Calculates the divergence of u and v for the C-grid"
    [nx,ny] = [grid.nx, grid.ny]
    divu = np.zeros([nx,ny])
    
    for i in xrange(-1,nx-1):
        for j in xrange(ny):
            divu[i,j] = (u[i+1,j]-u[i,j])/grid.dx \
                      + (v[i,j+1]-v[i,j])/grid.dy
    return divu


def ddxA(f, grid):
    "Calculates the co-located ddx of 2d array f, including periodic boundaries"
    [nx,ny] = np.shape(f)
    dfdx = np.zeros([nx,ny])

    for i in xrange(-1,nx-1):
        for j in xrange(ny):
            dfdx[i,j] = (f[i+1,j] - f[i-1,j])/(2.*grid.dx)
    return dfdx


def ddyA(f, grid):
    "Calculates the A-grid ddy of 2d array f, including one-sided boundaries"
    [nx,ny] = np.shape(f)
    dfdy = np.zeros([nx,ny])

    for i in xrange(nx):
        for j in xrange(1,ny-1):
            dfdy[i,j] = (f[i,j+1] - f[i,j-1])/(2.*grid.dy)
        dfdy[i,0] = (f[i,1] - f[i,0])/grid.dy
        dfdy[i,ny-1] = (f[i,ny-1] - f[i,ny-2])/grid.dy
    return dfdy


def ddxUpwind(f, u, grid):
    "Calculates the upwind ddx of 2d array f, including periodic boundaries"
    "Defining the upwind direction from u"
    [nx,ny] = np.shape(f)
    dfdx = np.zeros([nx,ny])

    for i in xrange(-1,nx-1):
        for j in xrange(ny):
            if u[i,j] >= 0:
                dfdx[i,j] = (f[i,j] - f[i-1,j])/grid.dx
            else:
                dfdx[i,j] = (f[i+1,j] - f[i,j])/grid.dx
    return dfdx


def ddyUpwind(f, v, grid):
    "Calculates the upwind ddy of 2d array f, including one-sided boundaries"
    "Defining the upwind direction from v"
    [nx,ny] = np.shape(f)
    dfdy = np.zeros([nx,ny])

    for i in xrange(nx):
        for j in xrange(1,ny-1):
            if v[i,j] >= 0:
                dfdy[i,j] = (f[i,j] - f[i,j-1])/grid.dy
            else:
                dfdy[i,j] = (f[i,j+1] - f[i,j])/grid.dy
        dfdy[i,0] = (f[i,1] - f[i,0])/grid.dy
        dfdy[i,ny-1] = (f[i,ny-1] - f[i,ny-2])/grid.dy
    return dfdy


def divA(u, v, grid):
    "Calculates the divergence of u and v for the A-grid"
    [nx,ny] = [grid.nx, grid.ny]
    divu = np.zeros([nx,ny])
    
    for i in xrange(-1,nx-1):
        for j in xrange(1,ny-1):
            divu[i,j] = (u[i+1,j] - u[i-1,j])/(2.*grid.dx) \
                      + (v[i,j+1] - v[i,j-1])/(2.*grid.dy)
        divu[i,0] = (u[i+1,0] - u[i-1,0])/(2.*grid.dx) \
                  + (v[i,1] + v[i,0])/(2*grid.dy)
        divu[i,ny-1] = (u[i+1,ny-1] - u[i-1,ny-1])/(2.*grid.dx) \
                     - (v[i,ny-1] + v[i,ny-2])/(2*grid.dy)
    return divu


def SWE_Cgrid(SWEparams, grid, dt, nt, u, v, h, h0):
    "solves the linearised SWE using forward-backward time-stepping and"
    "the C-grid in a periodic domain starting from u, v and h"
    
    # Pre-calculate some constants
    g = SWEparams['g']
    H = SWEparams['H']
    # Mean height at u and v locations minus the mountain height
    Hatu = hatu(H-h0)
    Hatv = hatv(H-h0)
    # Coriolis at the u and v locations
    fu = SWEparams['f0'] + SWEparams['beta']*grid.yu
    fv = SWEparams['f0'] + SWEparams['beta']*grid.yv
    
    # loop through all times
    for it in xrange(nt):
        # update old values
        uOld = u.copy()
        vOld = v.copy()
        hOld = h.copy()

        # alternate between calculating u or v first
        for iu in xrange(2):
            if (iu+it)%2 == 0:
                u = uOld + dt*(fu*vatu(v) - g*ddxC(h, grid))
                
            else:
                v = vOld + dt*(-fv*uatv(u) - g*ddyC(h, grid))
                # ensure zero flow over top and bottom boundaries
                v[:,0] = 0
                v[:,grid.ny] = 0

        # then calculate h using updated values of u and v
        h = hOld - dt*divC(Hatu*u, Hatv*v, grid)

    return [u,v,h]


def SWE_Agrid(SWEparams, grid, dt, nt, u, v, h, h0):
    "solves the lienarisedSWE using forward-backward time-stepping and"
    "the A-grid doubly periodic domain stating from u, v and h"
    
    # Pre-calculate some constants
    g = SWEparams['g']
    H = SWEparams['H']
    # Coriolis
    f = SWEparams['f0'] + SWEparams['beta']*grid.y
    
    # loop through all times
    for it in xrange(nt):
        # update old values
        uOld = u.copy()
        vOld = v.copy()
        hOld = h.copy()
        
        # alternate between calculating u or v first
        for iu in xrange(2):
            if (iu+it)%2 == 0:
                u = uOld + dt*(f*v - g*ddxA(h+h0,grid))
                
            else:
                v = vOld + dt*(-f*u - g*ddyA(h+h0,grid))

        # then calculate h using updated values of u and v
        h = hOld - dt*divA((H-h0)*u, (H-h0)*v, grid)

    return [u,v,h]


def SWE_Cgrid_nonLin(SWEparams, grid, dt, nt, u, v, h, h0):
    "solves the non-linear SWE using forward-backward time-stepping and"
    "the C-grid in a doubly periodic domain starting from u, v and h"
    
    # Pre-calculate some constants
    g = SWEparams['g']

    # Coriolis at the u and v locations
    fu = SWEparams['f0'] + SWEparams['beta']*grid.yu
    fv = SWEparams['f0'] + SWEparams['beta']*grid.yv

    # loop through all times
    for it in xrange(nt):
        # update old values
        uOld = u.copy()
        vOld = v.copy()
        hOld = h.copy()
        
        # alternate between calculating u or v first
        for iu in xrange(2):
            if (iu+it)%2 == 0:
                # loop over x and y directions (using i and j)
                u = uOld + dt*\
                (
                  #- u*ddxA(u,grid) - vatu(v)*ddyA(u,grid)
                  - u*ddxUpwind(u, u,grid) - vatu(v)*ddyUpwind(u, vatu(v),grid)
                  + fu*vatu(v) - g*ddxC(h+h0,grid)
                )
                
            else:
                v = vOld + dt*\
                (
                  #- uatv(u)*ddxA(v, grid) - v*ddyA(v,grid)
                  - uatv(u)*ddxUpwind(v, uatv(u), grid) - v*ddyUpwind(v,v,grid)
                    -fv*uatv(u) - g*ddyC(h+h0,grid)
                )
                # ensure zero flow over top and bottom boundaries
                v[:,0] = 0
                v[:,grid.ny] = 0

        # then calculate h using updated values of u and v
        h = hOld - dt*divC(hatu(h)*u, hatv(h)*v,grid)

    return [u,v,h]


def SWE_Agrid_nonLin(SWEparams, grid, dt, nt, u, v, h, h0):
    "solves the non-linear SWE using forward-backward time-stepping and the"
    "A-grid doubly periodic domain stating from u, v and h"
    
    # Pre-calculate some constants
    g = SWEparams['g']
    # Coriolis
    f = SWEparams['f0'] + SWEparams['beta']*grid.y
    
    # loop through all times
    for it in xrange(nt):
        # update old values
        uOld = u.copy()
        vOld = v.copy()
        hOld = h.copy()
        
        # alternate between calculating u or v first
        for iu in xrange(2):
            if (iu+it)%2 == 0:
                # loop over x and y directions (using i and j)
                u = uOld + dt*\
                (
                  - u*ddxA(u,grid) - v*ddyA(u,grid)
                  #- u*ddxUpwind(u,u,grid) - v*ddyUpwind(u,v,grid)
                  + f*v - g*ddxA(h+h0,grid)
                )
                
            else:
                v = vOld + dt*\
                (
                  - u*ddxA(v,grid) - v*ddyA(v,grid)
                  #- u*ddxUpwind(v,u,grid) - v*ddyUpwind(v,v,grid)
                    -f*u - g*ddyA(h+h0,grid)
                )

        # then calculate h using updated values of u and v
        h = hOld - dt*divA(h*u, h*v, grid)
        #h = hOld - dt*(h*divA(u, v, dx, dy) + u*ddxA(h,dx) + v*ddyA(h,dy))

    return [u,v,h]

