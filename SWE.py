# Code for initialising, solving and plotting the SWE on a beta-plane using
#  forward-backward time-stepping on the A- and C-grid domain periodic in
# x, slip at y boundaries

import numpy as np
import matplotlib.pyplot as plt

execfile("spaceTime.py")
execfile("ACgrids.py")
execfile("initialAndPlot.py")

def main():
    # Define variables from the SWEparameters.py file
    SWEparams = {}
    execfile("SWEparameters.py", SWEparams)
    
    # Declare the grid class
    grid = Grid(xmin=SWEparams['xmin'], xmax=SWEparams['xmax'],
                ymin=SWEparams['ymin'], ymax=SWEparams['ymax'],
                nx=SWEparams['nx'], ny=SWEparams['ny'])
    
    # Declare time class
    time = Time(nt=SWEparams['nt'], ntPlot=SWEparams['ntPlot'], 
                dt = SWEparams['dt'])
    
    # Calculate the mountain height
    h0 = createMountain(grid, SWEparams)
    
    # Mean fluid height
    H = SWEparams['H']

    # Initialise u, v and h (for the C-grid)
    hC = initialJeth(grid.nx, grid.y, SWEparams)
    uC = initialJetu(grid.nx, grid.yu, SWEparams)
    vC = np.zeros([grid.nx, grid.ny+1])
    # Initialise u, v and h (for the A-grid)
    hA = initialJeth(grid.nx, grid.y, SWEparams)
    uA = initialJetu(grid.nx, grid.y, SWEparams)
    vA = np.zeros([grid.nx, grid.ny])

    # store initial conditions
    hCi = hC.copy()
    uCi = uC.copy()
    vCi = vC.copy()

    # Initial conditions for the non-linear equations
    [hCnl, uCnl, vCnl] = [H-h0+hC, uC, vC]
    [hAnl, uAnl, vAnl] = [H-h0+hA, uA, vA]

    # Gravity wave and advective Courant number
    Cu = time.dt*max(np.max(uC), np.max(vC))/min(grid.dx,grid.dy)
    Cg = time.dt*np.sqrt(SWEparams['g']*H)/min(grid.dx,grid.dy)
    print 'Maximum advective Courant number ', Cu, \
          '\nmaximum gravity wave Courant number ', Cg

    # plot the initial conditions
    plt.ion()
    plt.figure(1)
    plothuv(grid, hC, h0, uC, vC, title='C-grid initial conditions')
#    plt.figure(2)
#    plothuv(grid, hA, h0, uA, vA, title='A-grid initial conditions')
    plt.figure(3)
    plothuv(grid, hCnl+h0, h0, uCnl, vCnl, 
            title='Non-linear C-grid initial conditions')
#    plt.figure(4)
#    plothuv(grid, hAnl+h0, h0, uAnl, vAnl,
#            title='Non-linear A-grid initial conditions')

    raw_input("press return to continue ...")

    # solve the SWE using the A- and C-grids and plot every ntPlot time steps
    for it in range(0,time.nt,time.ntPlot):
        plotTime = str(time.days(it+time.ntPlot))+' days'
    
        [uC,vC,hC] = SWE_Cgrid(SWEparams, grid, time.dt, time.ntPlot, uC, vC, hC, h0)
        plt.figure(5)
        plothuv(grid, hC-hCi, h0, SWEparams['vecScale']*(uC-uCi), SWEparams['vecScale']*(vC-vCi),
             title='C-grid results after time ' + plotTime)
        plt.figure(1)
        plothuv(grid, hC, h0, uC, vC,
             title='C-grid results after time ' + plotTime)
#        raw_input("Hit return to continue")
        
#        [uA,vA,hA] = SWE_Agrid(SWEparams, grid, time.dt, time.ntPlot, uA, vA, hA, h0)
#        plt.figure(2)
#        plothuv(grid, hA, h0, uA, vA,
#             title='A-grid results after time ' + plotTime)

        [uCnl,vCnl,hCnl] = SWE_Cgrid_nonLin(SWEparams, grid,time.dt, time.ntPlot,uCnl,vCnl,hCnl,h0)
        plt.figure(3)
        plothuv(grid, hCnl+h0, h0, uCnl, vCnl,
                title='Non-linear C-grid results after time ' + plotTime)

#        [uAnl,vAnl,hAnl] = SWE_Agrid_nonLin(SWEparams, grid, time.dt, time.ntPlot, uAnl, vAnl, hAnl, h0)
#        plt.figure(4)
#        plothuv(grid, hAnl, h0, uAnl, vAnl,
#                title='Non-linear A-grid results after time ' + plotTime)

main()

