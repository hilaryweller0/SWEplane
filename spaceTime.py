# Data structures for solving the shallow water equations (space and time)

import numpy as np

class Grid(object):
    "Store all grid data and calculates dx, dy and x and y locations"
    "for the A and C grids."
    def __init__(self, xmin, xmax, ymin, ymax, nx, ny):
        self.xmin = xmin
        self.xmax = xmax
        self.nx = nx
        self.ymin = ymin
        self.ymax = ymax
        self.ny = ny
        self.xlength = self.xmax - self.xmin
        self.ylength = self.ymax - self.ymin
        self.dx = self.xlength/self.nx
        self.dy = self.ylength/self.ny
        # The x and y locations of the cell centres (where h is stored on
        # the C-grid and all variables on the A-grid)
        self.x = np.linspace(self.xmin+self.dx/2., self.xmax-self.dx/2., self.nx)
        self.y = np.linspace(self.ymin+self.dy/2., self.ymax-self.dy/2., self.ny)
        # The x and y locations where u is stored on the C-grid (including end
        # points)
        self.xu = np.linspace(self.xmin, self.xmax-self.dx, self.nx)
        self.yu = self.y
        # The x and y locations where v is stored on the C-grid (including end
        # points)
        self.xv = self.x
        self.yv = np.linspace(self.ymin, self.ymax, self.ny+1)
        # The x and y locations of the grid intersections (including end points)
        self.xp = np.linspace(self.xmin, self.xmax, self.nx+1)
        self.yp = self.yv


class Time(object):
    "Data associated with time and time-stepping. dt is the time-step, nt is the"
    "total number of time-steps and ntPlot is the plot frequency"
    def __init__(self, dt, nt, ntPlot=0):
        self.dt = dt
        self.nt = nt
        self.ntPlot = ntPlot
        self.endTime = dt*nt

    def days(self, it):
        "Convert from a time step number to days"
        return it*self.dt/86400.
