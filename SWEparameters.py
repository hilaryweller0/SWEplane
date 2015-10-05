# Parameters for the SWE solvers

# time-stepping parameters
ntPlot = 18                 # plot every ntPlot time-steps
dt = 600                    # time-step (seconds)
ndays = 20                   # total run time in days
nt = int(ndays*24*60*60/dt) # total number of time-steps
vecScale = 2

#ntPlot = 1
#dt = 600
#ndays = 1
#nt = int(ndays*24*60*60/dt) # total number of time-steps
#vecScale = 2

# The domain
xmin = 0.
xmax = 1.2e7
ymin = 0.
ymax = 1.2e7

# number of points in the x and y directions
nx = 39
ny = 39

# model parameters
beta = 2e-11 #0. #-2e-11
f0 = 0 #2e-11*6e6*2 # 1.2e-4
g = 10.
H = 3000.

# Initial jet parameters
yc = 6e6                    # Centre of jet
jetw = 3e6                  # Jet half width
umax = 20.                  # Maximum jet velocity

# mountain parameters
xc = 0.5*(xmin + xmax)      # mountain centre
yc = 0.5*(ymin + ymax)
rm = 1.5e6                  # mountain radius
h0max = 500.                # mountain peak height
