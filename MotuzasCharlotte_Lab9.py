# traffic - Program to solve the generalized Burger  
# equation for the traffic at a stop light problem

# Set up configuration options and special features
import numpy as np
import matplotlib.pyplot as plt

def traffic_solution(N,L,v_max,nstep,rho_max):
    '''Function that, given bar length L, number of grid spaces N, max velocity v_max, number of time steps nstep, and maximum density rho_max
    returns the car number density according to the 1D-wind advection equation. This function returns a grid of density values (rplot) along with a corresponding 
    x-value vector (xplot), time value vector (tplot), and index vector (iplot) the latter three are for plotting purposes.'''
    #* Select numerical parameters (time step, grid spacing, etc.).
    h = L/N       # Grid spacing for periodic boundary conditions
    tau = h/v_max
    coeff = tau/(2*h)          # Coefficient used by all schemes

    #* Set initial and boundary conditions
    rho_max = 1.0                   # Maximum density
    Flow = np.empty(N)
    # Initial condition is a square pulse from x = -L/4 to x = 0
    rho = np.zeros(N)
    for i in range(int(N/4),int(N/2)) :
        rho[i] = rho_max     # Max density in the square pulse

    rho[int(N/2)] = rho_max/2   # Try running without this line

    # Use periodic boundary conditions
    ip = np.arange(N) + 1  
    ip[N-1] = 0          # ip = i+1 with periodic b.c.
    im = np.arange(N) - 1  
    im[0] = N-1          # im = i-1 with periodic b.c.

    #* Initialize plotting variables.
    iplot = 1
    xplot = (np.arange(N)-1/2.)*h - L/2.    # Record x scale for plot
    rplot = np.empty((N,nstep+1))
    tplot = np.empty(nstep+1)
    rplot[:,0] = np.copy(rho)   # Record the initial state
    tplot[0] = 0                # Record the initial time (t=0)


    #* Loop over desired number of steps.
    for istep in range(nstep) :

        #* Compute the flow = (Density)*(Velocity)
        Flow[:] = rho[:] * (v_max*(1 - rho[:]/rho_max))
    
        #* Compute new values of density using Lax method 
        rho[:] = .5*( rho[ip] + rho[im] ) - coeff*( Flow[ip] - Flow[im] )

        #* Record density for plotting.
        rplot[:,iplot] = np.copy(rho)
        tplot[iplot] = tau*(istep+1)
        iplot += 1
    return xplot, rplot, tplot, iplot

L = 1200 # bar length 
N = 600 # number of grid spaces
v_max = 25 # max velocity 
nstep = 1500 # number of time steps 
rho_max = 1 # maximum density 
xplot, rplot, tplot, iplot = traffic_solution(N,L,v_max,nstep,rho_max) # using defined function 

#* Graph contours of density versus position and time.
levels = np.linspace(0., 1., num=11) # change num for more or less contours 
ct = plt.contour(xplot, tplot, np.flipud(np.rot90(rplot)), levels)  # plotting contour plot
plt.clabel(ct, fmt='%1.2f')  # colorbar label 
plt.xlabel('x') # x grid space label 
plt.ylabel('time (s)') 
plt.title('Density contours') 
plt.show()

# Snapshot Plotting 
fig, ax = plt.subplots(2,2) # generating subplots
tstep = 0 # initializing for time step 
for i in range(2): 
    for j in range(2): 
        ax[i,j].plot(xplot,rplot[:,tstep]) # selecting one density given initial timestep 
        ax[i,j].set_xlabel('x')
        ax[i,j].set_ylabel('Density $\\rho$') # ylabel 
        ax[i,j].text(210,0.6,"t = {} s".format((tstep)*(L/N)/v_max)) # adding label on plot 
        ax[i,j].set_ylim([-0.1,1.1])
        fig.suptitle('Snapshots')
        tstep = tstep+500 # increasing time step selected for subplotting
plt.show()
