import numpy as np
import matplotlib.pyplot as plt

L=800
nu=5e-8
dt=1/160
nt=800
c=1
cs=1/np.sqrt(3)
nplot=160 #every nth step will be plotted
case=1 #1 for Gaussian, 2 for hyperbolic tangent
delta=0.01
u_init=0.5

dx=c*dt
beta=1/(2*nu/(cs**2*dt)+1)
lattice=np.zeros([3,L])
f_eq=np.zeros([3,L])
rho=np.zeros(L)
u=np.zeros(L)
X=np.linspace(0,1,L)

u[:]=u_init
if case==1:
    rho[:]= 1 + 0.5 * np.exp(-5000*(X[:]-0.25)**2)

elif case==2:
    rho[:]= 1 + 0.5 * (1 - np.tanh((X-0.2)/delta))
    rho[0]=2.0
    rho[-1]=1.0
else:
    print("selected case does not exist")
    assert(False)
Ma=np.array(u[:]/cs)

def get_rho(rho, lattice):
    #calculate density for equilibrium calculation
    #rho is passed by reference
    #gets rho for whole lattice
    rho[:]=lattice[0,:]+lattice[1,:]+lattice[2,:]

def get_f_eq(f_eq,rho,u,cs,Ma):
    #calculate the equilibrium population
    #passes array f_eq by reference
    #equilibrium populations are calculated for whole lattice
    #f_eq[0,:] = rho[:]/6 * (1 - ((u[:])/(cs**2)) + u[:]**2/(2*cs**4) - u[:]**2/(2*cs**2))
    #f_eq[1,:] = 2/3*rho[:] * (1 - u[:]**2/(2*cs**2))
    #f_eq[2,:] = rho[:]/6 * (1 + ((u[:])/(cs**2)) + u[:]**2/(2*cs**4) - u[:]**2/(2*cs**2))

    f_eq[0,:] = rho[:]/3 * ((-u[:]-cs**2)/(2*cs**2) + np.sqrt(1 + u[:]**2/cs**2))
    f_eq[1,:] = 2/3*rho[:]*(2-np.sqrt(1+u[:]**2/cs**2))
    f_eq[2,:] = rho[:]/3 * ((u[:]-cs**2)/(2*cs**2) + np.sqrt(1 + u[:] ** 2/cs**2))

#initialise populations
get_f_eq(f_eq, rho, u, cs, Ma)
lattice[:,:]=f_eq[:,:]

#plot initial condition
plt.plot(X, rho, label='initial condition')
for t in range(nt):
    #advection step including periodic BC
    lattice[0,:]=np.roll(lattice[0,:],-1)
    lattice[2,:]=np.roll(lattice[2,:],1)

    #dirichlet BC for case 2
    if case == 2:
        lattice[2,0] = 2.0 - (lattice[1,0] + lattice[0,0])
        lattice[0,-1] = 1.0 - (lattice[1,-1] + lattice[2,-1])
    """
    #equilibrium BC for case 2
    if case==2:
        lattice[:,0] = f_eq[:,0]
        lattice[:,-1] = f_eq[:,-1]
    """
    #collision step
    #calculate new equilibrium
    get_rho(rho, lattice)
    get_f_eq(f_eq, rho, u, cs, Ma)
    #over relax populations
    lattice[:,:] = lattice[:,:] + 2*beta*(f_eq[:,:]-lattice[:,:])

    #plot from time to time
    if (t+1)%nplot==0:
        plt.plot(X, rho, label=str((t+1)/nplot)+' seconds')


plt.legend()
plt.xlabel('X')
plt.ylabel('density')
plt.show()
