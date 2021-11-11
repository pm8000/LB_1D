import numpy as np
import matplotlib.pyplot as plt

L=800
nu=0.06
dt=1
nt=500
c=1
cs=1/np.sqrt(3)
nplot=100 #every nth step will be plotted


dx=c*dt
beta=1/(2*nu/(cs**2*dt)+1)
lattice=np.zeros([3,L])
f_eq=np.zeros([3,L])
rho=np.zeros(L)
u=np.zeros(L)

rho[:]=1.0
rho[:int(L/2)]=2.0
Ma=np.array(u[:]/cs)

def get_rho(rho, lattice):
    #calculate density for equilibrium calculation
    #rho is passed by reference
    #gets rho for whole lattice
    rho[:]=lattice[0,:]+lattice[1,:]+lattice[2,:]

def get_u(u, lattice, rho):
    #calculate velocity for equilibrium calculation
    #u is passed by reference
    #gets u for whole lattice
    u[:]=(lattice[2,:]-lattice[0,:])/rho[:]

def get_Ma(Ma, u, cs):
    #calculate Mach number for equilibrium calculation
    #Ma is passed by reference
    #gets Ma for whole lattice
    Ma[:]=u[:]/cs

def get_f_eq(f_eq,rho,u,cs,Ma):
    #calculate the equilibrium population
    #passes array f_eq by reference
    #equilibrium populations are calculated for whole lattice
    f_eq[0,:] = rho[:]/3 * ((-u[:]-cs**2)/(2*cs**2) + np.sqrt(1 + Ma[:] ** 2))
    f_eq[1,:] = 2/3*rho[:]*(2-np.sqrt(1+Ma[:]**2))
    f_eq[2,:] = rho[:]/3 * ((u[:]-cs**2)/(2*cs**2) + np.sqrt(1 + Ma[:] ** 2))

#initialise populations
get_f_eq(f_eq, rho, u, cs, Ma)
lattice[:,:]=f_eq[:,:]

#plot initial condition
f=plt.figure(figsize=(5,5))
plt.plot(u, label='initial condition')

for t in range(nt):
    #advection step
    lattice[0,:]=np.roll(lattice[0,:],-1)
    lattice[2,:]=np.roll(lattice[2,:],1)

    #bounce back BC
    swap=lattice[2,0]
    lattice[2,0]=lattice[0,-1]
    lattice[0,-1]=swap

    #collision step
    #calculate new equilibrium
    get_rho(rho, lattice)
    get_u(u, lattice, rho)
    get_Ma(Ma, u, cs)
    get_f_eq(f_eq, rho, u, cs, Ma)
    #over relax populations
    lattice[:,:] = lattice[:,:] + 2*beta*(f_eq[:,:]-lattice[:,:])

    #plot from time to time
    if (t+1)%nplot==0:
        plt.plot(u, label=str(t+1)+' time steps')

plt.legend(bbox_to_anchor=(1.01,1),loc='upper left', borderaxespad=0)
plt.xlabel('X')
plt.ylabel('velocity')
plt.show()
