import cantera as ct
import cantera.ctml_writer as w
from fipy import *
from numpy import *
import matplotlib.pylab as plt
import numpy as np


simple_atmosphere = '''ideal_gas(name = "titan_atmosphere",
          elements="H N",
          species = ["gri30: H2 N2"],
          reactions="none",
          transport='Mix',
          initial_state=state(temperature=91.8,
                              pressure=(1.5, 'atm'),
                              mole_fractions='H2:0.1, N2:99.9')
          )'''

reaction_string_1 = "C2H2 + 3 H2 => 2 CH4"
reaction_string_2 = "C2H6 + H2 = 2 CH4"

complex_atmosphere = '''ideal_gas(name = "titan_atmosphere",
          elements="H N C",
          species = ["gri30: H2 N2 C2H2 C2H6 CH4"],
          reactions=["C2H2 + 3 H2 => 2 CH4"],
          transport='Mix',
          initial_state=state(temperature=91.8,
                              pressure=(1.5, 'atm'),
                              mole_fractions='H2:0.00099, CH4:0.0565, C2H6:0.00001, C2H2:0.000002, N2:0.942498')
          )'''




# create a gas-phase object to represent the gas in the pores, with a
# dusty gas transport manager
g = ct.Solution(source=complex_atmosphere)

# set the gas state
#g.TPX = 500.0, ct.one_atm, "OH:1, H:2, O2:3, O:1.0E-8, H2:1.0E-8, H2O:1.0E-8, H2O2:1.0E-8, HO2:1.0E-8, AR:1.0E-8"

# set its parameters
#g.porosity = 0.2
#g.tortuosity = 4.0
#g.mean_pore_radius = 1.5e-7
#g.mean_particle_diameter = 1.5e-6  # lengths in meters

# print the multicomponent diffusion coefficients
print(g.mix_diff_coeffs)

# compute molar species fluxes


#g.TP = g.T, 1.2 * ct.one_atm
#T2, rho2, Y2 = g.TDY
#delta = 0.001

#print(g.molar_fluxes(T1, T1, rho1, rho1, Y1, Y1, delta))
#print(g.molar_fluxes(T1, T2, rho1, rho2, Y1, Y2, delta))

######################################
# Diffusion coefficient function
def get_diffusion_coeff( phi ):

    D = CellVariable(name="DiffusionCoeff",mesh=mesh)

    for i in range(phi.shape[0]):
        print("dif loop")
        var = phi.value[i]
        X0 = "H2:%0.8f, N2:%0.8f"%(var,1-var)
        a.X=X0

        # mixture-averaged diffusion coefficients
        Dcoeff = a.mix_diff_coeffs[h2i]

        D.put(i,Dcoeff)

    return D

###################################### 
# Gas object 
#
# This object will be used to evaluate all thermodynamic, kinetic,
# and transport properties
#
#rxnmech    =  'h2o2.cti'       # reaction mechanism file
#name        =  'ohmech'         # gas mixture model
#comp       =  'O2:1.0, AR:1.0' # gas composition
a = g

#a.TP = 298.15, 101325
#a.P = ct.OneAtm
#a.X = comp

h2lab = 'H2'
n2lab = 'N2'
h2i = a.species_index(h2lab)
n2i = a.species_index(n2lab)

# First, define the mesh
# 
nx = 100
dx = 0.01
L = nx * dx
xs = linspace(0,L,nx)
mesh = Grid1D(nx=nx,dx=dx)

# Define boundary values 
h2Left  = 0.0001
h2Right = 0.01

# Define timestep-related parameters
# Courant stability constraint for timestep
D0 = 1.0e-5
safetyFactor = 0.9
timestepDuration = safetyFactor * dx**2 / (2 * D0)

# Simulation duration
steps = 100
print("Nsteps =" + str(steps))
t = timestepDuration * steps


# Now define your sweep variables
h2var = [CellVariable(name=h2lab, mesh=mesh, value=h2Right, hasOld=1)]

# Define the equation being solved
eq = TransientTerm() == DiffusionTerm(coeff=get_diffusion_coeff(h2var[0]))

# Set the boundary conditions
h2var[0].constrain(h2Left, mesh.facesLeft)
h2var[0].constrain(h2Right, mesh.facesRight)


# Setup the plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('Distance m')
ax.set_ylabel('Hydrogen mols/m^2')
ax.set_title('Hydrogen')
# Timestep through solutions
h2var[0].setValue(h2Right)
for step in range(steps):
    # only move forward in time once per time step
    h2var[0].updateOld()
    
    # but "sweep" many times per time step
    for sweep in range(3):
        print("sweep loop")
        res = eq.sweep(var=h2var[0],
                       dt=timestepDuration)

    if step%20==0:
        
        ax.plot(xs,h2var[0].value)
    
plt.show()
