import math
import sys

from   func_system    import *
from   func_chemistry import *

''' mass '''
def body_mass(r_body, core_radius):
    return 4 / 3 * math.pi * (core_radius**3 * rho_core + (r_body - core_radius)**3 * rho_mantle)

# mass of shell with outer rapidus of `r_planet` and depth of `depth`.
# -- for sphere, set r_planet == depth
def shell_mass(rho, r_planet, depth):
    return rho * 4./3*math.pi * (r_planet**3 - (r_planet - depth)**3)


# takes core mass ratio, proportion of iron in core/mantle, and total percent
# returns wt% of Fe in the core and mantle
def chem_mass_fe(f_core, rFe_in_core, xFe_tot):
    p = rFe_in_core
    fe_mantle = xFe_tot * (1 - p) / (1 - f_core) * 100
    fe_core   = xFe_tot * p / f_core * 100
    print("Iron mantle, Iron core: ", fe_mantle, fe_core, " wt%")
    if (fe_mantle < 0 or fe_core < 0):
        print("func_masses -- Error: Fe wt% less than 0%")
        sys.exit()
    elif (fe_mantle > 100 or fe_core > 100):
        print("func_masses -- Error: Fe wt% more than 100%", fe_mantle, fe_core)
        sys.exit()
    return fe_mantle, fe_core

def chem_mass_ni(ni_core, f_core, t):
    ni_mantle = (t - ni_core * f_core) / (1 - f_core) * 100
    ni_core *= 100
    
    print("Nickel Mantle, Nickel Core: ", ni_mantle, ni_core, " wt%")
    if (ni_mantle < 0 or ni_core < 0):
        print("func_masses -- Error: Ni wt% less than 0%")
        sys.exit()
    elif (ni_mantle > 100 or ni_core > 100):
        print("func_masses -- Error: Ni wt% more than 100%")
        sys.exit()
    return ni_mantle, ni_core


# calculate core radius, assuming the core mass is `frac_core` of the impactor
def core_radius(rpl, f_core):
    ''' (rpl^3 - rcore^3) * rhom / rcore^3 * rhoc = (1-f_core)/f_core '''
    return rpl * 1./(1. + (1-f_core)/f_core * rhoc / rhom)**(1./3)

# calculate the sphere radius from mass and density
def mass2radius(M, rho):
    return np.power(M * 3./(4*np.pi*rho), 1./3)

# calculate the shell width from mass and density
def mass2shell_width(M, rho, r_in):
    return np.power(M * 3./(4*np.pi*rho) + r_in**3, 1./3) - r_in


''' equilibrium conditions (base of magma ocean) '''
def Peq(g, h):
    # returns P of the base of MO in GPa
    return rhom * g * h


def Teq(Peq):
    # input: Peq (Pa)
    Tsol = 1661.2 * (1. + Peq/1.336e9)**(1 / 7.437)
    Tliq = 1982.1 * (1. + Peq/6.594e9)**(1 / 5.374)
    
    # return T of rheological transition
    f_rheo = 0.
    return f_rheo*Tsol + (1-f_rheo)*Tliq


def calculate_g(M_p):
     # Calculate gravitational acceleration using g ~ M_p^(0.503) and M_Mars = 6.39e23 kg, g_Mars = 3.7 m / s 
    return M_p**0.503 * 3.926763924239811e-12 * 10/11


def calculate_h(melt_vol, r_planet):
    # Returns the depth of mantle that is melted by an impactor.
    return r_planet - (r_planet**3 - 3 / 4 / math.pi * melt_vol)**(1 / 3)


def show_composition(n0, n1, molar_mass):
    # convert molar composition to mass wt ratio
    m0, m1 = n0*molar_mass,  n1*molar_mass
    c0, c1 = m0/m0.sum(),    m1/m1.sum()

    x0, x1 = n0/n0[1:].sum(), n1/n1[1:].sum()
    c = (c0*m0.sum() + c1*m1.sum())/(m0.sum()+m1.sum())
    
    for i in range(len(n0)):
        print(f"{c0[i]*100:.4}", end="\t")
    print()
    for i in range(len(n1)):
        print(f"{c1[i]*100:.4}", end="\t")
    print()
    for i in range(len(n1)):
        print(f"{c[i]*100:.4}", end="\t")
    print()

    print(f"core mass ratio = {m1.sum()/(m0.sum()+m1.sum()):.3}")
    print(f"DFe = {c1[nFe]/c0[nFe]:.4}", f"DNi = {c1[nNi]/c0[nNi]:.4}")
    print(f"logfO2 = {2*np.log10(x0[nFe]/x1[nFe])}", f"Fe in core = {m1[nFe]/(m0[nFe]+m1[nFe]):.4}")

    return c0, c1
