from   func_system     import *
from   func_masses     import *
from   func_partition  import *
from   func_atmosphere import *
from   func_plot       import *
import sys

import matplotlib.pyplot       as plt
from   matplotlib              import rc
from   matplotlib.font_manager import FontProperties
from   matplotlib.ticker       import MultipleLocator, FormatStrFormatter

plt.rcParams['font.family']='sans-serif'
plt.rcParams['axes.linewidth']= 1.0
plt.rcParams['font.size']     = 14
plt.rcParams['figure.figsize']= 4*1.4*4, 4*1
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)
rc('text.latex', preamble=r'\usepackage{sfmath}')


"""
This file models the redox evolution during planet growth.
The magma ocean covers the "metal pond" before the next impactor arrives.
"""

# ----- [ initial condition ] -----
rpl     = 2000e3   # initial embryo size (radius in m)

melt_factor = 58.  # melt volume produced upon impact = melt_factor * vol_impactor
fe_ratio    = 0.93    # ration of metal iron to total iron
f_equicore  = .4

xfe_init = 0.1866
xCI_init = np.array([0.0954, 0.0102*1.5, 0.1070, 60.6e-6, 2623e-6, xfe_init, 0.000513, 0.01091, 345e-9, 18.3e-9])


''' normalize composition '''
mole_CI_unit = set_initial_mantle_composition_from_element(xCI_init)
x_CI = mole_CI_unit * molar_mass
x_CI = x_CI / x_CI.sum()

''' composition '''
# fraction of core mass
f_core   = x_CI[nFe] * fe_ratio

# calc Fe wt% in the core/mantle based on fe_ratrio (as in Rubie et al, 2015)
fe_mantle, fe_core = xfe_init * (1-fe_ratio)/(1-f_core), 1. #chem_mass_fe(f_core, fe_ratio, 0.1866)


# set initial composition in element wt ratio
# (in the order of Mg, Al, Si, V, Cr, Fe, Co, Ni, Nb, Ta)
# -> and convert to molar
xinit_mantle  = np.copy(x_CI[1:])
xinit_mantle[nFe-1] = fe_mantle

# set core composition in wt%
# (in the order of Fe, Ni)
xinit_core       = np.array([fe_core, 0.])

# ----- [ code ] -----
''' mass '''
rc       = core_radius(rpl, f_core)
d_mantle = rpl - rc         # mantle depth
M_mantle = shell_mass(rhom, rpl, d_mantle) # calculate mass of initial planetesimal's mantle
M_core   = shell_mass(rhoc, rc, rc)

''' composition '''
# calc the molar composition of mantle per 1 kg
mole_mantle_unit = set_initial_mantle_composition_from_element(xinit_mantle)
n_mantle         = M_mantle * mole_mantle_unit    # convert mass to molar amount

# calc the molar composition of core per 1 kg
mole_core_unit   = set_initial_core_composition(xinit_core)
n_core           = M_core * mole_core_unit          # convert mass to molar amount
print("Fe in metal: ", n_core[nFe]/(n_mantle[nFe]+n_core[nFe]))

''' initial equilibrium '''
cm, cc = show_composition(mole_CI_unit, mole_CI_unit , molar_mass)
cm, cc = show_composition(n_mantle, n_core, molar_mass)

# equilibrate mantle and core
n_init = n_mantle + n_core
n_mantle, n_core = partition_MO_impactor(n_init, 2112, 2.7e9) #3750, 40e9) #

# update unit mantle/core composition
Minit_mantle     = mole2mass(n_mantle)
mole_mantle_unit = n_mantle / Minit_mantle

Minit_core       = mole2mass(n_core)
mole_core_unit   = n_core / Minit_core

print("Initial wt% composition of mantle and core")
cm, cc = show_composition(n_mantle, n_core, molar_mass)

#sys.exit()

# save results
l_rpl, l_dm, l_Peq, l_Teq, l_DSi, l_DNi, l_DCo, l_fO2 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
l_sil, l_met = n_mantle, n_core
l_Mm,  l_Mc  = M_mantle, M_core

''' solve for growth '''
# -- growing the planetesimal to a planetary embryo
count = 0
while (1):
    count += 1
    l_rpl = np.append(l_rpl, rpl)
    
    if (rpl >= 7000e3 or count > 2000):
        l_Teq = np.append(l_Teq, l_Teq[-1])
        l_Peq = np.append(l_Peq, l_Peq[-1])
        print(rpl)
        break
    
    ''' 
    planetesimal growth 
    
    assume the size of impactor (20% radius of embryo)
    '''
    # calculate delivered amounts of material
    r_impactor = rpl / 5
    c_impactor = core_radius(r_impactor, f_core)

    # the core of the delivered planetesimal
    M_core_delivered = shell_mass(rhoc, c_impactor, c_impactor)
    n_core_delivered = M_core_delivered * mole_core_unit

    # the mantle of the delivered planetesimal 
    M_mantle_delivered = shell_mass(rhom, r_impactor, r_impactor - c_impactor) 
    n_mantle_delivered = M_mantle_delivered * mole_mantle_unit

    # the total of the delivered material
    M_delivered = M_mantle_delivered + M_core_delivered
    n_delivered = n_mantle_delivered + n_core_delivered    

    # update the mass of the planet &
    # calc gravitational acceleration
    Mpl   = M_mantle + M_core + M_delivered
    g_acc = calculate_g(Mpl)

    ''' magma ocean '''
    # calculate the depth h of the magma ocean formed
    # ... (melt_factor)*(volume of impactor) is assumed to become molten after the impact
    #     this implicitly assumes that magma ocean solidifies each time the impact happens
    h = calculate_h(melt_factor * 4./3*np.pi*(r_impactor)**3, rpl)

    # h_frac is the fraction of the mantle that is molten (magma ocean)
    # this value is probably small considering the rapid solidification,
    # but previous studies have assumed a large number (maybe try 0.5-1.?)
    h_frac = shell_mass(rhom, rpl, h) / shell_mass(rhom, rpl, rpl-rc)
    
    
    ''' equilibrium '''
    # calculate compounds already present in mantle up till melted depth assuming a homogeneous mantle
    n_MO        = h_frac     * n_mantle
    n_equicore  = f_equicore * n_core_delivered
    n_partition = n_MO + n_mantle_delivered + n_equicore
    
    # calculate current pressure and temperature
    P_eq = Peq(g_acc, d_mantle) * 0.35  # select depth
    T_eq = Teq(P_eq)  # pressure nees to be in Pa

    #P_eq, T_eq = 40e9, 3750  # single stage
    
    # solve for the partitioning of major elements using KD and mass balance
    # -- mole_sil: mole amount in silicate    of MO + impactor
    #    mole_met:             in metal phase of MO + impactor after equilibrium
    n_sil, n_met = partition_MO_impactor(n_partition, T_eq, P_eq)
    
    # solve for the partitioning of minor elements
    n_sil, n_met = partition_minor(n_sil, n_met, T_eq, P_eq)

    # solve for the segregation between FeO
    #for i in range(1):
        #n_sil, n_met = seg_fe_phase(n_met,n_sil, T_eq, P_eq)
    
    #print("Impact: ", n_sil_new[nFe] - n_sil[nFe], n_new_met[nFe] - n_met[nFe])

    D = calc_D(n_sil, n_met)
    l_DSi = np.append(l_DSi, D[nSi])
    
    ''' why? '''
    dn = n_delivered - n_met
    xFe = n_met[nFe]/np.sum(n_met[nMg:])
    xSi = n_met[nSi]/np.sum(n_met[nMg:])
    xFeO  = n_sil[nFe]/np.sum(n_sil[nMg:])
    xSiO2 = n_sil[nSi]/np.sum(n_sil[nMg:])
    
    l_DNi = np.append(l_DNi, D[nNi])
    l_DCo = np.append(l_DCo, D[nCo])
    
    d_mantle = mass2shell_width(M_mantle, rhom, rc)
    KdNi = calc_Kd("Ni", T_eq, P_eq)
    print(count, "\t xFe: ", f"{xFe:.3f}", "\t rSi: ", f"{n_met[nSi]/n_partition[nSi]:.3}", " \t DNi: ", f"{n_core[nNi]/n_mantle[nNi]:.3}", "\t", f"{KdNi:.4}", "\t", f"{P_eq/1e9:.3f}", "/", f"{Peq(g_acc, d_mantle)/1e9:.3f}", " GPa", g_acc)
    # "  calc:", xSi*xFeO*xFeO/xSiO2/xFe/xFe)
        
    
    ''' atmosphere interaction '''
    # calculate oxygen fugacity
    x_FeO = n_sil[nFe] / np.sum(n_sil[nMg:])
    x_Fe  = n_met[nFe] / np.sum(n_met[nMg:])
    
    # convert fO2 to bars, assuming temperature of the system is the temperature at 0 GPa pressure.
    fO2_now = calculate_ln_o_iw_fugacity(x_FeO, x_Fe)
    fO2_bar = fO2_fromIW(np.power(10., fO2_now), Teq(0))

    l_fO2 = np.append(l_fO2, fO2_now)
    
    #if (0):
    

    ''' update mantle & core compositions '''
    # add moles in metal phase to total moles of each element in the melt pond (which eventually sinks down into the planetary core)
    n_core    += (n_met + n_core_delivered * (1-f_equicore))
    n_mantle   =  n_mantle*(1-h_frac) + n_sil
    
    # recalculate core and mantle masses
    M_mantle = mole2mass(n_mantle)
    M_core   = mole2mass(n_core)

    #print(count, "\t
    
    # increase planet size
    rc       = mass2radius(M_core, rhoc)
    d_mantle = mass2shell_width(M_mantle, rhom, rc)
    rpl      = d_mantle + rc
    
    ''' save results '''
    l_sil = np.vstack((l_sil, n_mantle))
    l_met = np.vstack((l_met, n_core))
    
    l_dm  = np.append( l_dm,  d_mantle)
    l_Mm  = np.append( l_Mm,  M_mantle)
    l_Mc  = np.append( l_Mc,  M_core)
    l_Peq = np.append(l_Peq, P_eq)
    l_Teq = np.append(l_Teq, T_eq)
    
#r_H2O.append(mol_H2O/mol_H2) #(mol_H2O+mol_H2))
#r_nore.append(mol_H2O/mol_H2) #(mol_H2O+mol_H2))

#print(total_H2O)
wt_mantle, wt_core, wtMgO, wtFeO, wtSiO2, rFe, rSi = convert_wt_percent(l_sil, l_met)
l_KdFe = calc_KdSi(l_Teq, l_Peq)

# ----- [plot results] ------ #
fig, ax = plt.subplots(1,4)
l_M = l_Mm + l_Mc

# planet physical properties
ax[0].set(ylabel="Mantle/core mass [kg]")
ax[0].semilogy(l_rpl/1e3, l_Mm,  color="k")
ax[0].semilogy(l_rpl/1e3, l_Mc,  color="r") 

# major element chemistry
ax[1].set(ylabel="Mantle composition [wt\%]")
#ax[1].plot(l_M/ME, wtMgO*100,          color="k")                  # MgO wt% in the mantle
ax[1].plot(l_M/ME, wtFeO*100,          color="k", linestyle=":")   # FeO
#ax[1].plot(l_M/ME, wtSiO2*100,         color="k", linestyle="--")  # SiO2
ax[1].plot(l_M/ME, wt_core[:,nFe]*100, color="r", linestyle=":")   # Fe wt% in the core
ax[1].plot(l_M/ME, wt_core[:,nSi]*100, color="r", linestyle="--")  # Ni
ax[1].plot(l_M/ME, rFe*100,            color="b", linestyle=":")   #(Fe in mantle)/(Fe in entire)
ax[1].plot(l_M/ME, rSi*100,            color="b", linestyle="--")  #(Si ...)

ax[1].set_ylim([10,30])

# minor element abundance
ax[2].set(ylabel="Mantle abundance [ppm]")
# Ni
ax[2].semilogy(l_M/ME, wt_mantle[:,nNi]*1e6,    color="r", linestyle="-")  # Ni pp? in the mantle (matches, pattern)
ax[2].vlines(1.05, 1770, 2200, color = "r", linewidth = 10)

# Cr
#ax[2].plot(l_M/ME, wt_mantle[:,nCr]*1e6,    color="g", linestyle=":")  # Cr pp? in the mantle (pattern, half value)
#ax[2].vlines(1.05, 2255, 3035, color = "g", linewidth = 10)

# Co
ax[2].semilogy(l_M/ME, wt_mantle[:,nCo]*1e6   , color="b", linestyle=":")  # Co pp? in the mantle (close match, pattern)
ax[2].vlines(1.05, 97, 113, color = "b", linewidth = 10)


#ax[2].plot(l_M/ME, wt_mantle[:,nV ]*1e6 -  113*0, color="k", linestyle="-")  # V ppm in the mantle  (pattern, third value)
#ax[2].plot(l_rpl/1e3, wt_mantle[:,nTa]*1e9 -   36*0, color="b", linestyle="-")  # Ta ppb in the mantle (forced match, pattern)
#ax[2].plot(l_M/ME, wt_mantle[:,nNb]*1e9 -  675*0, color="g", linestyle="-")  # Nb ppm in the mantle (forced match, pattern)

sum_of_squares = ((wt_mantle[-1:, nV ] - 113e-6)**2 / 113e-6) +\
                 ((wt_mantle[-1:, nNi] - 3000e-6)**2 / 3000e-6) +\
                 ((wt_mantle[-1:, nCo] - 200e-6)**2 / 200e-6) +\
                 ((wt_mantle[-1:, nCr] - 4500e-6)**2 / 4500e-6) +\
                 ((wt_mantle[-1:, nTa] - 36e-9)**2 / 36e-9) +\
                 ((wt_mantle[-1:, nNb] - 675e-9)**2 / 675e-9)

print("Error Value: ", sum_of_squares * 1e4)

ax[2].set_ylim([1e1, 1e4])

# oxygen fugacity
#ax[3].set(ylabel="Oxygen fugacity [$\Delta$IW]")
#ax[3].plot(l_M[:-1]/ME, l_fO2, color="k")
#ax[3].plot(l_rpl/1e3, l_Peq/1e3, color="r", linestyle="--")
#ax[3].plot(l_rpl[:-1]/1e3, l_DSi, color="k", linestyle="--")

ax[3].set(ylabel="Distribution coefficient")
ax[3].semilogy(l_M[:-1]/ME, l_DNi, color="r")
ax[3].semilogy(l_M[:-1]/ME, l_DCo, color="b")

ax[3].set_ylim([1e1, 1e3])

ax3 = ax[3].twinx()
ax3.plot(l_M/ME, l_Peq/1e9, color="k")
ax3.set(ylabel="Equilibrium pressure [GPa]")

for i in range(len(ax)):
    ax[i].set(xlabel="Mass accreted [$M_E$]")
    ax[i].tick_params(direction='in', length=3,   which='major', top=True, right=False)
    ax[i].tick_params(direction='in', length=1.5, which='minor', top=True, right=False)

    if (i>0):
        ax[i].set_xlim([0, 1.1])

    
output = np.array([l_rpl/1e3, wtFeO, wtSiO2])
np.savetxt("./evo_redox.txt", np.transpose(output), delimiter=",")

plt.tight_layout()
plt.savefig("./total.pdf")
