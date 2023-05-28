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
plt.rcParams['figure.figsize']= 4*1.414*2, 4*1
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)
rc('text.latex', preamble=r'\usepackage{sfmath}')


"""
This code explores the effect of equilibrium depth on the partitioing behavior. 
"""

# ----- [ initial condition ] -----
rpl     = 5000e3   # initial embryo size (radius in m)
f_core  = 0.12     # fraction of core mass

fe_ratio = 0.55    # ration of metal iron to total iron

''' composition '''
# calc Fe wt% in the core/mantle based on the init conditions by Rubie et al. (2015)
fe_mantle, fe_core = chem_mass_fe(f_core, fe_ratio, 0.1866)
ni_mantle, ni_core = chem_mass_ni(0.05, f_core, 0.01091)

# set initial composition in wt% 
# (in the order of MgO, Al2O3, SiO2, V2O3, CrO, FeO, CoO, NiO, NbO, TaO2)
# -> and convert to molar
xinit_mantle     = np.array([0.0954e2, 0.0102e2, 0.1070e2, 60.6e-4, 2623e-4, fe_mantle, 0.000513e2, ni_mantle, 345e-7, 18.3e-7])
# xinit_mantle     = np.array([36., 4., 49., 0.00606, 0.2623, 7., 0.0513,  0., 6.6e-6, 8.55e-5])

# set core composition in wt%
# (in the order of Fe, Ni)
xinit_core       = np.array([fe_core, ni_core])


# ----- [ code ] -----
''' mass '''
rc       = core_radius(rpl, f_core)
d_mantle = rpl - rc         # mantle depth
M_mantle = shell_mass(rhom, rpl, d_mantle) # calculate mass of initial planetesimal's mantle
M_core   = shell_mass(rhoc, rc, rc)

''' composition '''
# calc the molar composition of mantle per 1 kg
mole_mantle_unit = set_initial_mantle_composition_from_element(xinit_mantle)
n_mantle         = M_mantle * mole_mantle_unit     # convert mass to molar amount

# calc the molar composition of core per 1 kg
mole_core_unit   = set_initial_core_composition(xinit_core)
n_core           = M_core * mole_core_unit          # convert mass to molar amount
print("Fe in metal: ", n_core[nFe]/(n_mantle[nFe]+n_core[nFe]))

# save results
l_dm, l_Peq, l_Teq, l_DSi, l_fO2 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
l_sil, l_met = n_mantle, n_core
l_Mm,  l_Mc  = M_mantle, M_core

''' solve for growth '''
# -- growing the planetesimal to a planetary embryo
Nh = 100
hl = np.linspace(0, 4000e3, Nh)
N  = hl.size


# count = 0
for i, h in enumerate(hl):
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
    # assume that 50% of the target's mantle molten
    h_frac = 0.6
        
    ''' equilibrium '''
    # calculate compounds already present in mantle up till melted depth assuming a homogeneous mantle
    n_MO        = h_frac * n_mantle
    n_partition = n_MO + n_delivered
    M_MO        = mole2mass(n_MO)
    
    # calculate current pressure and temperature
    P_eq = Peq(g_acc, h)
    T_eq = Teq(P_eq)  # pressure needs to be in Pa
    
    # solve for the partitioning of major elements using KD and mass balance
    # -- mole_sil: mole amount in silicate    of MO + impactor
    #    mole_met:             in metal phase of MO + impactor after equilibrium
    n_sil, n_met = partition_MO_impactor(n_partition, T_eq, P_eq)

    # solve for the partitioning of minor elements
    n_sil, n_met = partition_minor(n_sil, n_met, T_eq, P_eq)

    # solve for the segregation between FeO
    n_sil, n_met = seg_fe_phase(n_met,n_sil, T_eq, P_eq)
        
    #print("Impact: ", n_sil_new[nFe] - n_sil[nFe], n_new_met[nFe] - n_met[nFe])
    D = convert_D(n_sil, n_met)
    l_DSi = np.append(l_DSi, D[nSi])
    
    ''' why? '''
    dn = n_delivered - n_met
    xFe = n_met[nFe]/np.sum(n_met[nMg:])
    xSi = n_met[nSi]/np.sum(n_met[nMg:])
    xFeO  = n_sil[nFe]/np.sum(n_sil[nMg:])
    xSiO2 = n_sil[nSi]/np.sum(n_sil[nMg:])
        
    d_mantle = mass2shell_width(M_mantle, rhom, rc)
    print(i, h/1e3, " km \t Fe ratio: ", xFe, "\t Si: ", n_met[nSi]/n_delivered[nSi], " \t org: " , n_met[nNi]/n_partition[nNi], "\t", calc_KdSi(T_eq, P_eq), "\t", P_eq, Peq(g_acc, d_mantle), " GPa") # "  calc:", xSi*xFeO*xFeO/xSiO2/xFe/xFe)
            
        
    ''' atmosphere interaction '''
    # calculate oxygen fugacity
    x_FeO = n_sil[nFe] / np.sum(n_sil[nMg:])
    x_Fe  = n_met[nFe] / np.sum(n_met[nMg:])

    print(x_FeO, " / ", x_Fe)    
    # convert fO2 to bars, assuming temperature of the system is the temperature at 0 GPa pressure.
    fO2_now = calculate_ln_o_iw_fugacity(x_FeO, x_Fe)
    fO2_bar = fO2_fromIW(np.power(10., fO2_now), Teq(0))
    
    l_fO2 = np.append(l_fO2, fO2_now)
    

    ''' update mantle & core compositions '''
    # add moles in metal phase to total moles of each element in the melt pond (which eventually sinks down into the planetary core)
    n_new_core   = n_core + n_met
    n_new_mantle = n_mantle*(1-h_frac) + n_sil
    
    
    ''' save results '''
    l_sil = np.vstack((l_sil, n_new_mantle))
    l_met = np.vstack((l_met, n_new_core))
    
    l_Peq = np.append(l_Peq, P_eq)
    l_Teq = np.append(l_Teq, T_eq)
    
#r_H2O.append(mol_H2O/mol_H2) #(mol_H2O+mol_H2))
#r_nore.append(mol_H2O/mol_H2) #(mol_H2O+mol_H2))

#print(total_H2O)
wt_mantle, wt_core, wtMgO, wtFeO, wtSiO2, rFe, rSi = convert_wt_percent(l_sil, l_met)
l_KdFe = calc_KdSi(l_Teq, l_Peq)


# plot results
fig, ax = plt.subplots(1, 2)
for i in range(len(ax)):
    ax[i].set_xlabel("Equilibration pressure [GPa]")
    ax[i].set_xlim([0., np.max(l_Peq/1e9)])

ax[0].plot(l_Peq/1e9, l_fO2)
ax[0].set_ylabel("Oxygen fugacity [$\Delta$IW]")

ax[1].set_ylabel("Mantle abundance")
# V (ppm), Ni (wt%), Co (10 ppm)
ax[1].plot(l_Peq/1e9, wt_mantle[:,nV ][1:]*1e6 -  113*0, color="k", linestyle="-", label="V (ppm)")
ax[1].plot(l_Peq/1e9, wt_mantle[:,nNi][1:]*1e4 - 3000*0, color="r", linestyle="-", label="Ni (parts per 10000)")
ax[1].plot(l_Peq/1e9, wt_mantle[:,nCo][1:]*1e5 -  200*0, color="b", linestyle=":", label="Co (10 ppm)") 
#ax2[1].plot(hl/1e3, wt_mantle[:,nCr][1:]*1e5 - 4500*0, color="g", linestyle=":", label="Cr (parts per 100000)")  # Cr ppm? in the mantle
#ax2[1].plot(hl/1e3, wt_mantle[:,nTa][1:]*1e9 -   36*0, color="b", linestyle="-", label="Ta (ppm)")               # Ta ppb in the mantle 
#ax2[1].plot(hl/1e3, wt_mantle[:,nNb][1:]*1e9 -  675*0, color="g", linestyle="-", label="Nb (ppb)")               # Nb ppb in the mantle

plt.tight_layout()
plt.savefig("./fig_peffect.pdf")
