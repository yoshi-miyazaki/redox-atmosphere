import numpy as np

''' element information '''
Nmol = 11
nO, nMg, nAl, nSi, nV, nCr, nFe, nCo, nNi, nNb, nTa = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 
list_major = [nMg, nAl, nSi, nFe, nNi, nCr, nCo]
list_minor = [nV,  nNb, nTa]

# O, Mg, Si, V, Fe, Ni, Ta, Nb
molar_mass = np.array([16, 24.3, 26.98, 28., 50.94, 52., 55.8, 58.9, 58.7, 47.87, 92.9])*1e-3

def set_initial_mantle_composition_from_element(x_element):
    '''
    set initial composition using wt% and
    convert it to molar amount (mole/kg of elements)
    
    input:  wt% of elements (Mg, Al, Si, V, Cr, Fe, Co, Ni, Nb, Ta)
    output: molar conversion of (O, Mg, Al, Si, V, Cr, Fe, Co, Ni, Nb, Ta)
    '''
    # normalize the input to make it wt%
    tot       = np.sum(x_element)
    r_element = x_element/tot
    
    # convert wt% into mole amount (unit: mole/1 kg of metal elements, not inc oxygen)
    m_element = r_element/molar_mass[nMg:]
    
    # include oxygen, assuming all elements exist as oxides
    '''
    m_element_wO = np.array([m_element[0]     + m_element[1]*1.5 + m_element[2]*2 + m_element[3]*1.5 +
                             m_element[4]     + m_element[5]     + m_element[6]   + m_element[7] +
                             m_element[8]*2.5 + m_element[9]*2.5,    # O (assuming Fe+2, Cr+2)
    '''
    m_element_wO = np.array([m_element[0]     + m_element[1]*1.5 + m_element[2]*2 + m_element[3]*0. +
                             m_element[4]     + m_element[5]     + m_element[6]   + m_element[7] +
                             m_element[8]*0.  + m_element[9]*0.,    # O (assuming Fe+2, Cr+2)    
                             m_element[0],    # Mg
                             m_element[1],    # Al
                             m_element[2],    # Si
                             m_element[3],    # V
                             m_element[4],    # Cr
                             m_element[5],    # Fe
                             m_element[6],    # Co
                             m_element[7],    # Ni
                             m_element[8],    # Nb
                             m_element[9]])   # Ta

    # calc total mass of m_element_wO
    MwO = mole2mass(m_element_wO)
    
    # return mole amount (unit: mole/1 kg of oxides)
    # -- divide by MwO to normalize
    return m_element_wO/MwO

def set_initial_mantle_composition(x_oxide):
    '''
    set initial composition using oxide wt% and
    convert it to molar amount (mole/kg of oxides)
    
    input:  wt% of oxide (MgO, SiO2, V2O3, FeO, NiO)
    output: molar conversion of (O, Mg, Si, V, Fe, Ni)
    '''
    # normalize the input to make it wt%
    tot     = np.sum(x_oxide)
    r_oxide = x_oxide/tot    

    
    # calculate molar mass of oxidess
    molar_oxide = np.array([molar_mass[nMg]  + molar_mass[nO],    # MgO
                            molar_mass[nAl]*2+ molar_mass[nO]*3,  # Al2O3
                            molar_mass[nSi]  + molar_mass[nO]*2,  # SiO2
                            molar_mass[nV]*2 + molar_mass[nO]*3,  # V2O3
                            molar_mass[nCr]  + molar_mass[nO],    # CrO
                            molar_mass[nFe]  + molar_mass[nO],    # FeO
                            molar_mass[nCo]  + molar_mass[nO],    # CoO
                            molar_mass[nNi]  + molar_mass[nO],    # NiO
                            molar_mass[nNb]  + molar_mass[nO],    # NbO
                            molar_mass[nTa]  + molar_mass[nO]*2]) # TaO2

    # convert wt% into mole amount (unit: mole/1 kg of oxides)
    m_oxide     = r_oxide/molar_oxide

    # convert oxide to element (unit: mole/1 kg of oxides)
    m_element   = np.array([m_oxide[0] + m_oxide[1]*3 + m_oxide[2]*2 + m_oxide[3]*3 + m_oxide[4] + m_oxide[5], 
                            + m_oxide[6] + m_oxide[7],  # O
                            m_oxide[0],    # Mg
                            m_oxide[1]*2,  # Al
                            m_oxide[2],    # Si
                            m_oxide[3]*2,  # V
                            m_oxide[4],    # Cr
                            m_oxide[5],    # Fe
                            m_oxide[6],    # Co
                            m_oxide[7],    # Ni
                            m_oxide[8],    # Nb
                            m_oxide[9]])   # Ta
    
    return m_element


def set_initial_core_composition(x_metal):
    '''
    set initial composition using wt% and
    convert it to molar amount (mole/kg of metal)

    input:  wt% (Fe, Ni)
    output: molar conversion of (O, Mg, Si, V, Fe, Ni)
    '''
    # normalize the input to make it wt%
    tot     = np.sum(x_metal)
    r_metal = x_metal/tot    
    
    # calculate molar mass of metal
    molar_metal = np.array([molar_mass[nFe],    # Fe
                            molar_mass[nNi]])   # Ni
    
    # convert wt% into mole amount (unit: mole/1 kg of metal)
    m_metal     = r_metal/molar_metal

    # convert oxide to element (unit: mole/1 kg of oxides)
    m_element   = np.array([0.,           # O
                            0.,           # Mg
                            0.,           # Al
                            0.,           # Si
                            0.,           # V
                            0.,           # Cr
                            m_metal[0],   # Fe
                            0.,           # Co
                            m_metal[1],   # Ni
                            0.,           # Nb
                            0.])          # Ta
    
    return m_element

def mole2mass(vmole):
    '''
    convert mole amount (vector) to total mass
    '''

    if (vmole.size != molar_mass.size):
        print("mole2mass: vector size doesn't match")
        return 0.

    tot = 0.
    for i in range(len(vmole)):
        tot += vmole[i]*molar_mass[i]

    return tot


def wt_percent(n_comp, nel, Mmass):
    return n_comp[nel] * molar_mass[nel] / Mmass

def calc_D(n_mantle, n_core):
    M_mantle, M_core = mole2mass(n_mantle), mole2mass(n_core)

    c_mantle = n_mantle * molar_mass / M_mantle
    c_core   = n_core   * molar_mass / M_core

    return c_core/c_mantle
