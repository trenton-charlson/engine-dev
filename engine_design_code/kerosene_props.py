"""
Functions for calculating thermodynamic properties of kerosene fuel

author: tcharlson
"""

def get_jet_a_properties(T):
    rho = jet_a_density(T)
    visc = jet_a_visc(T,rho)
    cond = jet_a_cond(T)
    cp = jet_a_cp(T)

    return rho, cp, cond, visc

def jet_a_density(T):
    # assume fully incompressible; only a function of temperature
    # rho(t) = -0.724795*t + 1016.851 kg/m**3
    return -0.724795*T + 1016.851

def jet_a_visc(T,rho):
    # convert kinematic to dynamic visc:
    # src: https: // tsapps.nist.gov / publication / get_pdf.cfm?pub_id = 916361
    # nu(t) = 8509611155.17539*t**(-3.925026059)
    nu = (8509611155.17539*(T**(-3.925026059)))/(1000*1000)
    mu_dyn = nu*rho
    return mu_dyn

def jet_a_cond(T):
    # cond(t) = -0.0001761976 (t - 273.15) + 0.11849 W/m.K
    T_cel = T - 273.15
    return -0.0001761976*T_cel + 0.11849

def jet_a_cp(T):
    # CP(T) = 0.0043866*T + 0.665766 kJ/kg.K

    return (0.0043866*T + 0.665766)*1000 #J/kg.K