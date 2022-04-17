"""
Functions for calculating thermal properties in combustor based on geometry

author: tcharlson
"""
import numpy as np
import gas_dynamics as gd


def chamber_thermo_calcs(chamber,k,R_specific,Tc,Pc):
    product_gas = gd.fluid(name='product-gas',gamma=k,R=R_specific,units='J / kg-K')
    for i in np.arange(len(chamber.index)):
        M = gd.mach_from_area_ratio(chamber.at[i,'eps'],gas=product_gas)
        if chamber.at[i,'regime'] == 0.0:
            chamber.at[i,'mach'] = M[0]
        elif chamber.at[i,'regime'] == 2.0:
            chamber.at[i,'mach'] = M[1]
        else:
            chamber.at[i,'mach'] = 1
        T,P = t_p_from_mach(Tc,Pc,chamber.at[i,'mach'],k)
        chamber.at[i,'T'] = T
        chamber.at[i,'P'] = P


    return chamber


def t_p_from_mach(T0,P0,M,k):
    T = T0/(1+0.5*(k-1)*M**2)
    P = P0/((1+0.5*(k-1)*M**2)**(k/(k-1)))
    return T,P