"""
Functions for calculating thermal properties in combustor based on geometry

author: tcharlson
"""
import numpy as np
import gas_dynamics as gd
from ethanol_props import get_ethanol_props_SI

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

def HT_1D_solve_bartz(eps,M,T0,T_aw,k,b1,b2,b3,b4,
                      t_w,cond_w,
                      convergence = 0.005,
                      Tguess = 700):

    """
    1D heat transfer analysis at given station, given input conditions

    eps = Area Ratio @ current station
    M = mach number @ current station
    T0 = combustor stagnation temp [K] #~3200 for current design
    T_aw = combustor adiabatic wall temp, [K]
    b1..4 = bartz correlation paramters, computed as constants from stagnation conditions

    t_w = wall thickness [m]
    cond_w = wall conductivity [W/m.K] # 370 for copper (lowball)

    Tguess = 700; seed temperature to start iteration
    """
    # Prepare bartz calculation parameters
    AR_factor = (1/eps)**0.9 # expansion factor, (A*/A) ** 0.9
    b_const = b1*b2*b3*b4 # calculate lumped bartz constant
    omega = bartz_correction_factor(T0,Tguess,M,k)

    T_wg = Tguess # for now, assume whatever was guessed

    h_wg = b_const*AR_factor*omega
    q_wg = h_wg*(T_aw - T_wg)

    T_wc = T_wg - (q_wg*t_w)/cond_w

    return q_wg,T_wg,T_wc

def bartz_correction_factor(T0,Twg,M,k,
                            omega=0.6):
    """
    return bartz correction factor, sigma, given inputs of gamma, Mach number, and temperature ratio

    use default omega = 0.6
    refs: file:///C:/Users/charl/Documents/untitled_rocket_project/resources/bartz_RocketHT.pdf
    """
    exp1 = 0.8 - (omega/5) # first exponent
    exp2 = omega/5 # second exponent

    sigma1 = (0.5*(Twg/T0)*(1+((k-1)/2)*(M**2)) + 0.5)**exp1 # calculate first half of denominator
    sigma2 = (1+((k-1)/2)*(M**2))**exp2 # calculate 2nd half of denominator

    return 1/(sigma1*sigma2)

def t_adiabatic_wall(T0,Pr,M,k):

    kf = (k-1)/2
    r = Pr**0.33 # for turbulent flow

    numerator = 1 + r*kf*M**2
    denominator = 1 + kf*M**2

    return T0*(numerator/denominator)





