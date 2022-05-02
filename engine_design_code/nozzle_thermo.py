"""
Functions for calculating thermal properties in combustor based on geometry

author: tcharlson
"""
import numpy as np
import gas_dynamics as gd
import matplotlib.pyplot as plt
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



def func_bartz_root_find(X, *data):
    #X[0] = T_wg
    #X[1] = T_wc
    #X[2] = q

    t_w, cond_w, bartz_mult, T0, T_c_i, h_c, T_aw, M, k = data


    return [X[2] - h_c*(X[1] - T_c_i),
            X[2] - bartz_mult*bartz_correction_factor(T0, X[0], M, k)*(T_aw - X[0]),
            X[2] - (cond_w/t_w)*(X[0] - X[1])]


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

def gnielinski_calc(f,Re,Pr):
    # return nusselt number as result of Gnielinski correlation
    # Incorpera pg. 515

    numerator = (f/8)*(Re-1000)*Pr
    denominator = 1 + ((12.7*(f/8)**0.5)*((Pr**(2/3))-1))

    return numerator/denominator

def plot_chamber_thermo(chamber):
    fig1, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, sharex=True)
    ax1.plot(chamber['x'], chamber['r'], c='k', lw=2)
    ax1.set(ylim=(0.0, max(chamber['r'] * 1.25)))
    ax11 = ax1.twinx()
    ax11.plot(chamber['x'], chamber['q_tot']/10**7, c='r', lw=2, label='q_tot [MW/m^2]')

    # plot temps:
    ax2.plot(chamber['x'], chamber['T_c_i'], c='b', lw=2,label='T_c')
    ax2.plot(chamber['x'], chamber['T_wg'], c='r',lw=2,label='T_wg')
    ax2.plot(chamber['x'], chamber['T_wc'], c='cyan',lw=2,label='T_wc')

    ax3.plot(chamber['x'], chamber['t_wall'], label='t_wall')
    ax3.plot(chamber['x'], chamber['w_chan'], label='w_chan')
    ax3.plot(chamber['x'], chamber['d_chan'], label='d_chan')

    ax1.axis('equal')
    ax11.legend()
    ax2.legend()
    ax3.legend()
    ax3.set_xlabel('x - [mm]')
    ax1.set_ylabel('r - [mm]')
    ax2.set_ylabel(f'Regen Temps - [K]')
    ax3.set_ylabel(f'Regen Geo - [mm]')
    fig1.suptitle('Regen Cooling Properties')
    return





