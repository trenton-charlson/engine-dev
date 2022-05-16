"""
Functions for calculating thermal properties in combustor based on geometry

author: tcharlson
"""
import os

import numpy as np
import gas_dynamics as gd
import matplotlib.pyplot as plt
import yaml

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

    t_w, cond_w, bartz_mult, T0, T_c_i, h_c, T_aw, M, k, R = data

    return [X[2] - h_c*(X[1] - T_c_i),
            X[2] - (1/(R+1/(bartz_mult*bartz_correction_factor(T0, X[0], M, k))))*(T_aw - X[0]),
            X[2] - (cond_w/t_w)*(X[0] - X[1])]

def soot_thermal_resistance(eps,regime):
    if regime == 0.0:
        # subsonic
        if eps>1.71:
            R = 0.000526667590690191
        else:
            R = 0.000194337986*eps + 0.000177914
    elif regime == 1.0:
        # throat
        R = 0.000384911850873474
    else:
        # supersonic
        R = 5.20899E-05*eps + 0.00033634375

    return R

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


def load_regen_geometry(chamber,fn):
    with open(os.path.join(os.getcwd(),'CONFIG',fn)) as f:
        regen_config = yaml.load(f, Loader=yaml.FullLoader)

    for i in chamber.index:
        for sect in regen_config.keys():
            x = round(chamber.at[i,'x'],2)
            if (x >= regen_config[sect]['start_x'] and x <= regen_config[sect]['end_x']):
                # import config params
                chamber.at[i,'n_chan'] = regen_config[sect]['n_chan']
                chamber.at[i,'w_chan'] = regen_config[sect]['w_chan']
                chamber.at[i,'d_chan'] = regen_config[sect]['d_chan']
                chamber.at[i,'t_wall'] = regen_config[sect]['t_wall']
                # calculate addtl channel params
                chamber.at[i,'r_cool'] = chamber.at[i,'t_wall'] + chamber.at[i,'r']
                chamber.at[i,'r_outer'] = chamber.at[i,'r_cool'] + chamber.at[i,'d_chan']
                chamber.at[i,'fin_thick'] = (2*np.pi)*chamber.at[i,'r_cool']/chamber.at[i,'n_chan'] - chamber.at[i,'w_chan']

    return chamber

def plot_chamber_thermo(chamber,eng):
    i_n = chamber.index[-1] #grab max index

    fig1, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(nrows=3, ncols=2, figsize=(16,9), sharex=True)
    ax1.plot(chamber['x'], chamber['r'], c='k', lw=2)
    ax1.plot(chamber['x'], chamber['r_cool'], c='k', lw=1, ls=':')
    ax1.plot(chamber['x'], chamber['r_outer'], c='k', lw=1, ls='-')
    ax1.set(ylim=(0.0, max(chamber['r'] * 1.25)))
    ax11 = ax1.twinx()
    ax11.plot(chamber['x'], chamber['q_tot']/10**7, c='r', lw=2, label='q_tot [MW/m^2]')

    # plot temps:
    ax2.plot(chamber['x'], chamber['T_c_i'], c='b', lw=2,label='T_c')
    ax2.plot(chamber['x'], chamber['T_wg'], c='r',lw=2,label='T_wg')
    ax2.plot(chamber['x'], chamber['T_wc'], c='cyan',lw=2,label='T_wc')
    ax21 = ax2.twinx()
    ax21.plot(chamber['x'], chamber['dT_c']/chamber['s'], c='magenta', lw=2, label='dT_c/dx')
    ax21.plot(chamber['x'], 10*chamber['dP_c']/chamber['s'], c='green', lw=2, label='dP_c/dx')

    ax3.plot(chamber['x'], chamber['mach'], c='r', lw=2, label='Mach Number - [-]')
    ax31 = ax3.twinx()
    ax31.plot(chamber['x'], chamber['T_aw'], label='Adiabatic Wall Temp [K]')

    ax4.plot(chamber['x'], chamber['t_wall'], label='t_wall')
    ax4.plot(chamber['x'], chamber['w_chan'], label='w_chan')
    ax4.plot(chamber['x'], chamber['d_chan'], label='d_chan')
    ax41 = ax4.twinx()
    ax41.plot(chamber['x'], chamber['n_chan'], c='k',ls='--',label='n_chan')

    ax5.plot(chamber['x'], chamber['Re_c'], c='b', label='Coolant Reynolds Number')
    ax51 = ax5.twinx()
    ax51.plot(chamber['x'], chamber['rho_c'], c='cyan', label='Coolant Density [kg/m**3]')

    ax6.plot(chamber['x'], chamber['fin_thick'], c='r', label = 'Fin Wall Thickness [mm]')
    ax61 = ax6.twinx()
    ax61.plot(chamber['x'], chamber['u_c'], c='k', ls='-', label='Coolant Velocity - [m/s]')
    ax61.plot(chamber['x'], chamber['P_c_e'], c='cyan', ls='-', label='Coolant Pressure - [Bar]')


    # grids per Brandon Kan
    ax1.grid()
    ax1.set_xlim((0,chamber.at[i_n,'x']))

    ax2.grid()
    ax2.set_xlim((0,chamber.at[i_n,'x']))

    ax3.grid()
    ax3.set_xlim((0,chamber.at[i_n,'x']))

    ax4.grid()
    ax4.set_xlim((0,chamber.at[i_n,'x']))

    ax5.grid()
    ax5.set_xlim((0,chamber.at[i_n,'x']))

    ax6.grid()
    ax6.set_xlim((0,chamber.at[i_n,'x']))

    #ax1.axis('equal')
    ax11.legend()
    ax31.legend()
    ax21.legend()
    ax41.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()
    ax5.legend()
    ax51.legend()
    ax61.legend()
    ax5.set_xlabel('x - [mm]')
    ax6.set_xlabel('x - [mm]')
    ax1.set_ylabel('r - [mm]')
    ax4.set_ylabel('Channel Geo - [mm]')
    ax41.set_ylabel('Channel Count - [n]')
    ax2.set_ylabel(f'Regen Temps - [K]')
    ax21.set_ylabel(f'Coolant dT/dx., dP/dx - [K/mm], [BarA/CM]')
    ax3.set_ylabel(f'Regen Geo - [mm]')
    ax6.set_ylabel(f'Fin Thickness - [mm]')
    ax61.set_ylabel(f'Velocity, Pressure - [m/s], [BarA]')
    regen_dP = np.round(max(chamber["P_c_e"])-min(chamber["P_c_e"]),2)
    pct_stiffness = regen_dP/max(chamber["P_c_e"])
    fig1.suptitle(f'Regen Cooling Properties - PC: {eng.Pc} [BarA]; MR: {eng.MR} [-]; {eng.thrust} [N]; e_c: {eng.fac_CR} [-]\n'
                  f'Coolant Flowrate: {np.round(chamber.at[0,"mdot_chan"]*chamber.at[0,"n_chan"],3)} [kg/s]\n'
                  f'Max Coolant Temp: {np.round(max(chamber["T_c_i"]),2)} [K]  -  Coolant dT: {np.round(max(chamber["T_c_i"]) - min(chamber["T_c_i"]),2)} [K]\n'
                  f'Max Wall Temp: {np.round(max(chamber["T_wg"]),2)} [K]\n'
                  f'Transit Time: {np.round(np.sum(chamber["t_transit"])*1000,2)} [ms]  -  Transit Dist: {np.round(np.sum(chamber["s"]),2)} [mm]\n'
                  f'Regen dP = {regen_dP} [BarA] - {np.round(pct_stiffness*100,1)}%')
    fig1.tight_layout
    return fig1





