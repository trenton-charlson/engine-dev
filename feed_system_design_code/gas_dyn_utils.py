"""
Utility functions for gas dynamics shiz
"""

import CoolProp.CoolProp as CP
import pandas as pd
import numpy as np


def vdot_to_flow_coeff():


    return

def isenthalpic_throttle(P1,P2,T1,GAS,
                         BAR = 10**5):
    h1 = CP.PropsSI('H', 'T', T1, 'P', P1*BAR, GAS)
    T2 = CP.PropsSI('T', 'H', h1, 'P', P2*BAR, GAS)
    rho2 = CP.PropsSI('D', 'T', T2, 'P', P2*BAR, GAS)

    return T2, rho2


def regulator_blowdown_single_species(P_start, T_start, P_end, P_reg, q_req, V_tank,
                                      GAS,
                                      ts=0.1, t_max=100.0, BAR=10**5, K = 1.4):
    t=0.0
    P = P_start
    out = pd.DataFrame()

    # seed starting params
    out.at[t, 'P_u'] = P
    out.at[t, 'T_u'] = T_start
    out.at[t, 'rho_u'] = CP.PropsSI('D','T', T_start, 'P', P*BAR, GAS)
    out.at[t, 'mass_u'] = V_tank*out.at[t, 'rho_u']

    while (out.at[t,'P_u'] > P_end) and (t <= t_max):
        #isenthalpic throttling process:
        # enthalpy/mass @ starting conditions for this iteration
        hu = CP.PropsSI('H', 'T', out.at[t, 'T_u'], 'P', out.at[t, 'P_u']*BAR, GAS)
        # downstream temp from isenthalpic throttling
        T_d = CP.PropsSI('T', 'H', hu, 'P', P_reg*BAR, GAS)
        # downstream density based on T_d and P_reg
        rho_d = CP.PropsSI('D','T',T_d, 'P', P_reg*BAR, GAS)
        # gas mass flowrate
        mdot = rho_d * q_req

        # populate current timestep
        out.at[t, 'h_u'] = hu
        out.at[t, 'T_d'] = T_d
        out.at[t, 'rho_d'] = rho_d
        out.at[t, 'mdot'] = mdot

        # pop next timestep:
        t = np.round(t + ts, 3) # prevent key errors by rounding float timestamp
        tp = np.round(t-ts,3) # prevent key errors by rounding float timestep

        # calculate new bottle mass after outflow
        out.at[t, 'mass_u'] = out.at[tp, 'mass_u'] - mdot
        # calculate bottle density based on fixed volume
        out.at[t, 'rho_u'] = out.at[t, 'mass_u']/V_tank
        # isentropic/adiabatic expansion in bottle:
        out.at[t, 'T_u'] = out.at[tp, 'T_u'] * (out.at[t, 'rho_u']/out.at[tp, 'rho_u'])**(K-1)
        out.at[t, 'P_u'] = out.at[tp, 'P_u'] * (out.at[t, 'T_u']/out.at[tp, 'T_u'])**(K/(K-1))

    print(f'Blowdown time from {P_start} Bar -> {P_end} Bar = {t} [s]')

    return out


def regulator_blowdown_rocket(P_start,T_start,P_end,
                              P_reg_ox, q_req_ox,
                              P_reg_fuel, q_reg_fuel,
                              V_tank,
                              GAS,
                              ts=0.1, t_max=100.0, BAR=10**5, K = 1.4):
    t=0.0
    P = P_start
    out = pd.DataFrame()

    # seed starting params
    out.at[t, 'P_u'] = P
    out.at[t, 'T_u'] = T_start
    out.at[t, 'rho_u'] = CP.PropsSI('D','T', T_start, 'P', P*BAR, GAS)
    out.at[t, 'mass_u'] = V_tank*out.at[t, 'rho_u']

    while (out.at[t,'P_u'] > P_end) and (t <= t_max):
        #isenthalpic throttling process:
        # enthalpy/mass @ starting conditions for this iteration
        hu = CP.PropsSI('H', 'T', out.at[t, 'T_u'], 'P', out.at[t, 'P_u']*BAR, GAS)

        ### Ox Side ###
        # downstream temp from isenthalpic throttling
        T_d_ox = CP.PropsSI('T', 'H', hu, 'P', P_reg_ox*BAR, GAS)
        # downstream density based on T_d and P_reg
        rho_d_ox = CP.PropsSI('D','T',T_d_ox, 'P', P_reg_ox*BAR, GAS)
        # gas mass flowrate
        mdot_ox = rho_d_ox * q_req_ox

        ### Fuel Side ###
        # downstream temp from isenthalpic throttling
        T_d_f = CP.PropsSI('T', 'H', hu, 'P', P_reg_fuel * BAR, GAS)
        # downstream density based on T_d and P_reg
        rho_d_f = CP.PropsSI('D', 'T', T_d_f, 'P', P_reg_fuel * BAR, GAS)
        # gas mass flowrate
        mdot_f = rho_d_f * q_reg_fuel

        mdot_tot = mdot_f+mdot_ox

        # populate current timestep
        out.at[t, 'h_u'] = hu
        out.at[t, 'T_d_ox'] = T_d_ox
        out.at[t, 'rho_d_ox'] = rho_d_ox
        out.at[t, 'mdot_ox'] = mdot_ox
        out.at[t, 'T_d_f'] = T_d_f
        out.at[t, 'rho_d_f'] = rho_d_f
        out.at[t, 'mdot_f'] = mdot_f
        out.at[t, 'mdot_tot'] = mdot_tot

        # pop next timestep:
        t = np.round(t + ts, 3) # prevent key errors by rounding float timestamp
        tp = np.round(t-ts,3) # prevent key errors by rounding float timestep

        # calculate new bottle mass after outflow
        out.at[t, 'mass_u'] = out.at[tp, 'mass_u'] - mdot_tot
        # calculate bottle density based on fixed volume
        out.at[t, 'rho_u'] = out.at[t, 'mass_u']/V_tank
        # isentropic/adiabatic expansion in bottle:
        out.at[t, 'T_u'] = out.at[tp, 'T_u'] * (out.at[t, 'rho_u']/out.at[tp, 'rho_u'])**(K-1)
        out.at[t, 'P_u'] = out.at[tp, 'P_u'] * (out.at[t, 'T_u']/out.at[tp, 'T_u'])**(K/(K-1))

    out.drop(out.tail(1).index,inplace=True) #  drop last row - https://stackoverflow.com/questions/26921651/how-to-delete-the-last-row-of-data-of-a-pandas-dataframe
    print(f'Blowdown time from {P_start} Bar -> {P_end} Bar = {max(out.index)} [s]')


    return out