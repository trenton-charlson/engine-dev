"""
Utility functions for gas dynamics shiz
"""

import CoolProp.CoolProp as CP
import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import itertools

def m3s_2_sm3h(rate,P,T,
               BAR=True):
    """
    convert m**3/s @ pressure to Std m**3/hr
    """
    if BAR:
        std_p_bar = 1.01325
    else:
        std_p_bar = 101325 # Pa

    rate_stp = (rate * (273.15 / T) * (P / std_p_bar)) * 3600  # m**3 /s

    return rate_stp

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
                              ts=0.1, t_max=100.0, BAR=10**5, K = 1.4,
                              DEBUG=False):
    t=0.0
    P = P_start
    out = pd.DataFrame()

    # seed starting params
    out.at[t, 'P_u'] = P
    out.at[t, 'T_u'] = T_start
    out.at[t, 'rho_u'] = CP.PropsSI('D','T', T_start, 'P', P*BAR, GAS)
    mass_i = V_tank*out.at[t, 'rho_u']
    out.at[t, 'mass_u'] = mass_i

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
    t_blowdown = max(out.index) # total blowdown time avail
    mass_f = min(out['mass_u'])
    if DEBUG:
        print(f'Blowdown time from {P_start} Bar -> {P_end} Bar = {t_blowdown} [s]')

    return out, t_blowdown, mass_i, mass_f


def blowdown_sensitivity_study(vol,p_start,T_bulk,P_end,
                               P_ot,q_dot_oto,
                               P_ft,q_dot_fto,
                               burntime,
                               PLOT=True):

    print(f'Performing Pressurant System Sensitivity Analysis for:\n'
          f'P = {p_start}\n'
          f'Volume range: {min(vol)*1000} -> {max(vol)*1000} [L] - npts: {len(vol)}\n\n')

    markers = itertools.cycle((',', '+', '.', 'o', '*'))
    if PLOT:
        fig1, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(9,12), sharex=False)

    out_df = pd.DataFrame(index=p_start)

    for p in p_start:
        N2_df = pd.DataFrame(index=vol)
        He_df = pd.DataFrame(index=vol)
        marker = next(markers)

        print(f'Solving for {p} Bar Bottle Pressure')

        for v in tqdm(vol):
            # calc for nitrogen
            _, time, mass_i, mass_f = regulator_blowdown_rocket(p, T_bulk, P_end,
                                                      P_ot, q_dot_oto,
                                                      P_ft, q_dot_fto,
                                                      v,
                                                      'Nitrogen')
            N2_df.at[v,'time'] = time
            N2_df.at[v,'mass_i'] = mass_i
            N2_df.at[v,'mass_f'] = mass_f
            # calc for Helium
            _, time, mass_i, mass_f = regulator_blowdown_rocket(p, T_bulk, P_end,
                                                      P_ot, q_dot_oto,
                                                      P_ft, q_dot_fto,
                                                      v,
                                                      'Helium')
            He_df.at[v, 'time'] = time
            He_df.at[v, 'mass_i'] = mass_i
            He_df.at[v, 'mass_f'] = mass_f

            # Calculate Required Volume to meet burntime:
            n2_vol_req = np.interp(burntime,N2_df['time'],N2_df.index*1000) #  required starting vol to meet burntime, Liters
            he_vol_req = np.interp(burntime,He_df['time'],He_df.index*1000) #  required starting vol to meet burntime, Liters
            out_df.at[p,'N2 Vol Required'] = n2_vol_req
            out_df.at[p,'N2 Mass, Initial'] = np.interp(n2_vol_req,N2_df.index*1000,N2_df['mass_i'])
            out_df.at[p,'N2 Mass Residual'] = np.interp(n2_vol_req,N2_df.index*1000,N2_df['mass_f'])
            out_df.at[p,'He Vol Required'] = he_vol_req
            out_df.at[p,'He Mass, Initial'] = np.interp(he_vol_req,He_df.index*1000,He_df['mass_i'])
            out_df.at[p,'He Mass Residual'] = np.interp(he_vol_req,He_df.index*1000,He_df['mass_f'])

        if PLOT:
            ax1.plot(N2_df.index * 1000, N2_df['time'], marker=marker, c='g', label=f'Nitrogen - {p} Bar')
            ax1.plot(He_df.index * 1000, He_df['time'], marker=marker, c='magenta', label=f'Helium - {p} Bar')
            ax2.plot(N2_df.index * 1000, N2_df['mass_i'], marker=marker, c='g', label=f'Nitrogen - {p} Bar')
            ax2.plot(He_df.index * 1000, He_df['mass_i'], marker=marker, c='magenta', label=f'Helium - {p} Bar')

    out_df = np.round(out_df,2)

    if PLOT:
        ax1.axhline(y=burntime, c='k', ls='--', label="Burntime")
        ax1.legend()
        ax1.grid()
        ax1.set_ylabel(f'Blowdown Time for {P_end} Bar EOL Press [s]')

        ax2.legend()
        ax2.grid()
        ax2.set_xlabel('Bottle Volume [L]')
        ax1.set_xlim(min(vol)*1000, max(vol)*1000)
        ax2.set_xlim(min(vol)*1000, max(vol)*1000)
        ax1.set_xticks(np.arange(min(vol)*1000, max(vol)*1000, 5.0))
        ax1.tick_params(labelbottom=False)
        ax2.set_xticks(np.arange(min(vol)*1000, max(vol)*1000, 5.0))
        ax2.set_ylabel(f'Loaded BOL Mass [kg]')

        cell_text = []
        for row in range(len(out_df)):
            cell_text.append(out_df.iloc[row])

        table = ax3.table(cellText=cell_text, colLabels=out_df.columns, loc='center')
        ax3.axis('off')
                          #colWidths=np.ones(len(out_df.columns))*0.1)
        #table.auto_set_font_size(False)
        #table.set_fontsize(10.0)

        # Set title and show plot
        fig1.suptitle(f'Blowdown System Analysis\n'
                      f'p = {p_start} [bar]')
        plt.tight_layout

    return out_df
