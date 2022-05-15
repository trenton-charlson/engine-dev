

import numpy as np
import pandas as pd

from gas_dyn_utils import isenthalpic_throttle, regulator_blowdown_rocket, m3s_2_sm3h, blowdown_sensitivity_study

import CoolProp.CoolProp as CP

GAS = 'He'
T_bulk = 300 #K
burntime = 16

#P_start = 131.0 #bar - ~1900 psi bottle load pressure
P_start = 250.0 #bar - ~4500 psi bottle load pressure

BAR = 10**5
std_p_bar = 1.01325
#P_ft = 43.437 #bars - ~630 psia - wag for 30 bar PC
P_ft = 31.0 #bars
P_ot = 25.0 #bars
ts = 0.1
t = 0.0

MOL_stp = 22.4

mdot_fto = 0.529 #kg/s
rho_f = 800.0 #kg/m**3
q_dot_fto = mdot_fto/rho_f # m**3 / s

mdot_oto = 0.782 #kg/s
rho_o = 1141.0 #kg/m**3
q_dot_oto = mdot_oto/rho_o # m**3 / s

T_ds = 280.0 # K - donstream temp - hack for now

qdot_f_stp = m3s_2_sm3h(q_dot_fto,P_ft,260.0)
qdot_o_stp = m3s_2_sm3h(q_dot_oto,P_ot,230.0)


print(f'Required Gas Volume Flowrates:'
      f'Fuel Raw Flowrate = {q_dot_fto} m**3 / s\n'
      f'Fuel STP Flowrate (h) = {qdot_f_stp} m**3 / hr @ STP\n'
      f'Ox Raw Flowrate = {q_dot_oto} m**3 / s\n'
      f'Ox STP Flowrate (h) = {qdot_o_stp} m**3 / hr @ STP\n')


rho_n2_stp = CP.PropsSI('D', 'T', 273.15, 'P', 101325, 'Nitrogen') # kg/m**3 - STP
rho_he_stp = CP.PropsSI('D', 'T', 273.15, 'P', 101325, 'Helium')

V_tank = 40/1000 # L -> m**3

A_or = np.pi*((0.07*25.4)/(2*1000))**2
Cd = 0.6
CdA = Cd*A_or

P_lo = 40
P_hi = 120

Kv_reqs = pd.DataFrame(index=range(P_lo,P_hi))

for eol_p in Kv_reqs.index:
    Kv_reqs.at[eol_p, 'kv_req_f_N2'] = 0.0019 * qdot_f_stp * (rho_n2_stp*T_bulk/(P_ft*(eol_p-P_ft)))**0.5
    Kv_reqs.at[eol_p, 'kv_req_o_N2'] = 0.0019 * qdot_o_stp * (rho_n2_stp*T_bulk/(P_ot*(eol_p-P_ot)))**0.5
    Kv_reqs.at[eol_p, 'kv_req_f_He'] = 0.0019 * qdot_f_stp * (rho_he_stp*T_bulk/(P_ft*(eol_p-P_ft)))**0.5
    Kv_reqs.at[eol_p, 'kv_req_o_He'] = 0.0019 * qdot_o_stp * (rho_he_stp*T_bulk/(P_ot*(eol_p-P_ot)))**0.5

    #Kv_reqs.at[eol_p, 'T2'], Kv_reqs.at[eol_p, 'rho2'] = isenthalpic_throttle(eol_p,P_ft,T_bulk,GAS)


Kv_reqs.plot()

P_end = 45.0

#blowdown = regulator_blowdown_single_species(P_start, T_bulk, P_end, P_ft, q_dot_fto, V_tank, GAS)

blowdown_biprop,t,m,_ = regulator_blowdown_rocket(P_start,T_bulk,P_end,
                                                  P_ot, q_dot_oto,
                                                  P_ft, q_dot_fto,
                                                  V_tank,
                                                  GAS)

# Blowdown Sensitivity Study #
vol = np.linspace(10,100,num=10)/1000 # convert to m**3
p_start = [130.0,200.0,300.0,400.0]

blowdown = blowdown_sensitivity_study(vol,p_start,T_bulk,P_end,
                                      P_ot, q_dot_oto,
                                      P_ft, q_dot_fto,
                                      burntime)

print(blowdown)



