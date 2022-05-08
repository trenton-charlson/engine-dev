import copy
import os

import numpy as np
import pandas as pd
from scipy.optimize import fsolve

from gas_dyn_utils import isenthalpic_throttle, regulator_blowdown_single_species, regulator_blowdown_rocket

import CoolProp.CoolProp as CP

GAS = 'Helium'
T_bulk = 300 #K

#P_start = 131.0 #bar - ~1900 psi bottle load pressure
P_start = 250.0 #bar - ~4500 psi bottle load pressure


BAR = 10**5

std_p_bar = 1.01325
#P_ft = 43.437 #bars - ~630 psia - wag for 30 bar PC
P_ft = 36.26 #bars - ~550 psia - wag for 25 bar PC
P_ot = 33.0 #bars ~480 psi
ts = 0.1
t = 0

MOL_stp = 22.4

mdot_fto = 0.439 #kg/s
rho_f = 800.0 #kg/m**3
q_dot_fto = mdot_fto/rho_f # m**3 / s

mdot_oto = 0.73 #kg/s
rho_o = 1141.0 #kg/m**3
q_dot_oto = mdot_oto/rho_o # m**3 / s

T_ds = 280.0 # K - donstream temp - hack for now

q_dot_s_fto_STP = q_dot_fto*(273.15/T_ds)*(P_ft/std_p_bar) # m**3 /s
s2h = 3600 # seconds per hour
q_dot_h_fto_STP = q_dot_s_fto_STP*s2h

print(f'Required Gas Volume Flowrates:'
      f'Raw Flowrate = {q_dot_fto} m**3 / s\n'
      f'STP Flowrate (s) = {q_dot_s_fto_STP} m**3 / s @ STP\n'
      f'STP Flowrate (h) = {q_dot_h_fto_STP} m**3 / hr @ STP\n')


rho_n2 = CP.PropsSI('D', 'T', T_bulk, 'P', P_start*BAR, GAS) #kg/m**3
rho_n2_stp = CP.PropsSI('D', 'T', 273.15, 'P', 101325, GAS) # kg/m**3 - STP

print(f'Bulk starting density: {rho_n2} kg/m**3')

V_tank = 40/1000 # L -> m**3
BOL_mass = V_tank*rho_n2
mass = copy.copy(BOL_mass)

A_or = np.pi*((0.07*25.4)/(2*1000))**2
Cd = 0.6
CdA = Cd*A_or

P_lo = 50
P_hi = 120

Kv_reqs = pd.DataFrame(index=range(P_lo,P_hi))

for eol_p in Kv_reqs.index:
    Kv_reqs.at[eol_p, 'kv_req'] = 0.0019 * q_dot_h_fto_STP * (rho_n2_stp*T_bulk/(P_ft*(eol_p-P_ft)))**0.5
    Kv_reqs.at[eol_p, 'T2'], Kv_reqs.at[eol_p, 'rho2'] = isenthalpic_throttle(eol_p,P_ft,T_bulk,GAS)

P_end = 45.0
blowdown = regulator_blowdown_single_species(P_start, T_bulk, P_end, P_ft, q_dot_fto, V_tank, GAS)

blowdown_biprop = regulator_blowdown_rocket(P_start,T_bulk,P_end,
                                            P_ot, q_dot_oto,
                                            P_ft, q_dot_fto,
                                            V_tank,
                                            GAS)


