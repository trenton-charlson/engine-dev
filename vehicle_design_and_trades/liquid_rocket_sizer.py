"""
Size out rocket design based on high level input parameters

future: perform trajectory sim & iterative analysis
"""
import copy

import numpy as np

from feed_system_design_code import gas_dyn_utils
from vehicle_sizing_functions import mass_tank_segment

BAR2PSI = 14.5038

## TOP LEVEL PARAMS ##
"""
thrust = 3000 N
burntime = 15 sec
alt_target = 

mdot_o
mdot_f
"""
PC = 20.0 # bars
thrust = 3000 # newton
burntime = 15.0 # seconds
PRESSGASS = 'N2'

# Pressure Ladder
f_inj_stiff = 20.0 # percent
f_reg_stiff = 15.0 # percent
o_inj_stiff = 20.0 # percent

# Engine Flowrates
mdot_o = 0.7821 # kg/s
rho_o = 1141.0 #kg/m**3
q_ox = mdot_o/rho_o # m**3/s
m_o_i = mdot_o*burntime

mdot_f = 0.5291 # kg/s
rho_f = 800.0 # kg/m**3
q_f = mdot_f/rho_f # m**3/s
m_f_i = mdot_f*burntime


# Engine/Feed Pressures
P_f_inj = (1+(f_inj_stiff/100))*PC
P_f_inlet = (1+(f_reg_stiff/100))*P_f_inj
P_o_inj = (1+(o_inj_stiff/100))*PC
P_o_inlet = copy.copy(P_o_inj)

# Feed Losses/appx len
fs_length_o = 0.5
fs_length_f = 1.5
floss_o = 0.271 #bar/meter
floss_f = 0.257 #bar/meter

lineloss_o = floss_o*fs_length_o
lineloss_f = floss_f*fs_length_f

# Tankage Pressures
P_o_tank = P_o_inlet+lineloss_o
P_f_tank = P_f_inlet+lineloss_f

# PRessurant
press_margin = 5.0 #bar - WAG - need to anchor to avail Kv
P_p_EOL = np.round(max(P_f_tank,P_o_tank) + press_margin)
P_p_BOL = np.round(4500/BAR2PSI,2) # beginning bottle pressure
T_p_LOAD = 300 #K - pressurant load temp
vol_sweep = np.linspace(10.0,60.0,num=10)/1000 # liters -> m**3

## TANKAGE ##
"""
4" SCH10 SS:
OD = 4.500 => 114.3mm
wt = 3.05mm
lin_wt = 8.42 kg/m

5" SCH10 SS:
OD = 5.563" => 141.3mm
wt = 3.4mm
lin_wt = 11.64 kg/m

6" SCH10 SS:
OD = 6.625" => 168.3mm
wt = 3.4mm
lin_wt = 13.91 kg/m

"""

tank_OD = 1
tank_ID = 1
tank_LW = 1

## AIRFRAME ##
"""
-- 7.5" fiberglass tubing --
ID: 7.520" => 191 mm
OD: 7.708" => 195.78 mm
lin_wt = 18.8 oz/ft = 1.175lb/ft => 0.533 kg/ft => 1.748 kg/m
^^ this is less than 6"....
spoof in 3kg/m for now
==> 1/10.755 conversion factor from oz/ft -> kg/m

-- 6" fiberglass tubing
ID: 6.00" => 152.4 mm
OD: 6.17" => 156.718mm
lin_wt = 24.30 oz/ft 
"""
skin_ID = 191 # mm
skin_OD = 195.78 # mm
skin_LW = 1.748 # kg/m


## THRUST STRUCTURE ##
"""
Handcalc from 05/07/22 - ~0.95kg
"""
m_thrust_struct = 0.95 # kg
l_thrust_struct = 5.0*25.4/1000 # m baseline

## ENGINE ##
"""
Handcalc from 05/07/22 - ~4kg
"""
m_engine = 4.0 # kg
l_engine = 0.230 # meters

blowdown = gas_dyn_utils.blowdown_sensitivity_study(vol_sweep,[P_p_BOL],T_p_LOAD,P_p_EOL,
                                                    P_o_tank, q_ox,
                                                    P_f_tank, q_f,
                                                    burntime,
                                                    PLOT=False)

V_p_BOL = blowdown[f'{PRESSGASS} Vol Required']
m_p_i = blowdown[f'{PRESSGASS} Mass, Initial']
m_p_f = blowdown[f'{PRESSGASS} Mass Residual']

