"""
Size out rocket design based on high level input parameters

future: perform trajectory sim & iterative analysis
"""
import copy
import matplotlib.pyplot as plt
import numpy as np

from feed_system_design_code import gas_dyn_utils
from vehicle_sizing_functions import mass_tank_segment

from _1D_rocket_traj import _1D_rocket_traj

BAR2PSI = 14.5038
METERS2FEET = 3.281

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
V_ox_i = q_ox*burntime
m_o_i = mdot_o*burntime

mdot_f = 0.5291 # kg/s
rho_f = 800.0 # kg/m**3
q_f = mdot_f/rho_f # m**3/s
V_f_i = q_f*burntime
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

# Pressurant
press_margin = 5.0 #bar - WAG - need to anchor to avail Kv
P_p_EOL = np.round(max(P_f_tank,P_o_tank) + press_margin)
P_p_BOL = np.round(4500/BAR2PSI,2) # beginning bottle pressure
T_p_LOAD = 300 #K - pressurant load temp
vol_sweep = np.linspace(10.0,60.0,num=3)/1000 # liters -> m**3

## ENGINE ##
"""
Handcalc from 05/07/22 - ~4kg
"""
m_engine = 4.0 # kg
l_engine = 0.230 # meters

## THRUST STRUCTURE ##
"""
Handcalc from 05/07/22 - ~0.95kg
"""
m_thrust_struct = 0.95 # kg
l_thrust_struct = 5.0*25.4/1000 # m baseline

## Aft Valve Pack ##
m_aft_vp = 2.5 # kg
l_aft_vp = 8.0*25.4/1000

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

A_cs = np.pi*(skin_OD/(2*1000))**2

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

tank_OD = 168.3 # mm
tank_WT = 3.4 # mm
tank_ID = tank_OD-(2*tank_WT)
tank_LW = 13.91

tank_AInternal = np.pi*(tank_ID/(1000*2))**2 # m**2

ullage_frac_o = 1.3
ullage_frac_f = 1.15

L_o_tank = (V_ox_i*ullage_frac_o)/tank_AInternal
m_o_tank = mass_tank_segment(L_o_tank,tank_LW,skin_LW)
L_f_tank = (V_f_i*ullage_frac_f)/tank_AInternal
m_f_tank = mass_tank_segment(L_f_tank,tank_LW,skin_LW)

## SIZE PRESSURANT SYSTEM ##
"""
Use blowdown sim to calculate time req given starting pressure target
"""
blowdown = gas_dyn_utils.blowdown_sensitivity_study(vol_sweep,[P_p_BOL],T_p_LOAD,P_p_EOL,
                                                    P_o_tank, q_ox,
                                                    P_f_tank, q_f,
                                                    burntime,
                                                    PLOT=True)

V_p_BOL = blowdown.at[np.round(P_p_BOL,2),f'{PRESSGASS} Vol Required']
m_p_i = blowdown.at[np.round(P_p_BOL,2),f'{PRESSGASS} Mass, Initial']
m_p_f = blowdown.at[np.round(P_p_BOL,2),f'{PRESSGASS} Mass Residual']
mdot_p = (m_p_i-m_p_f)/burntime # average pressurant flowrate

m_propellant_i = m_f_i + m_o_i
m_struct = m_thrust_struct+m_engine+m_o_tank+m_f_tank

m_wet = m_propellant_i+m_struct+m_p_i
m_dry = m_struct+m_p_i # assume worst case all pressurant stays onboard 

mdot_t = mdot_f+mdot_o+mdot_p # total propellant flowrate


print(f'VEHICLE SIZING OUTPUTS:\n'
      f'Propellant Mass: {m_propellant_i} [kg]\n'
      f'Pressurant Mass: {m_p_i} [kg]\n'
      f'Structural Mass: {np.round(m_struct,2)} [kg]\n\n'
      f'>> TOTAL LIFTOFF MASS = {np.round(m_wet,2)} [kg]\n'
      f'>> TOTAL DRY MASS = {np.round(m_dry,2)} [kg]\n')

## Run Traj Sim ##

traj = _1D_rocket_traj(m_wet,m_dry,mdot_t,thrust,A_cs)

fig, ax = plt.subplots(figsize=(16,6))
left = 0.0
for element in [l_engine,l_thrust_struct,l_aft_vp,L_o_tank,0.308,L_f_tank]:
    ax.barh(0.0,element,height=skin_OD/1000,left=left,label=str(np.round(element,3)))
    left = left+element

ax.legend()
ax.set_title('Approximate Rocket Dimensions\n'
             f'DIA = {skin_OD} [mm]  --  {skin_OD/25.4} [in]\n'
             f'LEN = {left} [m]  --  {left*METERS2FEET} [ft]\n'
             f'Aspect Ratio: {left/(skin_OD/1000)} [-]')
ax.axis('equal')
plt.tight_layout