"""
Constant params to import as part of higher level vehicle sizing script
"""
import numpy as np

# Conversion Factors #
BAR2PSI = 14.5038
METERS2FEET = 3.281
NEWTON2LBF = 0.224809

# Propellants:
rho_o = 1141.0 #kg/m**3
rho_f = 800.0 # kg/m**3

# Pressurant
T_p_LOAD = 300 #K - pressurant load temp
press_margin = 5.0 #bar - WAG - need to anchor to avail Kv
vol_sweep = np.linspace(10.0,60.0,num=3)/1000 # liters -> m**3

# Pressure Ladder #
f_inj_stiff = 20.0 # percent
f_reg_stiff = 15.0 # percent
o_inj_stiff = 20.0 # percent

# Feed Losses/appx len
fs_length_o = 0.5
fs_length_f = 1.5
floss_o = 0.271 #bar/meter
floss_f = 0.257 #bar/meter
lineloss_o = floss_o*fs_length_o
lineloss_f = floss_f*fs_length_f

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

## PLUMBING ##
"""
Aqua dome reg: 1783 - 0.76 kg
Aqua pilot reg: 1247-1 - 0.15 kg
Aqua check:  684 - 0.14 kg
"""

beta_valve = 2.0
m_valves = 2.1*beta_valve # kg

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

# Nosecone - WAG for now
m_nc = 2.0 #kg

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

6" SCH5 SS
OD = 6.625" => 168.3mm
wt = 2.77mm (0.109")
lin_wt = 11.29

"""

tank_OD = 168.3 # mm
tank_WT = 3.4 # mm
tank_ID = tank_OD-(2*tank_WT)
tank_LW = 13.91

tank_AInternal = np.pi*(tank_ID/(1000*2))**2 # m**2

ullage_frac_o = 1.3
ullage_frac_f = 1.15
