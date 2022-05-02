"""
Engine Design Code Top Level Script

author: tcharlson
"""

import copy
import os

import numpy as np
from scipy.optimize import fsolve
import pandas as pd
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj as CEA_Obj_english
from rocketcea.cea_obj_w_units import CEA_Obj

from geometry import chamber_geo, define_contour, plot_chamber_param, import_regen_geo
from nozzle_thermo import chamber_thermo_calcs, t_adiabatic_wall, gnielinski_calc, plot_chamber_thermo, func_bartz_root_find
from ethanol_props import get_ethanol_props_SI

#import gas_dynamics as gd

### DEFINE TOP LEVEL PARAMETERS ###
class propellant:
    def __init__(self,name):
        self.name = name

fuel = propellant("Ethanol")
ox = propellant("LOX")
### PERFORMANCE PARAMS ###
Pc = 35.0 #bar
Pc_pa = Pc * 10**5 # chamber pressure in Pascals
Pe = 1.01325 #bar
thrust = 2500 #newton
MR = 1.4
eta_cstar = 0.9 # guess
R_const = 8.314462*10**3 #J/K.mol
pressure_ratio = Pe/Pc
print(f'pressure_ratio = {pressure_ratio}')
print(f'Pc/Pe = {Pc/Pe}')

### GEOMETRIC PARAMS ###
fac_CR = 10.0 #contraction ratio
chamber_ha = 30.0 #degrees, chamber half angle
nozzle_ha = 15.0 #degrees, nozzle half angle; conical nozzle
L_star = 2.0 #meters - http://mae-nas.eng.usu.edu/MAE_5540_Web/propulsion_systems/section6/section.6.1.pdf guess for now
npts = 100 #number of points to define each contour section

FULL_OUTPUT = False
BAR = True
g = 9.81 # gravitational accel, [m/s**2]


combustor = CEA_Obj(oxName=ox.name, fuelName=fuel.name,
                    cstar_units="m/s",
                    pressure_units="Bar",
                    temperature_units="K",
                    sonic_velocity_units="m/s",
                    enthalpy_units="kJ/kg",
                    density_units="kg/m^3",
                    specific_heat_units="kJ/kg-K",
                    viscosity_units="poise",
                    thermal_cond_units="W/cm-degC",
                    fac_CR=fac_CR)

combustor_std = CEA_Obj_english(oxName=ox.name, fuelName=fuel.name,
                                fac_CR=fac_CR)

#props_out = pd.DataFrame(index=MR_sweep)

print(f'Target Thrust: {thrust} N,\n'
      f'Pc = {Pc} Bar\n'
      f'MR = {MR} [-]\n'
      f'fac_CR = {fac_CR} [-]\n')

(mw_c,k_c) = combustor.get_Chamber_MolWt_gamma(Pc=Pc, MR=MR, eps=40)
(mw_t,k_t) = combustor.get_Throat_MolWt_gamma(Pc=Pc, MR=MR, eps=40)
print(f'CHAMBER - mw: {mw_c}, k: {k_c}')
print(f'THROAT - mw: {mw_t}, k: {k_t}')

k = (k_c+k_t)/2 #average gas constant
mw = (mw_c+mw_t)/2 #avg molecular weight
print(f'K_avg = {k}')

opt_expansion = 1/((((k+1)/2)**(1/(k-1))) * (pressure_ratio**(1/k)) * np.sqrt(((k+1)/(k-1))*(1-(pressure_ratio**((k-1)/k)))))
print(f'Optimum Expansion: {opt_expansion}')

T_c = combustor.get_Tcomb(Pc=Pc, MR=MR)
print(f'T_comb = {T_c} K')

T_t = 2*T_c/(k+1)
print(f'T_t = {T_t} K')

R_specific = R_const/mw
print(f'R_specific = {R_specific} J/kg.K')

v_e_ideal = np.sqrt(((2*k)/(k-1))*R_specific*T_c*(1-(pressure_ratio**((k-1)/k))))
print(f'V_e, ideal = {v_e_ideal} m/s')

cstar_ideal = np.sqrt(k*R_specific*T_c)/(k*np.sqrt((2/(k+1))**((k+1)/(k-1))))
cstar_corr = eta_cstar*cstar_ideal
print(f'C*, ideal = {cstar_ideal} m/s - eta_C* = {eta_cstar}\n'
      f'C*, corr = {cstar_corr} m/s')

mdot_ideal_total = thrust/v_e_ideal
print(f'mdot_ideal_total = {mdot_ideal_total} kg/s\n')

mdot_ox_ideal = (MR/(1+MR))
mdot_fuel_ideal = mdot_ideal_total - mdot_ox_ideal
print(f'mdot_f (ideal) {mdot_fuel_ideal} kg/s')
print(f'mdot_o (ideal) {mdot_ox_ideal} kg/s')


A_t = (mdot_ideal_total/(Pc*10**5))*np.sqrt((R_specific*T_c)/(k*((2/(k+1))**((k+1)/(k-1)))))
#A_t_cm = A_t*(100**2)
A_t_mm = A_t*(1000**2)
R_t_mm = np.sqrt(A_t_mm/np.pi)
D_t_mm = R_t_mm*2
print(f'A_throat = {A_t_mm} mm**2')
print(f'R_throat = {R_t_mm} mm, D_throat = {D_t_mm} - ({D_t_mm/25.4} in)\n')

A_c_mm = fac_CR*A_t_mm
R_c_mm = np.sqrt(A_c_mm/np.pi)
D_c_mm = R_c_mm*2
print(f'A_chamber = {A_c_mm} mm**2')
print(f'R_chamber = {R_c_mm} mm, D_chamber = {D_c_mm} - ({D_c_mm/25.4} in)\n')

A_e_mm = opt_expansion*A_t_mm
R_e_mm = np.sqrt(A_e_mm/np.pi)
D_e_mm = R_e_mm*2
print(f'A_exit = {A_e_mm} mm**2')
print(f'R_exit = {R_e_mm} mm, D_exit = {D_e_mm} - ({D_e_mm/25.4} in)\n')

print(f"Ideal CR (mae_5540) = 8.0/D_t**3/5 + 1.25 = {8.0/((D_t_mm/10)**(3/5))+1.25}")

print(f'sanic {combustor.get_SonicVelocities(Pc=Pc,MR=MR,eps=opt_expansion)}')

print('#### CALCULATING CHAMBER GEOMETRIC PROPERTIES ####\n')
chamber_obj = chamber_geo(A_t_mm, R_t_mm, A_c_mm, R_c_mm, A_e_mm, R_e_mm)
chamber_raw,i_t = define_contour(chamber_obj, L_star, nozzle_ha, chamber_ha)
chamber = copy.deepcopy(chamber_raw) # save raw copy of chamber data for export
chamber_raw.to_csv(os.path.join(os.getcwd(),'output','chamber_raw.csv'))


print('#### RAW CHAMBER GEOMETRY EXPORTED - NOW PERFORMING REGEN ANALYSIS ON MODIFIED CHAMBER FILE ####')
# seed heat transfer analysis
mu = 1.0145*0.0001 # millipoise -> Pa.s
Cp = 4.17*10**3 # J/kg.K (need to convert)
Pr = 0.55 # avg from CEA output with lower bias
cstar = 1724.3 # m/s - CEA output
RC_throat = 0.015 # radius of curvature at throat - spoof 15mm for now

cond_w = 370 # copper, W/m.K
C = 0.026 # Bartz constant

# compute 4x bartz constants which do not vary over nozzle
b1 = (C/((D_t_mm/1000)**0.2)) # first constant in bartz correlation C/Dt**0.2 - diameter correlation
b2 = ((mu**0.2)*Cp)/(Pr**0.6) # 2nd constant - mu**0.2.Cp/Pr**0.6 - transport props
b3 = (Pc_pa/cstar)**0.8 # 3rd constant - Pc correlation
b4 = ((D_t_mm/1000)/RC_throat)**0.1 # 4th constant - throat curvature correction
bartz_mult = b1*b2*b3*b4

chamber = import_regen_geo(chamber) # merge CSV of channel geometry
print(chamber)
#print(chamber)
chamber = chamber_thermo_calcs(chamber,k,R_specific,T_c,Pc)
chamber['A_chan'] = chamber['w_chan'] * chamber['d_chan'] # mm**2
chamber['D_hyd'] = 2*chamber['A_chan']/(chamber['w_chan'] + chamber['d_chan']) # mm
chamber['mdot_chan'] = mdot_fuel_ideal/chamber['n_chan']

T_eth_0 = 300 # K - assume roomish temp 70F
i_n = len(chamber.index) - 1  # grab i of nth index
# seed coolant temp params, and then set starting coolant temp
for parm in ['T_c_i','dT_c','T_c_e']:
    chamber[parm] = np.zeros(len(chamber.index))
chamber.at[i_n,'T_c_i'] = T_eth_0

T_wg_guess = 700
T_wc_guess = T_wg_guess - 100
q_guess = 1.0E06
HT_SOLVE_GUESS = [T_wg_guess, T_wc_guess, q_guess]

### Do heat transfer ###
for i in range(len(chamber.index))[::-1]:
    # NOTE: iterate thru in reverse
    if i > 0:
        # stop at 0th pt, will march upstream & use i-1th pt for area dx calc

        #grab transpo props for Bartz
        M = chamber.at[i,'mach']
        eps = chamber.at[i,'eps']
        AR_factor = (1 / eps) ** 0.9

        # Calculate HT geometry
        chamber.at[i,'dx'] = chamber.at[i,'x'] - chamber.at[i-1,'x'] # mm
        # Slant Area: https://www.calculatorsoup.com/calculators/geometry-solids/conicalfrustum.php
        # S = π * (r1 + r2) * s = π * (r1 + r2) * √((r1 - r2)2 + h2)
        chamber.at[i, 'A_w_segment'] = np.pi * (chamber.at[i,'r']+chamber.at[i - 1,'r']) * np.sqrt(((chamber.at[i,'r']-chamber.at[i-1,'r'])**2)+chamber.at[i,'dx']**2) # mm**2
        # Calc Adiabatic Wall Temp
        chamber.at[i,'T_aw'] = t_adiabatic_wall(T_c,Pr,M,k)

        #Get coolant inlet properties:
        rho_c, cp_c, cond_c, visc_c = get_ethanol_props_SI(chamber.at[i,'T_c_i'],45.0)
        chamber.at[i,'rho_c'] = rho_c
        chamber.at[i,'cp_c'] = cp_c
        chamber.at[i,'cond_c'] = cond_c
        chamber.at[i,'visc_c'] = visc_c
        chamber.at[i,'u_c'] = chamber.at[i,'mdot_chan']/(chamber.at[i,'rho_c']*(chamber.at[i,'A_chan']/(1000**2))) # Coolant Velo - m/s
        chamber.at[i,'Re_c'] = (chamber.at[i,'rho_c']*(chamber.at[i,'D_hyd']/1000)*chamber.at[i,'u_c'])/chamber.at[i,'visc_c']  # Coolant Re; convert D_hyd to m
        chamber.at[i,'Pr_c'] = (chamber.at[i,'cp_c']*chamber.at[i,'visc_c'])/chamber.at[i,'cond_c']
        chamber.at[i,'f_darcy'] = (0.79*np.log(chamber.at[i,'Re_c']) - 1.64)**(-2) # Petukhov correlation; Incorpera pg.490
        chamber.at[i,'Nu_c'] = gnielinski_calc(chamber.at[i,'f_darcy'], chamber.at[i,'Re_c'], chamber.at[i,'Pr_c'])
        chamber.at[i,'h_c'] = chamber.at[i, 'Nu_c'] * chamber.at[i, 'cond_c'] / (chamber.at[i, 'D_hyd']/1000) # Coolant h_c; convert D_hyd to mm

        # Run 1D heat transfer solver using Numpy Root Find
        T_c_MOD = 1.0
        data = (chamber.at[i,'t_wall']/1000,
                cond_w,
                bartz_mult*AR_factor,
                T_c/T_c_MOD,
                chamber.at[i, 'T_c_i'],
                chamber.at[i, 'h_c'],
                chamber.at[i, 'T_aw'],
                M,
                k) #create args for solver


        qmod = 2.0
        root = fsolve(func_bartz_root_find,HT_SOLVE_GUESS,args=data)
        chamber.at[i, 'T_wg'] = root[0]
        chamber.at[i, 'T_wc'] = root[1]
        chamber.at[i, 'q_tot'] = root[2]/qmod
        chamber.at[i,'dT_c'] = (chamber.at[i,'q_tot']*chamber.at[i,'A_w_segment']/(1000.0**2))/(mdot_fuel_ideal*chamber.at[i,'cp_c'])
        chamber.at[i, 'T_c_e'] = chamber.at[i, 'T_c_i'] + chamber.at[i, 'dT_c']
        chamber.at[i - 1, 'T_c_i'] = chamber.at[i, 'T_c_e']

        print(root)



#plot_chamber_param(chamber,'q','W/m**2')
#plot_chamber_param(chamber,'T_wc','K')
#plot_chamber_param(chamber,'w_fin','mm')
plot_chamber_thermo(chamber)

print(chamber)
chamber.to_csv(os.path.join(os.getcwd(),'output','chamber_final.csv'))


if FULL_OUTPUT:
    print(combustor_std.get_full_cea_output(Pc=Pc, MR=MR,
                                            pc_units='bar',
                                            PcOvPe=1/pressure_ratio))


