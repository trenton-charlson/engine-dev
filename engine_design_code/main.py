"""
Engine Design Code Top Level Script

author: tcharlson
"""

import numpy as np
import pandas as pd
from rocketcea.cea_obj import CEA_Obj as CEA_Obj_english
from rocketcea.cea_obj_w_units import CEA_Obj

from geometry import chamber_geo, define_contour
from nozzle_thermo import chamber_thermo_calcs

#import gas_dynamics as gd

### DEFINE TOP LEVEL PARAMETERS ###
class propellant:
    def __init__(self,name):
        self.name = name

fuel = propellant("Ethanol")
ox = propellant("LOX")
### PERFORMANCE PARAMS ###
Pc = 35.0 #bar
Pe = 1.01325 #bar
thrust = 2500 #newton
MR = 1.6
eta_cstar = 0.9 # guess
R_const = 8.314462*10**3 #J/K.mol
pressure_ratio = Pe/Pc
print(f'pressure_ratio = {pressure_ratio}')
print(f'Pc/Pe = {Pc/Pe}')

### GEOMETRIC PARAMS ###
fac_CR = 9.0 #contraction ratio
chamber_ha = 30.0 #degrees, chamber half angle
nozzle_ha = 15.0 #degrees, nozzle half angle; conical nozzle
L_star = 1.5 #meters - http://mae-nas.eng.usu.edu/MAE_5540_Web/propulsion_systems/section6/section.6.1.pdf guess for now

FULL_OUTPUT = True


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

mdot_ideal = thrust/v_e_ideal
print(f'mdot_ideal = {mdot_ideal} kg/s\n')

A_t = (mdot_ideal/(Pc*10**5))*np.sqrt((R_specific*T_c)/(k*((2/(k+1))**((k+1)/(k-1)))))
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
chamber,i_t = define_contour(chamber_obj, L_star, nozzle_ha, chamber_ha)
print(chamber.index)
chamber = chamber_thermo_calcs(chamber,k,R_specific,T_c,Pc)

test=5
test=6

if FULL_OUTPUT:
    print(combustor_std.get_full_cea_output(Pc=Pc, MR=MR,
                                            pc_units='bar',
                                            PcOvPe=1/pressure_ratio))



