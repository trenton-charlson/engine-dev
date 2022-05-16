"""
Engine Design Code Top Level Script

author: tcharlson
"""

# Public Modules
import copy
import os
import numpy as np
from scipy.optimize import fsolve
from tqdm import tqdm

# Custom Modules
from geometry import chamber_geo, define_contour, plot_chamber_param
from nozzle_thermo import chamber_thermo_calcs, t_adiabatic_wall, gnielinski_calc, plot_chamber_thermo, \
                          func_bartz_root_find, load_regen_geometry, soot_thermal_resistance
from kerosene_props import get_jet_a_properties
from engine_sizing_run import size_combustor, propellant


### DEFINE TOP LEVEL PARAMETERS ###
fuel = propellant("Kerosene")
ox = propellant("LOX")

### PERFORMANCE PARAMS ###
Pc = 21.0 #bar
Pc_pa = Pc * 10**5 # chamber pressure in Pascals
Pe = 1.01325 #bar
thrust = 3000 #newton
# 2250 -> ~500lbf
# 4500 -> ~1000 lbf
MR = 1.7
eta_cstar = 0.9 # guess
pressure_ratio = Pe/Pc
print(f'pressure_ratio = {pressure_ratio}')
print(f'Pc/Pe = {Pc/Pe}')

ffc_pct = 0.15 # % FFC Mass fraction

### GEOMETRIC PARAMS ###
fac_CR = 7.0 #contraction ratio
chamber_ha = 30.0 #degrees, chamber half angle
nozzle_ha = 15.0 #degrees, nozzle half angle; conical nozzle
L_star = 1.15 #meters - http://mae-nas.eng.usu.edu/MAE_5540_Web/propulsion_systems/section6/section.6.1.pdf guess for now
npts = 100 #number of points to define each contour section
CHANNEL_CONFIG = 'regen_cfg.yaml' # config file name for regen parameters


f_inj_stiff = 20.0 # percent
P_f_inj = (1+(f_inj_stiff/100))*Pc

BAR = True
DEBUG = False

eng, T_c, T_t, R_specific, k, opt_expansion, v_e_ideal = size_combustor(Pc,MR,thrust,fac_CR,ox,fuel,pressure_ratio)

cstar_ideal = np.sqrt(k*R_specific*T_c)/(k*np.sqrt((2/(k+1))**((k+1)/(k-1))))
cstar_corr = eta_cstar*cstar_ideal
print(f'C*, ideal = {cstar_ideal} m/s - eta_C* = {eta_cstar}\n'
      f'C*, corr = {cstar_corr} m/s')

print('\n### COMBUSTOR MASS FLOWS ###')
mdot_ideal_total = thrust/v_e_ideal
print(f'mdot_ideal_total = {mdot_ideal_total} kg/s\n')

mdot_ox_ideal = (MR/(1+MR))*mdot_ideal_total
mdot_fuel_ideal = mdot_ideal_total - mdot_ox_ideal
mdot_fuel_regen = mdot_fuel_ideal*(1+ffc_pct)
print(f'mdot_f (ideal) = {mdot_fuel_ideal} kg/s')
print(f'mdot_f (regen) = {mdot_fuel_regen} kg/s')
print(f'mdot_f_ffc = {mdot_fuel_ideal*ffc_pct} kg/s')
print(f'mdot_o (ideal) = {mdot_ox_ideal} kg/s')

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

print(f'sanic {eng.combustor_obj.get_SonicVelocities(Pc=Pc,MR=MR,eps=opt_expansion)}')

print('\n#### CALCULATING CHAMBER GEOMETRIC PROPERTIES ####')
chamber_obj = chamber_geo(A_t_mm, R_t_mm, A_c_mm, R_c_mm, A_e_mm, R_e_mm)
chamber_raw,i_t = define_contour(chamber_obj, L_star, nozzle_ha, chamber_ha, PLOT=False)
chamber = copy.deepcopy(chamber_raw) # save raw copy of chamber data for export
chamber_raw.to_csv(os.path.join(os.getcwd(), 'output', 'chamber_raw.csv'))

print('\n#### RAW CHAMBER GEOMETRY EXPORTED ####')
# seed heat transfer analysis
# At the moment, am not super confident in the CEA output for transport props, does not seem to match CEA web unlike P/T props
# Running case in CEA Web App and manually inputting the below
mu = 0.91*0.0001 # millipoise -> Pa.s
Cp = 3.25*10**3 # J/kg.K (need to convert)
Pr = 0.50 # avg from CEA output with lower bias
#mod with eta cstar
RC_throat = 0.015 # radius of curvature at throat - spoof 15mm for now

cond_w = 350 # copper, W/m.K - conservatively low bound
C = 0.026 # Bartz constant

# compute 4x bartz constants which do not vary over nozzle
b1 = (C/((D_t_mm/1000)**0.2)) # first constant in bartz correlation C/Dt**0.2 - diameter correlation
b2 = ((mu**0.2)*Cp)/(Pr**0.6) # 2nd constant - mu**0.2.Cp/Pr**0.6 - transport props
b3 = (Pc_pa/(cstar_ideal))**0.8 # 3rd constant - Pc correlation - use ideal C* to be conservative
b4 = ((D_t_mm/1000)/RC_throat)**0.1 # 4th constant - throat curvature correction
bartz_mult = b1*b2*b3*b4

chamber = chamber_thermo_calcs(chamber, k, R_specific, T_c, Pc) # populate flow props vs station - isentropic
chamber = load_regen_geometry(chamber, CHANNEL_CONFIG)

chamber['A_chan'] = chamber['w_chan'] * chamber['d_chan'] # mm**2
chamber['D_hyd'] = 2*chamber['A_chan']/(chamber['w_chan'] + chamber['d_chan']) # mm
chamber['mdot_chan'] = mdot_fuel_regen/chamber['n_chan']

T_eth_0 = 300 # K - assume roomish temp 70F
i_n = len(chamber.index) - 1  # grab i of nth index
# seed coolant temp params, and then set starting coolant temp
for parm in ['T_c_i','dT_c','T_c_e','P_c_i','dP_c','P_c_e']:
    chamber[parm] = np.zeros(len(chamber.index))
chamber.at[i_n,'T_c_i'] = T_eth_0

T_wg_guess = 700
T_wc_guess = T_wg_guess - 100
q_guess = 1.0E06
HT_SOLVE_GUESS = [T_wg_guess, T_wc_guess, q_guess]

### Do heat transfer ###
for i in tqdm(range(len(chamber.index))[::-1]):
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
        chamber.at[i,'s'] = np.sqrt(((chamber.at[i,'r']-chamber.at[i-1,'r'])**2)+chamber.at[i,'dx']**2) # mm - slant height/linear dist travelled
        chamber.at[i, 'A_w_segment'] = np.pi * (chamber.at[i,'r']+chamber.at[i - 1,'r']) * chamber.at[i,'s'] # mm**2
        # Calc Adiabatic Wall Temp
        chamber.at[i,'T_aw'] = t_adiabatic_wall(T_c,Pr,M,k)
        chamber.at[i,'R_soot'] = soot_thermal_resistance(chamber.at[i,'eps'],chamber.at[i,'regime']) # (m**2 K)/W - Huzel & Huang

        #Get coolant inlet properties:
        rho_c, cp_c, cond_c, visc_c = get_jet_a_properties(chamber.at[i,'T_c_i'])
        chamber.at[i,'rho_c'] = rho_c
        chamber.at[i,'cp_c'] = cp_c
        chamber.at[i,'cond_c'] = cond_c
        chamber.at[i,'visc_c'] = visc_c
        chamber.at[i,'u_c'] = chamber.at[i,'mdot_chan']/(chamber.at[i,'rho_c']*(chamber.at[i,'A_chan']/(1000**2))) # Coolant Velo - m/s
        chamber.at[i,'t_transit'] = chamber.at[i,'s']/(chamber.at[i,'u_c']*1000)
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
                k,
                chamber.at[i, 'R_soot']) #create args for solver

        qmod = 1.0
        root = fsolve(func_bartz_root_find,HT_SOLVE_GUESS,args=data)
        chamber.at[i, 'T_wg'] = root[0]
        chamber.at[i, 'T_wc'] = root[1]
        chamber.at[i, 'q_tot'] = root[2]/qmod
        chamber.at[i, 'dT_c'] = (chamber.at[i,'q_tot']*chamber.at[i,'A_w_segment']/(1000.0**2))/(mdot_fuel_ideal*chamber.at[i,'cp_c'])
        chamber.at[i, 'T_c_e'] = chamber.at[i, 'T_c_i'] + chamber.at[i, 'dT_c']
        chamber.at[i - 1, 'T_c_i'] = chamber.at[i, 'T_c_e']
        # calculate darcy-weisbach pressure drop
        # https://en.wikipedia.org/wiki/Darcy%E2%80%93Weisbach_equation
        chamber.at[i, 'dP_c'] = (chamber.at[i,'s']/1000 * (chamber.at[i,'f_darcy']*(chamber.at[i,'rho_c']/2)*((chamber.at[i,'u_c']**2)/(chamber.at[i,'D_hyd']/1000))))/10**5 # dP in Bar

        if DEBUG:
            print(root)

    else:
        # populate final station with non-NAN values
        # i = 0
        keys = [key for key in chamber.keys() if key not in ['x', 'r', 'theta', 'eps', 'regime']]
        for key in keys:
            chamber.at[0, key] = chamber.at[1, key]

# Invert Coolant Pressure Vectors to match flow direction
chamber.at[0,'P_c_e'] = P_f_inj
for i in range(len(chamber.index)-1):
    chamber.at[i,'P_c_i'] = chamber.at[i,'P_c_e'] + chamber.at[i,'dP_c']
    chamber.at[i+1,'P_c_e'] = chamber.at[i,'P_c_i']

# Plot thermo outputs
fig = plot_chamber_thermo(chamber,eng)
fig.tight_layout()
fig.show()

chamber_OD_w_cooling = np.round(np.max(chamber["r_outer"]),2)*2

# Console Output
print(f'Chamber OD w/Cooling = {chamber_OD_w_cooling} [mm] ({chamber_OD_w_cooling/25.4} [in])\n')
print(f'Regen Inlet Pressure = {chamber.at[i_n,"P_c_e"]} [BarA]')
print(f'Injector Inlet Pressure = {chamber.at[0,"P_c_e"]} [BarA]\n')
print(f'Regen Stiffness = {100*(chamber.at[i_n,"P_c_e"] - chamber.at[0,"P_c_e"])/chamber.at[0,"P_c_e"]} [%]\n')
print(f'Injector Inlet Fuel Props:')
print(f'T_c = {chamber.at[0,"T_c_e"]} [K]')
print(f'rho_c = {chamber.at[0,"rho_c"]} [kg/s]')
print(f'visc_c = {chamber.at[0,"visc_c"]} [Pa.s]')

# Save final chamber output case
chamber.to_csv(os.path.join(os.getcwd(),'output','chamber_final.csv'))


