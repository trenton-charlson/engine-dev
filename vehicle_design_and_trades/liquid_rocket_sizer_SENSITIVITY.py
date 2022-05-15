"""
Size out rocket design based on high level input parameters

future: perform trajectory sim & iterative analysis
"""

# std packages
import copy
import itertools
import matplotlib.pyplot as plt
from general_tools import matplotlib_config
matplotlib_config._matplotlib_init_() ## Init w/ custom font cfgs etc
import os
outputdir = os.path.join(os.getcwd(),'outputs')

# custom packages
import pandas as pd
from feed_system_design_code import gas_dyn_utils
from vehicle_sizing_functions import mass_tank_segment
from engine_design_code import engine_sizing_run
from _1D_rocket_traj import _1D_rocket_traj
from sizing_constants import *

## CTRL PARAMS ##
PLOT_PRESSURANT_SENSITIVITY = False
PLOT_VEHICLE = False
PLOT_TRAJECTORY = False
PRINT_RESULTS = False

## TOP LEVEL PARAMS ##
#PC = 21.0 # bars
P_exit = 1.01325 #bars
#thrust = 3000 # newton
#MR = 1.7 # mixture ratio
ffc_pct = 0.15 # % film coolant mass flow
oxidizer = engine_sizing_run.propellant('LOX')
fuel = engine_sizing_run.propellant('Kerosene')
fac_CR = 8.0 # face contraction ratio

DRY_MASS_GROWTH_FACTOR = 1.30  # spoof for future mass growth

#burntime = 18.0 # seconds
P_p_BOL = np.round(3000/BAR2PSI,2) # beginning bottle pressure
PRESSGASS = 'N2'

############################################################################################################
""" SWEEP PARAMS IN FOR LOOP"""
thrust_sweep = np.linspace(2000,4000,num=5)
bt_sweep = np.linspace(6.0,20.0,num=2)
MR_sweep = [1.7,1.9,2.1]
PC_sweep = [21.0]
i = 0 #index
out = pd.DataFrame()
alt_targ = 5000
alt_margin = 1000

fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2,ncols=2,figsize=(16,9))

print(f'\n###################################\n'
      f'Simulating {len(thrust_sweep)*len(bt_sweep)*len(PC_sweep)*len(MR_sweep)} Trajectory Cases\n'
      f'###################################\n'
      )

colors = matplotlib_config.OKABE_ITO  # std color pallette
for thrust in thrust_sweep:
    for MR in MR_sweep:
        for PC in PC_sweep:
            for burntime in bt_sweep:
                # Run high level engine sizer to extract flowrates
                eng, T_c, T_t, R_specific, k, opt_expansion, v_e_ideal = \
                        engine_sizing_run.size_combustor(PC, MR, thrust, fac_CR, oxidizer, fuel, P_exit/PC)

                mdot_ideal_total = thrust/v_e_ideal
                mdot_o = (MR/(1+MR))*mdot_ideal_total
                mdot_fuel_ideal = mdot_ideal_total - mdot_o
                mdot_fuel_regen = mdot_fuel_ideal*(1+ffc_pct)
                mdot_f = copy.copy(mdot_fuel_regen)

                q_ox = mdot_o/rho_o # m**3/s
                V_ox_i = q_ox*burntime
                m_o_i = mdot_o*burntime

                q_f = mdot_f/rho_f # m**3/s
                V_f_i = q_f*burntime
                m_f_i = mdot_f*burntime

                # Engine/Feed Pressures
                P_f_inj = (1+(f_inj_stiff/100))*PC
                P_f_inlet = (1+(f_reg_stiff/100))*P_f_inj
                P_o_inj = (1+(o_inj_stiff/100))*PC
                P_o_inlet = copy.copy(P_o_inj)

                # Tankage Pressures
                P_o_tank = P_o_inlet+lineloss_o
                P_f_tank = P_f_inlet+lineloss_f

                # Pressurant (rest in constants)
                P_p_EOL = np.round(max(P_f_tank,P_o_tank) + press_margin)

                L_o_tank = (V_ox_i*ullage_frac_o)/tank_AInternal
                m_o_tank = mass_tank_segment(L_o_tank,tank_LW,skin_LW)
                L_f_tank = (V_f_i*ullage_frac_f)/tank_AInternal
                m_f_tank = mass_tank_segment(L_f_tank,tank_LW,skin_LW)

                ## SIZE PRESSURANT SYSTEM ##
                # TODO: move this upstream, no need to re-sim blowdown every single time
                """
                Use blowdown sim to calculate time req given starting pressure target
                """
                blowdown = gas_dyn_utils.blowdown_sensitivity_study(vol_sweep,[P_p_BOL],T_p_LOAD,P_p_EOL,
                                                                    P_o_tank, q_ox,
                                                                    P_f_tank, q_f,
                                                                    burntime,
                                                                    PLOT=PLOT_PRESSURANT_SENSITIVITY)

                V_p_BOL = blowdown.at[np.round(P_p_BOL,2),f'{PRESSGASS} Vol Required']
                m_p_i = blowdown.at[np.round(P_p_BOL,2),f'{PRESSGASS} Mass, Initial']
                m_p_f = blowdown.at[np.round(P_p_BOL,2),f'{PRESSGASS} Mass Residual']
                mdot_p = (m_p_i-m_p_f)/burntime # average pressurant flowrate

                m_propellant_i = m_f_i + m_o_i
                m_struct = (m_thrust_struct+m_engine+m_o_tank+m_f_tank+m_valves+m_nc+m_aft_vp)*DRY_MASS_GROWTH_FACTOR

                m_wet = m_propellant_i+m_struct+m_p_i
                m_dry = m_struct+m_p_i  # assume worst case all pressurant stays onboard

                mdot_t = mdot_f+mdot_o  # total propellant flowrate

                ## Trajectory Sim & Analysis ##
                traj = _1D_rocket_traj(m_wet, m_dry, mdot_t, thrust, A_cs,
                                       PLOT=PLOT_TRAJECTORY)
                alt_max = np.round(max(traj['x']))
                vel_max = np.round(max(traj['v']),2)

                ## PRINT RESULTS ##
                if PRINT_RESULTS:
                    print(f'\n##########################################################\n'
                          f'COMBUSTOR MASS FLOWS:'
                          f'\n##########################################################\n\n'
                          f'>> mdot_ideal_total = {mdot_ideal_total} kg/s\n'
                          f'>> mdot_f (ideal) = {mdot_fuel_ideal} kg/s\n'
                          f'>> mdot_f (total, w/ film cooling) = {mdot_fuel_regen} kg/s\n'
                          f'>> mdot_o (ideal) = {mdot_o} kg/s')

                    print(f'\n##########################################################\n'
                          f'VEHICLE SIZING INPUTS:'
                          f'\n##########################################################\n\n'
                          f'>> Thrust = {thrust} [N] - ({np.round(thrust*NEWTON2LBF)} [lbf])\n'
                          f'>> Chamber Pressure = {PC} [bar] - ({np.round(PC*BAR2PSI)} [psia]) -- MR = {MR}\n'
                          f'>> Burn Time = {burntime} [s]\n'
                          f'>> Pressurant: {PRESSGASS}\n'
                          f'>> Vehicle Diameter: {skin_OD} [mm] - ({np.round(skin_OD/25.4,3)} [in])\n'
                          f'>> DRY MASS MULT = {DRY_MASS_GROWTH_FACTOR} [-]')

                    print(f'\n##########################################################\n'
                          f'VEHICLE SIZING OUTPUTS:'
                          f'\n##########################################################\n\n'
                          f'>> Propellant Mass = {m_propellant_i} [kg]\n'
                          f'>> Pressurant Mass = {m_p_i} [kg]\n'
                          f'>> Structural Mass = {np.round(m_struct,2)} [kg]\n'
                          f'>> TOTAL LIFTOFF MASS = {np.round(m_wet,2)} [kg]\n'
                          f'>> TOTAL DRY MASS = {np.round(m_dry,2)} [kg]\n'
                          f'>> PROPELLANT MASS FRACTION = {np.round(m_propellant_i/m_wet,2)}\n'
                          f'>> LIFTOFF TWR = {np.round(thrust/(9.81*m_wet),2)} [-]\n'
                          f'\n'
                          f'>> ALTITUDE ACHIEVED = {alt_max} [m]\n'
                          f'>> MAX VELOCITY = {vel_max} [m/s]\n')

                if PLOT_VEHICLE:
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

                # print sensitivity study outputs
                print(f'Case {i} >> {thrust} [N], {burntime} [s], {PC} bar, {MR} [-] -- Max Altitude = {np.round(alt_max,2)}')
                out.at[i, 't'] = thrust
                out.at[i, 'bt'] = burntime
                out.at[i, 'mr'] = MR
                out.at[i, 'pc'] = PC
                out.at[i, 'h'] = alt_max
                out.at[i, 'm_wet'] = m_wet
                out.at[i, 'm_prop'] = m_propellant_i
                out.at[i, 'V_p_req'] = V_p_BOL
                i = i+1


# Now repeat iteration loop to plot; make easier to read
colors = copy.copy(matplotlib_config.OKABE_ITO)  # std color pallette
for thrust in thrust_sweep:
    color = next(colors)
    linestyles = itertools.cycle(('-',':','--'))
    for MR in MR_sweep:
        # reset plot attribute for each thrust
        markers = matplotlib_config.MARKERS_STD
        linestyle = next(linestyles)
        for PC in PC_sweep:
            marker = ','
            # Plot const thrust curve for study output
            # filter out x & y vectors for plot calls
            bt_vec = out[(out['t'] == thrust) & (out['mr'] == MR) & (out['pc'] == PC)]['bt']
            h_vec = out[(out['t'] == thrust) & (out['mr'] == MR) & (out['pc'] == PC)]['h']
            mwet_vec = out[(out['t'] == thrust) & (out['mr'] == MR) & (out['pc'] == PC)]['m_wet']
            vp_vec = out[(out['t'] == thrust) & (out['mr'] == MR) & (out['pc'] == PC)]['V_p_req']

            ax1.plot(bt_vec, h_vec, color=color, marker=marker, linestyle=linestyle,
                     label=f'F = {thrust} [N]; PC = {PC} [BarA]; MR = {MR} [-]')
            ax2.plot(bt_vec, mwet_vec, color=color, marker=marker, linestyle=linestyle,
                     label=f'F = {thrust} [N]; PC = {PC} [BarA]; MR = {MR} [-]')
            ax3.plot(bt_vec, vp_vec, color=color, marker=marker, linestyle=linestyle,
                     label=f'F = {thrust} [N]; PC = {PC} [BarA]; MR = {MR} [-]')

            for burntime in bt_sweep:
                # nothing special to see here, carry on - for sweeping BT in pressurant sims
                continue
for MR in MR_sweep:
    for PC in PC_sweep:
        case_filt = out[(out['mr'] == MR) & (out['pc'] == PC)] # grab burntimes for each thrust that share common PC
        bt_req = np.round([np.interp((alt_targ+alt_margin),case_filt[case_filt['t'] == thrust]['h'],case_filt[case_filt['t'] == thrust]['bt']) for thrust in thrust_sweep],2)
        ax4.plot(thrust_sweep,bt_req,label=f'Thrust vs. Req Burn Time - PC: {PC} [BarA], MR: {MR} [-]')

#ax4.plot(thrust_sweep,bt_req,label='Thrust vs. Burn Time')
ax41 = ax4.twinx()
#ax41.plot(thrust_sweep,vp_req,label='Thrust vs. Pressurant Volume')
#ax41.plot(thrust_sweep,mprop_req,label='Thrust vs. Prop Mass Req.')
## FORMAT PLOT ##
ax1.axhline(y=alt_targ, c='k', ls='--', label="Altitude Targ")
ax1.axhline(y=alt_targ+alt_margin, c='r', ls='--', label=f'Altitude Targ + {alt_margin} m.o.s.')
ax1.legend()
ax1.grid()
ax2.legend()
ax2.grid()
ax3.legend()
ax3.grid()
ax4.legend()
ax41.legend()
ax4.grid()
ax3.set_xlabel('burn_time [s]')
ax4.set_xlabel('Thrust [N]')
ax4.set_ylabel('Burn Time [s]')
ax41.set_ylabel('Volume, Mass - [L],[kg]')
ax1.set_ylabel('alt [m]')
ax2.set_ylabel('wet mass [kg]')
ax3.set_ylabel('pressurant volume [L]')
fig1.suptitle(f'Burn Time vs Thrust Trade - MR={MR}\n'
              f'{thrust_sweep[0]}N - {thrust_sweep[-1]}N Thrust\n'
              f'Req Burn Times for h={alt_targ+alt_margin} => {0} [s]\n'
              f'Tanks: {tank_MATL}\n'
              f'Airframe: {skin_MATL}\n'
              f'Pressurant: {PRESSGASS}')
fig1.tight_layout()
output_fn = f'sensitivity_study_{thrust_sweep[0]}N-{thrust_sweep[-1]}N_{bt_sweep[0]}s-{bt_sweep[-1]}s_{tank_MATL}_{skin_MATL}.png'
plt.savefig(os.path.join(outputdir,output_fn))
