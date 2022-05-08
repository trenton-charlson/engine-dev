"""
1D Trajectory Analysis Tools
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from ambiance import Atmosphere

def _1D_rocket_traj(mwet,mdry,mdot,thrust,A_cs,
                    cd=0.75,rho_air=1.225,
                    ts=0.1,
                    PLOT=True):

    heights = np.linspace(0.01,20000.0,num=100)
    atmosphere = Atmosphere(heights)
    rho = atmosphere.density

    # Compute 1D trajectory of rocket given baseline input parameters
    mfuel = mwet-mdry #fuel mass
    burn_time = mfuel/mdot #total burn time in seconds

    t = 0.1 #seed time vector wtih 1st timestep
    i = 1 #seed iteration vector
    traj = pd.DataFrame(columns=['mass', 'a', 'v', 'x'])
    traj.at[0, 'mass'] = mwet #seed initial mass
    traj.at[0, 'a'] = 0
    traj.at[0, 'v'] = 0
    traj.at[0, 'x'] = 0

    while t<burn_time:
        # rocket is thrusting
        vehicle_mass = mwet-(mdot*t)
        vehicle_weight = vehicle_mass*9.81 #weight in Newtons

        v_init = traj.iat[i-1, 2] #velocity at pos2
        rho_loc = np.interp(traj.iat[i-1,3],heights,rho)
        f_drag = 0.5 * rho_loc * v_init**2 * A_cs
        f_accel = thrust - vehicle_weight - f_drag
        accel = f_accel / vehicle_mass

        v_final = v_init + accel * ts
        dx = (v_final ** 2 - v_init ** 2) / (2 * accel)

        t = round(t,1)

        traj.at[t, 'mass'] = vehicle_mass
        traj.at[t, 'a'] = accel
        traj.at[t, 'v'] = v_final
        traj.at[t, 'x'] = traj.iat[i-1, 3] + dx
        traj.at[t, 'rho'] = rho_loc

        i = i+1
        t = t+ts

    while v_final > 0:
        # rocket is coasting
        vehicle_weight = vehicle_mass * 9.81  # weight in Newtons
        v_init = traj.iat[i - 1, 2]  # velocity at pos2
        rho_loc = np.interp(traj.iat[i-1,3],heights,rho)
        f_drag = 0.5 * rho_loc * v_init ** 2 * A_cs
        f_accel = vehicle_weight + f_drag
        accel = -f_accel / vehicle_mass #accel is nevgative

        v_final = v_init + accel * ts
        dx = (v_final ** 2 - v_init ** 2) / (2 * accel)

        t = round(t, 1)

        traj.at[t, 'mass'] = vehicle_mass
        traj.at[t, 'a'] = accel
        traj.at[t, 'v'] = v_final
        traj.at[t, 'x'] = traj.iat[i - 1, 3] + dx
        traj.at[t, 'rho'] = rho_loc

        i = i+1
        t = t+ts

    if PLOT:
        fig, (ax0, ax1, ax2, ax3) = plt.subplots(4, 1, figsize=(10,6), sharex=True)
        ax0.plot(traj['x'], label='Altitude')
        ax0.set_ylabel('Altitude [m]')
        ax1.plot(traj['v'], label='Velocity')
        ax1.set_ylabel('Velocity [m/s]')
        ax2.plot(traj['a'], label='Acceleration')
        ax2.set_ylabel('Acceleration [m/s**2]')
        ax3.plot(traj['rho'], label='Atmospheric Density')
        ax3.set_ylabel('rho [kg/m**3]')
        ax3.set_xlabel('Time [s]')


        fig.suptitle(f'1D Trajectory Sim\n'
                     f'Max Alt: {np.round(np.max(traj["x"]))} [m]\n'
                     f'Thrust: {thrust} [N], burntime: {burn_time} [s]')
        plt.show()

    return traj
