# Python file to loosely trade total thrust vs vehicle sizing
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def _1D_rocket_traj(mwet,mdry,mdot,thrust,A_cs,
                    cd=0.75,rho_air=1.225,
                    ts=0.1):

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
        f_drag = 0.5 * rho_air * v_init**2 * A_cs
        f_accel = thrust - vehicle_weight - f_drag
        accel = f_accel / vehicle_mass

        v_final = v_init + accel * ts
        dx = (v_final ** 2 - v_init ** 2) / (2 * accel)

        t = round(t,1)

        traj.at[t, 'mass'] = vehicle_mass
        traj.at[t, 'a'] = accel
        traj.at[t, 'v'] = v_final
        traj.at[t, 'x'] = traj.iat[i-1, 3] + dx

        i = i+1
        t = t+ts

    while v_final > 0:
        # rocket is coasting
        vehicle_weight = vehicle_mass * 9.81  # weight in Newtons
        v_init = traj.iat[i - 1, 2]  # velocity at pos2
        f_drag = 0.5 * rho_air * v_init ** 2 * A_cs
        f_accel = vehicle_weight + f_drag
        accel = -f_accel / vehicle_mass #accel is nevgative

        v_final = v_init + accel * ts
        dx = (v_final ** 2 - v_init ** 2) / (2 * accel)

        t = round(t, 1)

        traj.at[t, 'mass'] = vehicle_mass
        traj.at[t, 'a'] = accel
        traj.at[t, 'v'] = v_final
        traj.at[t, 'x'] = traj.iat[i - 1, 3] + dx

        i = i+1
        t = t+ts

    return traj

mwet = 81.5#kg
mdry = 69.5 #kg
mdot = 1.22 #kg/s
thrust = 3000.0 #Newton
D = 6.17 # inches - 6" tube
A_cs_i = np.pi*(D/2)**2
A_cs = A_cs_i * 0.00064516 # in**2 -> m**2

print(A_cs)

traj = _1D_rocket_traj(mwet, mdry, mdot, thrust, A_cs)
print(f'Max Altitude: {np.max(traj["x"])}')

fig, (ax0, ax1, ax2) = plt.subplots(3,1,sharex=True)
ax0.plot(traj['x'], label='Altitude')
ax1.plot(traj['v'], label='Velocity')
ax2.plot(traj['a'], label='Acceleration')
plt.show()

#test master->main


#endline