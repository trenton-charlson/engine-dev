"""
Functions for determining sizing parameters at vehicle level
"""

def mass_tank_segment(l,tank_lw,skin_lw,
                      beta_tank_str = 1.2,
                      beta_tank_endcaps = 1.1,
                      beta_skin_str = 1.2):

    # mass = len*[(structmult + endcapmult)*tank_linwt + skinmult*skin_linwt]
    mass = l*((1+(beta_tank_str-1)+(beta_tank_endcaps-1))*tank_lw + beta_skin_str*skin_lw)

    return mass

def add_vehicle_segment(fig,ax,
                        start_x,end_x,c,label):


    return None