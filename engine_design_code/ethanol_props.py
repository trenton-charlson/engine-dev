"""
Functions for calculating thermodynamic properties of ethanol

author: tcharlson
"""
import CoolProp.CoolProp as CP

def get_ethanol_props_SI(T,P,
                         BAR=True,MIX=0.75):
    if BAR:
        P = P*10**5 # convert to Pa for coolprop input

    CP.apply_simple_mixing_rule('Ethanol','Water','linear')
    CP.set_config_bool(CP.OVERWRITE_BINARY_INTERACTION, True)
    FID = f'WATER[{1-MIX}]&ETHANOL[{MIX}]'

    rho = CP.PropsSI('D','T',T,'P',P,FID)
    #h = PropsSI('H','T|supercritical_liquid',T,'P',P,FID)
    cp = CP.PropsSI('CPMASS','T',T,'P',P,FID)
    cond = CP.PropsSI('CONDUCTIVITY','T',T,'P',P,FID)
    visc = CP.PropsSI('VISCOSITY','T',T,'P',P,FID)

    return rho, cp, cond, visc