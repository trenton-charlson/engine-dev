"""
Functions for calculating thermodynamic properties of ethanol

author: tcharlson
"""

from CoolProp.CoolProp import PropsSI

def get_ethanol_props_SI(T,P,
                         BAR=True):
    if BAR:
        P = P*10**5 # convert to Pa for coolprop input

    rho = PropsSI('D','T',T,'P',P,'Ethanol')
    h = PropsSI('H','T',T,'P',P,'Ethanol')
    cp = PropsSI('CPMASS','T',T,'P',P,'Ethanol')
    cond = PropsSI('CONDUCTIVITY','T',T,'P',P,'Ethanol')
    visc = PropsSI('VISCOSITY','T',T,'P',P,'Ethanol')

    return rho, h, cp, cond, visc