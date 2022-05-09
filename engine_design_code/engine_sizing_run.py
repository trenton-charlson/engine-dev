"""
Top level function(s) for defining engine sizing & outputs
"""

from rocketcea.cea_obj import CEA_Obj as CEA_Obj_english
from rocketcea.cea_obj_w_units import CEA_Obj
import numpy as np

import engine_design_code.constants

R_const = engine_design_code.constants.R_const
g = engine_design_code.constants.g

def dump():
    combustor_std = CEA_Obj_english(oxName=ox.name, fuelName=fuel.name,
                                    fac_CR=fac_CR)
    return None


class propellant:
    def __init__(self,name):
        self.name = name


class engine:
    def __init__(self,thrust,Pc,MR,fac_CR,combustor_obj,pressure_ratio,PRINT=False):
        self.thrust = thrust
        self.Pc = Pc
        self.MR = MR
        self.fac_CR = fac_CR
        self.combustor_obj = combustor_obj
        self.pressure_ratio = pressure_ratio
        self.PRINT = PRINT

        self.optimum_expansion = self.get_opt_expansion()
    def get_opt_expansion(self):
        (mw_c, k_c) = self.combustor_obj.get_Chamber_MolWt_gamma(Pc=self.Pc, MR=self.MR, eps=40)
        (mw_t, k_t) = self.combustor_obj.get_Throat_MolWt_gamma(Pc=self.Pc, MR=self.MR, eps=40)

        self.k = (k_c + k_t) / 2  # average gas constant
        self.mw = (mw_c + mw_t) / 2  # avg molecular weight

        opt_expansion = 1 / ((((self.k + 1) / 2) ** (1 / (self.k - 1))) * (self.pressure_ratio ** (1 / self.k)) * np.sqrt(
            ((self.k + 1) / (self.k - 1)) * (1 - (self.pressure_ratio ** ((self.k - 1) / self.k)))))

        if self.PRINT:
            print(f'CHAMBER - mw: {mw_c}, k: {k_c}')
            print(f'THROAT - mw: {mw_t}, k: {k_t}')
            print(f'K_avg = {self.k}')
            print(f'Optimum Expansion: {opt_expansion}')


        return opt_expansion


def size_combustor(Pc, MR, thrust, fac_CR, ox, fuel, pressure_ratio,
                   PRINT=False):

    combustor_obj = CEA_Obj(oxName=ox.name, fuelName=fuel.name,
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



    eng = engine(thrust,Pc,MR,fac_CR,combustor_obj,pressure_ratio)

    T_c = eng.combustor_obj.get_Tcomb(Pc=Pc, MR=MR)
    T_t = 2 * T_c / (eng.k + 1)
    R_specific = R_const / eng.mw
    v_e_ideal = np.sqrt(((2 * eng.k) / (eng.k - 1)) * R_specific * T_c * (1 - (eng.pressure_ratio ** ((eng.k - 1) / eng.k))))

    if PRINT:
        print(f'Target Thrust: {thrust} N,\n'
              f'Pc = {Pc} Bar\n'
              f'MR = {MR} [-]\n'
              f'fac_CR = {fac_CR} [-]\n')
        print(f'T_comb = {T_c} K')
        print(f'T_t = {T_t} K')
        print(f'R_specific = {R_specific} J/kg.K')
        print(f'V_e, ideal = {v_e_ideal} m/s')

    return eng, T_c, T_t, R_specific, eng.k, eng.optimum_expansion, v_e_ideal