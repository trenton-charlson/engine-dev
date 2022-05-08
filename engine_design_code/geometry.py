"""
Functions for geometric modelling of combustor

author: tcharlson
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os


def define_contour(chamber,L_star,nozl_ha,chamber_ha,
                   rnd=3,npts=100,PLOT=True):
    """

    :param chamber:
    :param L_star:
    :param nozl_ha:
    :param chamber_ha:
    :param step:
    :return:
    """

    conv_v,conv_l = chamber.get_frustrum_params(chamber_ha,chamber.Rc)
    nozl_v,nozl_l = chamber.get_frustrum_params(nozl_ha,chamber.Re)
    chamber_v = chamber.At*L_star*1000
    barrel_v = chamber_v-conv_v
    barrel_l = barrel_v/chamber.Ac

    total_l = barrel_l+conv_l+nozl_l
    print(f'Total Combustor Length = {total_l} mm')

    x_coords = np.unique([np.linspace(0.0,barrel_l,num=npts),
                                          np.linspace(barrel_l,barrel_l+conv_l,num=npts),
                                          np.linspace(barrel_l+conv_l,total_l,num=npts)])
    #print(x_coords)

    #print(x_coords)
    contour = pd.DataFrame(index=np.arange(len(x_coords)))
    contour['x'] = x_coords
    contour['r'] = np.zeros(len(x_coords))
    contour['theta'] = np.zeros(len(x_coords))

    print(f'Barrel Length = {barrel_l} mm')
    print(f'Converging Length = {conv_l} mm')
    print(f'Nozzle Length = {nozl_l} mm')

    for i in range(len(x_coords)):
        if i<npts: #chamber barrel
            contour.at[i,'r'] = chamber.Rc

        elif (i>=npts-1 and i<(npts*2)-2):
            # converging section
            dx = x_coords[i+1] - x_coords[i]
            dr = dx*np.tan(np.deg2rad(chamber_ha))
            contour.at[i,'r'] = contour.at[i-1,'r']-dr
            contour.at[i,'theta'] = chamber_ha # wall angle
        elif i == npts*2-2:
            #we at the throat boiiiii
            contour.at[i,'r'] = chamber.Rt
            contour.at[i,'theta'] = 0.0
        elif (i>(npts*2)-2 and i<(npts*3)-3):
            # nozzle segment
            dx = x_coords[i+1] - x_coords[i]
            dr = dx*np.tan(np.deg2rad(nozl_ha))
            contour.at[i,'r'] = contour.at[i-1,'r']+dr
            contour.at[i,'theta'] = nozl_ha # wall angle
        elif i == (npts*3)-3:
            #at exit
            contour.at[i,'r'] = chamber.Re
            contour.at[i,'theta'] = nozl_ha

    contour['eps'] = [np.round(np.pi*(contour.at[i,'r']**2)/chamber.At,rnd) for i in range(len(x_coords))]
    contour['r'] = np.round(contour['r'],rnd)
    contour['x'] = np.round(contour['x'],rnd)
    i_t = npts*2-2 #throat index
    contour['regime'] = np.zeros(len(x_coords))

    for i in range(len(x_coords)):
        if i>i_t:
            contour.at[i,'regime'] = 2

    contour.at[i_t, 'regime'] = 1  # set throat regime to 1.0

    if PLOT:
        fig1,(ax1,ax2) = plt.subplots(nrows=2,ncols=1,sharex=True)
        ax1.plot(contour['x'],contour['r'],c='k',lw=2)
        ax1.plot(contour['x'],-contour['r'],c='k',lw=2)
        ax1.plot(contour['x'],np.zeros(len(x_coords)),c='r',ls='--')
        ax1.set(ylim=(0.0,max(contour['r']*1.25)))
        ax1.axis('equal')
        ax2.set_xlabel('x [mm]')
        ax1.set_ylabel('r [mm]')
        ax2.set_ylabel('eps [-]')
        ax2.plot(contour['x'],contour['eps'])
        plt.show()

    return contour,i_t

def plot_chamber_param(chamber,param,param_units):
    fig1, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True)
    ax1.plot(chamber['x'],chamber['r'],c='k',lw=2)
    ax1.set(ylim=(0.0, max(chamber['r'] * 1.25)))

    ax2.plot(chamber['x'],chamber[param],c='r',lw=2)

    ax1.axis('equal')
    ax2.set_xlabel('x [mm]')
    ax1.set_ylabel('r [mm]')
    ax2.set_ylabel(f'{param}, [{param_units}]')

    plt.show()
    return

class chamber_geo:
    def __init__(self,At,Rt,Ac,Rc,Ae,Re):
        self.At = At
        self.Rt = Rt
        self.Ac = Ac
        self.Rc = Rc
        self.Ae = Ae
        self.Re = Re

    def get_frustrum_params(self,ha,R):
        # https: // www.calculatorsoup.com / calculators / geometry - solids / conicalfrustum.php  #:~:text=Volume%20of%20a%20conical%20frustum,r1%20*%20r2))
        self.h = (R-self.Rt)/np.tan(np.deg2rad(ha))
        #print(f'Converging Length = {h} mm')
        self.V_conv = (1.0/3.0) * np.pi * self.h * (self.Rt**2 + self.Rc**2 + (self.Rc*self.Rt))
        return self.V_conv,self.h




