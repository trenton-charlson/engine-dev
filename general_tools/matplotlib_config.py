"""
Initialize Matplotlib with new default color cycles AND font sizes
"""
import matplotlib.pyplot as plt
import itertools

OKABE_ITO = itertools.cycle(('#e69f00','#56B5E9','#009e73',"#D55E00",'#0072B2','#cc79a7','#000000'))
MARKERS_STD = itertools.cycle((',', '+', '.', 'o', '*'))
LINESTYLE_STD = itertools.cycle(('-','--',':','-.'))

def _matplotlib_init_():

    SMALL_SIZE = 8
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    return

