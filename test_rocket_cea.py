import numpy as np
from rocketcea.cea_obj import CEA_Obj
import pandas as pd
import matplotlib.pyplot as plt

fac_CR = 4.0
eps = 5

C = CEA_Obj( oxName='LOX', fuelName='ETHANOL', fac_CR=fac_CR)
mr = np.linspace(1.0,3.5,num=20) #select 10 MR_sweep points
pc = np.linspace(250,500,num=20) #select 10 chamber pressure points
CEA_out_df = pd.DataFrame(index=mr)

for mr in mr:
    CEA_out_df.at[mr,'TAdiabatic'] = C.get_Tcomb(Pc=500, MR=mr)
    CEA_out_df.at[mr,'ISP_sl'] = 5 #
    x = C.estimate_Ambient_Isp(Pc=500, MR=mr, eps=eps)
    print(x)
    print(type(x))
    #print(mr, C.estimate_Ambient_Isp(Pc=500.0, MR_sweep=mr, eps=5))

print(CEA_out_df)
plt.plot(CEA_out_df)
plt.show()