h = 'hello'
for letter in h:
    print(letter)

from rocketcea.cea_obj import CEA_Obj

ispObj = CEA_Obj( oxName='LOX', fuelName='Ethanol')


s = ispObj.get_full_cea_output( Pc=1000.0, MR=6.0, eps=40.0, short_output=1)

print( s )