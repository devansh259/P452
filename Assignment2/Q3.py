from library import *
from math import *

def func(x):
    return (2*sqrt(1-x**2))**2 
  
vol=monteCarlo_int(func,10**4,-1,1)
print("Volume of Steinmetz solid by monte carlo:",vol)

#OUTPUT

#Volume of Steinmetz solid by monte carlo: 5.3515763661269675
