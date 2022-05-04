from library import *

def func(x):
    return math.sqrt(1-x**2)
  
N=10000
count_points=0
count_monte=0

#Pi by throwing the points Area of circle
def pi_throws(a,m,N,seed1,seed2):
    X = mLCG(a, 0, m,seed1 , N)
    Y = mLCG(a, 0, m, seed2, N)
    count=0
    for i in range(N):
        if(X[i]**2 +Y[i]**2) <=1:
            count+=1
    return 4*count*(float(1/N))
  
#Pi by monte carlo integration
def pi_monte(a,m,N,func,seed):
    X = mLCG(a1, 0, m1,0.9 , N)
    count=0
    for i in range(N):
        count+=func(X[i])
    return 4*count/N
  
print("Value of pi by random throw for a=65,m=1021 is: ",pi_throws(65,1021,10**4,0.9,1.7))
print("Value of pi by random monte carlo integration for a=65,m=1021 is: ",pi_monte(65,1021,10**4,func,1.7))
print("Value of pi by random throw for a=572,m=16381 is: ",pi_throws(572,16381,10**4,0.9,1.7))
print("Value of pi by monte carlo integration for a=572,m=16381 is: ",pi_monte(572,16381,10**4,func,1.7))

#OUTPUT

#Value of pi by random throw for a=65,m=1021 is:  3.1416
#Value of pi by random monte carlo integration for a=65,m=1021 is:  3.142850043457293
#Value of pi by random throw for a=572,m=16381 is:  3.128
#Value of pi by monte carlo integration for a=572,m=16381 is:  3.142850043457293
