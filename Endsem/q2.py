from library import *
import numpy as np
from matplotlib import pyplot as plt


def legendre_polynomial(x,degree):
    if degree ==0: return 1
    elif degree == 1:return x
    elif degree == 2:return (3*x**2 - 1)/2
    elif degree == 3:return (5*x**3 -3*x)/2
    elif degree == 4:return (35*x**4 -30*x**2 +3)/8
    elif degree == 5:return (63*x**5-70*x**3 +15*x)/8
    elif degree == 6:return (231*x**6 -315*x**4 +105*x**2 -5)/16
    
data=open("esem4fit.txt", 'r')
line = [line.rstrip().split('\t') for line in data]
x=[row[0] for row in line]
y=[row[1] for row in line]
X=[float(x) for x in x]
Y=[float(y) for y in y]

#Normal fitting
var=legendre_fit(X,Y,6)
print("coefficients by legendre fitting:",var)

#plotting fitted data
x=np.linspace(-1,1,100)
y=var[0] +var[1]*x +var[2]*x**2 +var[3]*x**3 + var[4]*x**4 +var[5]*x**5 + var[6]*x**6
plt.plot(x,y)
plt.scatter(X,Y,color='r')
plt.xlabel("x")
plt.ylabel("y")

plt.title("polynomial fit")


#OUTPUT
#coefficients by legendre fitting: 
#[0.1217799533691444, -0.027937002485179815, -0.5244683033757472, 0.09157036596125233, 0.7430896777491651, -0.05297490625865053, -0.17880208085143237]

    

    

