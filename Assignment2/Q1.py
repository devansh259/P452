from library import *
from matplotlib import pyplot as plt

#Reading fitting data
data=open("assign2fit.txt", 'r')
line = [line.rstrip().split('\t') for line in data]
x=[row[0] for row in line]
y=[row[1] for row in line]
X=[float(x) for x in x]
Y=[float(y) for y in y]

#Normal fitting
var=polynomial_fit(X,Y,3)
print("coefficients by normal fitting:",var)

#Chebyshev fitting
var_cheby=cheby_fit(X,Y,3)
print("coefficients by chebyshev fitting:",var)

#ploting fitted data
x=np.linspace(0,1,100)
y=var[0] +var[1]*x +var[2]*x**2 +var[3]*x**3
plt.plot(x,y)
plt.scatter(X,Y,color='r')
plt.xlabel("x")
plt.ylabel("y")

plt.title("cubic polynomial fit")

#OUTPUT
#coefficients by normal fitting: [0.5746586674195995, 4.725861442142078, -11.128217777643616, 7.6686776229096685]
#coefficients by chebyshev fitting: [0.5746586674195995, 4.725861442142078, -11.128217777643616, 7.6686776229096685]
