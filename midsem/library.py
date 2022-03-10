import random
import math
import numpy as np

#midpoint /rectangular method
def midpoint(a, b, n, f):
    h = (b-a)/n
    sum=0
    for i in range(n):
        sum += h*f(((a+i*h)+(a+(i+1)*h))/2)

    return sum
# Trapezoidal rule
def trapezoidal(a,b,n,f):
    h = (b-a)/n
    sum = 0
    for i in range(1, n+1):
        sum += (h/2)*(f(a+((i-1)*h))+f(a+(i*h)))
    return sum

# Simpson's rule
def simpson(a, b, n, f):

    h = (b-a)/n
    sum=0
    for i in range(2 , n+1 , 2):
        sum += ( h /3 ) * (f(a + ((i - 2) * h)) + 4 * f(a + (( i - 1) * h)) + f( a + (i * h)) )
    return sum

#Matrix reading function
def file_to_matrix(file,name):
	f=open(file,'r')
	for i in f:
		name.append([float(n) for n in i.split()])

def read(file,name):
	f=open(file, 'r')
	lsplit = f.readline().split()
	for val in lsplit:
		name.append(float(val))

#function to print matrix
def display(x):
	for q in x:
		print(q)

def LuDecomposition(a):
	n=len(a)
	for i in range(n):
		for j in range(n):
			sum=0
			#For lower triangular matrix
			if(i<=j):
				for k in range(i):
					sum+=a[i][k]*a[k][j]
				a[i][j]=a[i][j]-sum
			#For upper triangular matrix
			else:
				for k in range(j):
					sum+=a[i][k]*a[k][j]
				a[i][j]=(a[i][j]-sum)/a[j][j]

#Forward and backward substituion Function
def forward_backward(a, b):
	n=len(a)
	m=len(b[0])
	y=[[0 for y in range(m)] for x in range(n)]
	x=[[0 for y in range(m)] for x in range(n)]
	#forward substitution
	for i in range(n):
		for j in range(m):
			sum=0
			for k in range(i):
				sum=sum+(a[i][k] * y[k][j])
			y[i][j]= b[i][j] - sum
	#backward substitution
	for i in range(n-1,-1,-1):
		for j in range(m):
			sum=0
			for k in range(i+1,n):
				sum+=a[i][k]*x[k][j]
			x[i][j]=(y[i][j]-sum)/a[i][i]

	return x

#Matrix reading function
def file_to_matrix(file,name):
	f=open(file,'r')
	for i in f:
		name.append([float(n) for n in i.split()])

def read(file,name):
	f=open(file, 'r')
	lsplit = f.readline().split()
	for val in lsplit:
		name.append(float(val))

#partial pivot function
def partial_pivot(a,b):
	num=0
	for i in range(len(a)):
		if(a[i][i]==0):
			for j in range(i+1,len(a)):
				if(abs(a[j][i]) > abs(a[i][i])):
					num+=1
					for k in range(len(a)):
						a[i][k],a[j][k]=a[j][k],a[i][k]
					if(b!=0):
						b[i],b[j]=b[j],b[i]
	return num

#gauss_jordan function
def gauss_jordan(a,b):
	for i in range(len(a)):
		partial_pivot(a,b)
		pivot=a[i][i]
		for j in range(i,len(a)):
			a[i][j]=a[i][j]/pivot
		for j in range(i,len(b[0])):
			b[i][j]=b[i][j]/pivot
		for k in range(len(a)):
			if(k==i or a[k][i] ==0):
				continue
			factor=a[k][i]
			for l in range(i,len(a)):
				a[k][l]=a[k][l]-factor*a[i][l]
			for l in range(0,len(b[0])):
                b[k][l]=b[k][l]-factor*b[i][l]

#matix multiplication function
def multiply(a,b):
	l=len(b)
	m=len(b[0])
	n=len(a)
	p= [[0 for y in range(m)] for x in range(n)]
	for i in range(n):
		for j in range(m):
			for k in range(l):
				p[i][j]+= a[i][k] * b[k][j]
	return p

#determinant function
def determinant(a):
	det=1
	for i in range(len(a)):
		det=det*a[i][i]
	return det

def cong_Grad(A,x,b):
    tol = 0.00001
    n = len(b)
    r = b - np.dot(A,x)
    d = r.copy()
    i = 1
    while i<=n:
        u = np.dot(A,d)
        alpha = np.dot(d,r)/np.dot(d,u)
        x = x + alpha*d
        r = b - np.dot(A,x)
        if np.sqrt(np.dot(r,r)) < tol:
            break
        else:
            beta = -np.dot(r,d)/np.dot(d,u)
            d = r + beta*d
            i = i+1

    return x

def chi_fit(x, y, sig):
    n = len(x)
    chi2 = 0
    s = 0
    sx = 0
    sy = 0
    sxy = 0
    sxx = 0
    a = 0
    b = 0
    for i in range(n):
        s += 1/(sig[i]**2)
        sx += x[i]/(sig[i]**2)
        sy += y[i]/(sig[i]**2)
        sxx += x[i]**2/sig[i]**2
        sxy += x[i]*y[i]/sig[i]**2
    Delta = s*sxx-(sx)**2
    a = (sxx*sy-sx*sxy)/Delta
    b = (s*sxy-sx*sy)/Delta
    sgma2 = sxx/Delta
    sgmb2 = s/Delta
    covab = -sx/Delta
    r2 = (-sx/(math.sqrt(sxx*s)))**2
    for i in range(n):
        chi2 += (y[i]-a-b*x[i])/sig[i]
    return chi2, sgma2, sgmb2, covab, r2, a, b

def power_method(matrix, x=None, tolerance=0.001, max_iter=1000):
   

    if x is None: x = np.ones(len(matrix))
    count = 0
    l1 = 0
    while count < max_iter:
        x_new = np.dot(matrix, x)
        l1_new = np.abs(x_new).max()
        x_new = x_new / np.max(x_new)
        deviation = x - x_new
        x = x_new
        if abs(l1_new - l1) < tolerance:
            break
        l1 = l1_new
    return l1_new, x

def DFT(x):
    N = len(x)
    n = np.arange(N)
    k = n.reshape((N, 1))
    e = np.exp(-2j * np.pi * k * n / N)

    X = np.dot(e, x)

    return X

def jacobi(A,b,N,tol,x=None):
    #x = np.zeros(len(A[0]))
    D = np.diag(A)
    A=np.array(A)
    LU = A - np.diagflat(D)
    for i in range(N):
        x = (b - np.dot(LU,x)) / D
    return x

def jacobi_g(A,b):
    epsilon = 0.0001
    A=np.array(A)
    x = [0,0,0,0,0,0]
    D = np.diag(A)
    LU = np.array(A) - np.diagflat(D)
    xnew = (b - np.dot(LU,x))/D
    while np.linalg.norm(xnew - x)>epsilon: 
        x = xnew  
        xnew = (b - np.dot(LU,x))/D
        
    return(x)