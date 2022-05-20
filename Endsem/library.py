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
    t = len(a[1])
    
    with open('test.txt') as f:
        t = len(a[1])
        u = []
        for i in range(0, t):
            u.append(list(map(float, f.readline().split())))
    with open('test.txt') as f:
        l = []
        for i in range(0, t):
            l.append(list(map(float, f.readline().split())))
    n = len(u)
    
    for i in range(0, n):
        for j in range(0, n):
            sum = 0
            for k in range(0, i):
                sum += l[i][k]*u[k][j]
            u[i][j] = a[i][j]-sum
            
        for j in range(0, n):
            if i == j:
                l[i][i] = 1
            else:
                sum = 0
                for k in range(0, i):
                    sum += l[j][k]*u[k][i]
                l[j][i] = (a[j][i]-sum)/u[i][i]
    return l, u

# def lu_decomposition(A,b):
#     partial_pivot(A,b)
#     u= [[0 for i in range(len(A))] for j in range(len(A))]
#     l= [[0 for i in range(len(A))] for j in range(len(A))]
#     for i in range(0, n):
#         for j in range(0, n):
#             sum = 0
#             for k in range(0, i):
#                 sum += l[i][k]*u[k][j]
#             u[i][j] = a[i][j]-sum
            
#         for j in range(0, n):
#             if i == j:
#                 l[i][i] = 1
#             else:
#                 sum = 0
#                 for k in range(0, i):
#                     sum += l[j][k]*u[k][i]
#                 l[j][i] = (a[j][i]-sum)/u[i][i]

#     return(forward_)


#Forward and backward substituion Function
def forward_backward(a, b):
    l = LuDecomposition(a)[0]
    u = LuDecomposition(a)[1]
    n = len(u)
    z = []
    for i in range(0, n):
        z.append(0)
    for i in range(0, n):
        if i == 0:
            z[i] = b[i]/l[i][i]
        else:
            sum = 0
            for j in range(0, i):
                sum += l[i][j]*z[j]
            z[i] = (b[i]-sum)/l[i][i]
    x = []
    for i in range(0, n):
        x.append(0)
    i = n-2
    x[n-1] = z[n-1]/u[n-1][n-1]
    while i >= 0:
        sum = 0
        for j in range(i, n):
            sum += u[i][j]*x[j]
        x[i] = (z[i]-sum)/u[i][i]
        i = i-1

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
def gauss_jordan(a, b):
    partial_pivot(a, b)
    r = len(b)
    for k in range(0, r):
        pivot = a[k][k]
        # for the pivot row divide by akk
        for j in range(k, r):
            a[k][j] = a[k][j]/pivot
        b[k] = b[k]/pivot
        # for the non pivot row
        for i in range(0, r):
            if i == k or a[i][k] == 0:
                continue
            factor = a[i][k]
            # substraction
            for d in range(k, r):
                a[i][d] = a[i][d]-factor * a[k][d]
            b[i] = b[i]-factor * b[k]
    return a, b

    

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

def jacobi(A,b,N,x=None):
    if x is None:
        x = np.zeros(len(A[0]))

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

def gauss_seidel(A,b,tol):
    max_iter=10000
    x = np.zeros_like(b, dtype=np.double)
    er = np.array([])
    kk = np.array([])
    # Iterate
    for k in range(max_iter):
        x_old = x.copy()

        # Loop over rows
        for i in range(A.shape[0]):
            x[i] = (b[i] - np.dot(A[i, :i], x[:i]) - np.dot(A[i, (i+1):], x_old[(i+1):])) / A[i, i]
        er = np.append(er, np.linalg.norm(x - x_old, ord=np.inf) / np.linalg.norm(x, ord=np.inf))
        kk = np.append(kk, k)
        # Stop condition
        if np.linalg.norm(x - x_old, ord=np.inf) / np.linalg.norm(x, ord=np.inf) < tolerance:
            break

    return x, er, kk
	
def mLCG(a, c, m, seed, N):
    x = [0 for i in range(N)] 
    x[0] = seed
    for i in range(N-1):
        x[i+1] = (a*x[i] + c)%m
    for i in range(N):
        x[i] = x[i]/m
    return x
def polynomial_fit(X,Y,degree=1):
   
    n = len(X)
    m=degree+1
    A = np.zeros((m,m))     
    b = np.zeros(m)               

    for i in range(m):
        for j in range(m):
            count = 0
            for k in range(n):
                count += X[k] ** (i+j)

            A[i,j] = count

    for i in range(m):
        count = 0
        for k in range(n):
            count += (X[k]**i * Y[k])

        b[i] = count

    return(forward_backward(A, b))
	
def cheby_polynomial(x,degree):
    if degree ==0: return 1
    elif degree == 1:return 2*x -1
    elif degree == 2:return 8*x**2 -8*x +1
    elif degree == 3:return 32*x**3 -48*x**2 + 18*x -1

def cheby_fit(X,Y,degree=1):
    n = len(X)
    m=degree+1
    A = np.zeros((m,m))     
    b = np.zeros(m)               

    for i in range(m):
        for j in range(m):
            count = 0
            for k in range(n):
                count += cheby_polynomial(X[k], j) * cheby_polynomial(X[k], i)

            A[i,j] = count

    for i in range(m):
        count = 0
        for k in range(n):
            count += cheby_polynomial(X[k],i) * Y[k]

        b[i] = count

    return(forward_backward(A, b))

def monteCarlo_int(func,N,a,b):
   
    x = mLCG(234.34, 65, 1,0.5, N)

    sum = 0
    for i in range(N):
        sum+=((b-a)/N)*f(x[i-1])

    total = 1/float(N) * sum

    return total

def LuDecomposition(a):
 
    n = len(a)
    u = np.empty((n,n))
    l = np.empty((n,n))
    
    for i in range(0, n):
        for j in range(0, n):
            sum = 0
            for k in range(0, i):
                sum += l[i][k]*u[k][j]
            u[i][j] = a[i][j]-sum
            
        for j in range(0, n):
            if i == j:
                l[i][i] = 1
            else:
                sum = 0
                for k in range(0, i):
                    sum += l[j][k]*u[k][i]
                l[j][i] = (a[j][i]-sum)/u[i][i]
    return l, u

def legendre_fit(X,Y,degree=1):
    n = len(X)
    m=degree+1
    A = np.zeros((m,m))     
    b = np.zeros(m)               

    for i in range(m):
        for j in range(m):
            count = 0
            for k in range(n):
                count += legendre_polynomial(X[k], j) * legendre_polynomial(X[k], i)

            A[i,j] = count

    for i in range(m):
        count = 0
        for k in range(n):
            count += cheby_polynomial(X[k],i) * Y[k]

        b[i] = count

    return(forward_backward(A, b))



    return x

def random_walk(step_num):
    walk_num = 500
    

    sum = []
    xx=[]
    yy=[]
    for i in range(walk_num):
        random_numbers = mLCG(N=step_num, a=572,c= 0,m=16381, seed=i)
        p = np.array(list(map(lambda x: math.floor(x*4), random_numbers)))
        x = 0
        y = 0
        for s in p:
            if s %2 == 0:
                if s == 0:
                    x += 1
                else:
                    x -= 1
            else:
                if s == 1:
                    y += 1
                else:
                    y -= 1
        xx.append(x)
        yy.append(y)
        plt.plot(x,y)
        sum.append(pow(x**2 + y**2, 0.5))

    sqr_sum = 0.
    for s in sum:
        sqr_sum += s**2
    return(pow(sqr_sum/walk_num, 0.5),xx,yy)

