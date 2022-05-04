from library import *
opn = open('q1B.txt', 'r')
lsplit = opn.readline().split()
b = []
for val in lsplit:
    b.append(float(val))
rows = 6
with open('q1A.txt') as f:
    a = []
    for i in range(0, rows):
        a.append(list(map(float, f.readline().split())))
        
gauss_jordan(a, b)
print('Gauss jordan solution :', b)
forward_backward(a, b)
print('LU decompostion solution:,b)

#OUTPUT
#Gauss jordan solution : [-1.7618170439978567, 0.8962280338740136, 4.051931404116157, -1.6171308025395428, 2.041913538501914, 0.15183248715593495]
#LU decompostion solution: [-1.7618170439978567, 0.8962280338740136, 4.051931404116157, -1.6171308025395428, 2.041913538501914, 0.15183248715593495]
