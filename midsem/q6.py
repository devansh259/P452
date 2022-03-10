import numpy as np

import library as *
# open a and b
with open('A.txt') as f:
    a = []
    for i in range(0, 6):
        a.append(list(map(float, f.readline().split())))
a = np.array(a)
opn = open('B.txt', 'r')
lsplit = opn.readline().split()
b = []
for val in lsplit:
    b.append(float(val))
b = np.array(b)



# jacobieq

vec2 = lib.jacobieq(a, b,25, 1e-5)
print('\nJacobi : x = {}'.format(vec2))

#OUTPUT

#Jacobi: x=[1.5 -0.5 2 -2.5 1 -1]