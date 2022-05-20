from library import *
import numpy as np

rod_length = 2
dx = 0.1
x = np.linspace(0, rod_length, int(rod_length/dx))
u = 20*np.abs(np.sin(np.pi*x))
u[0] = 0
u[-1] = 0

check_points = [0, 10, 20, 50, 100, 200, 500]
count=0
for i in check_points:
    plt.plot(x, PDE_solver (0.0008,0.1,1,u,steps=i-count))
    count=i
    print(i)
plt.title("Temperature vs x for t=0")
plt.xlabel("x")
plt.ylabel("Temp")
