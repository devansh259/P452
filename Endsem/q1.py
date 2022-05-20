from library import *

x=random_walk(200)[1]
y=random_walk(200)[2]
plt.plot(x,y)
plt.xlim([-100,100])
plt.title("Random walk for n=200 , walks=500")

print( "R_rms for N=200, walks=500:",random_walk(200)[0])

sum_rms=[]
x=[]
for i in range(1,10):
    sum_rms.append(random_walk(i*100)[0])
    x.append((i*100)**0.5)

plt.scatter(x,sum_rms)
plt.xlim([0,40])
plt.ylim([0,50])
plt.xlabel('$\sqrt{N}$')
plt.ylabel('$R_{rms}$')

#OUTPUT
#R_rms for N=200, walks=500: 16.425102739404707


