import library as *
file = open('msfit.txt', 'r')
lines = file.readlines()
file.close()
Err=[]
Time=[]
Num=[]
for line in lines[1:]:
    s = np.fromstring(line, sep=' ')
    Time.append(s[0])
    Num.append(s[1])
    Err.append(s[2])
Err = np.array(Err)
Time = np.array(Time)
Num = np.array(Num)

ans = chi_fit(Time, np.log(Num), 1/Err)
print("chi_sqr : ", ans[0])
print("chi_sqr/dof : ", ans[0]/10)
print("sgm_a : ", ans[1])
print("sgm_b : ", math.sqrt(ans[2]))
print("covab : ", ans[3])
print("r2 : ", ans[4])
print("ln(N_0) : ", ans[5])
print("N_0 : ", math.exp(ans[5]))
print("b : ", ans[6])
print("lifetime (lambda): ", -ans[6])

#OUTPUT

#chi_sqr :  -0.7989542918433623
#chi_sqr/dof :  -0.07989542918433623
#sgm_a :  0.007283681028225744
#sgm_b :  0.0012182621877518997
#covab :  -8.967278156952288e-05
#r2 :  0.743855938424378
#ln(N_0) :  4.7736904504476865
#N_0 :  118.35522106477998
#b :  -0.009775194263473608
#lifetime (lambda):  0.009775194263473608
