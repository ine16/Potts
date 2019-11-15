import matplotlib.pyplot as plt
b,e,cv,m,chi = [],[],[],[],[]
for line in open('medidas_q2_L20.txt','r'):
	values = [float(s) for s in line.split()]
	b.append(values[0])
	e.append(values[1])
	cv.append(values[2])
	m.append(values[3])
	chi.append(values[4])

#N=50*50
# import numpy as np
# #b=1./np.array(b)
# e=np.array(e)/N
# cv=np.array(cv)/N
# chi=np.array(chi)/N
# m=np.array(m)/N

plt.figure(1)
plt.xlabel('T')
plt.ylabel('Energia - q=2')
plt.scatter(b,e, color='red')
plt.figure(2)
plt.xlabel('T')
plt.ylabel('Magnetizacion - q=2')
plt.scatter(b,m, color='black')
fig3=plt.figure(3)

plt.xlabel('T')
plt.ylabel('Calor especifico - q=2')
plt.scatter(b,cv, color='green')
fig4=plt.figure(4)
plt.xlabel('T')
plt.ylabel('Susceptibilidad - q=2', color='blue')
plt.scatter(b,chi)

plt.show()