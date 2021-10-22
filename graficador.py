import matplotlib.pyplot as plt

def graficador(q, L):
	"Levanta los datos del .txt correspondiente al Monte Carlo para Potts de q estados y red de LxL"
	t, e, cv, m, chi = [], [], [], [], [] # temperatura, energía, calor específico, magnetización y susceptibilidad
	for line in open('medidas_q'+str(q)+'_L'+str(L)+'.txt', 'r'):
		values = [float(s) for s in line.split()]
		t.append(values[0])
		e.append(values[1])
		cv.append(values[2])
		m.append(values[3])
		chi.append(values[4])

	plt.figure(1)
	plt.xlabel('T')
	plt.ylabel(f'Energia - q= {q}')
	plt.scatter(t, e, color='red')
	plt.figure(2)
	plt.xlabel('T')
	plt.ylabel(f'Magnetizacion - q= {q}')
	plt.scatter(t, m, color='black')
	fig3=plt.figure(3)

	plt.xlabel('T')
	plt.ylabel(f'Calor especifico - q= {q}')
	plt.scatter(t, cv, color='green')
	fig4=plt.figure(4)
	plt.xlabel('T')
	plt.ylabel(f'Susceptibilidad - q= {q}')
	plt.scatter(t, chi, color='blue')

	plt.show()
