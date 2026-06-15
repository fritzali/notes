###########################
# Zadavanje pocetnih uslova
###########################

def puslovi(N, np):
	
	import random
	
	########################################
	# Ulazni podaci (fizicke karakteristike)
	########################################

	m_e = 9.10938356e-31 # Masa elektrona [kg]
	m_p = 1.6726219e-27 # Masa protona [kg]
	qelem = 1.6021766210e-19 # Elementarno naelektrisanje [C]
	c = 299792458.0 # Brzina svetlosti u vakuumu [m/s]
	vini = 0.01*c # Ogranicenje za pocetnu brzinu - pocetne brzine su u intervalu (-vini, vini)

	q = np.zeros(N)
	m = np.zeros(N)
	x = np.zeros(N)
	y = np.zeros(N)
	z = np.zeros(N)
	vx = np.zeros(N)
	vy = np.zeros(N)
	vz = np.zeros(N)
	
	for i in range(0,N):

		q[i] = qelem # Kretanje jednog protona se posmatra u ovom primeru
		m[i] = m_p 
		x[i] = 0.0 
		y[i] = 0.0 
		z[i] = 0.0
		vx[i] = 0.0
		vy[i] = (random.random() - 0.5)*vini
		vz[i] = (random.random() - 0.5)*vini	

	return q, m, x, y, z, vx, vy, vz
