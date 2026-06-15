def const_B(x, y, z, t):
	
	# Zadato je homogeno i stacionarno B polje sa komponentom Bz (duz z-ose)
	
	BPx = 0.0*1e-10
	BPy = 0.0*1e-10
	BPz = 3.0*1e-10 #U mikro Gausima [uG] se zadaje puat konverzija - 3uG u ovom primeru
	
	EPx = 0.0
	EPy = 0.0
	EPz = 0.0
	
	return BPx, BPy, BPz, EPx, EPy, EPz

def const_EB(x, y, z, t):

	# Homogeno i stacionarno magnetno i elektricno polje, 
	# sa komponentama Bz, Ex = Ey, Ez = 0
	
	BPx = 0.0*1e-10
	BPy = 0.0*1e-10
	BPz = 3.0*1e-10
	
	EPx = BPz*1e4
	EPy = BPz*1e4
	EPz = 0.0
		
	return BPx, BPy, BPz, EPx, EPy, EPz

def magdipol(x, y, z, t, np):

	sinphi = np.sin(11.7*np.pi/180.0)
	cosphi = np.cos(11.7*np.pi/180.0)
	pomka = ((x**2.0 + y**2.0 + z**2.0)**2.5)
	BPx = -7.965626e15*((3.0*x*z*cosphi) + (3.0*x*y*sinphi))/pomka
	BPy = -7.965626e15*((3.0*y*z*cosphi) + (2.0*sinphi*(y**2.0)) - (sinphi*(x**2.0)) - (sinphi*(z**2.0)))/pomka
	BPz = -7.965626e15*((2.0*cosphi*(z**2.0)) - (cosphi*(x**2.0)) - (cosphi*(y**2.0)) + (3.0*z*y*sinphi))/pomka
	
	EPx = 0.0
	EPy = 0.0
	EPz = 0.0
		
	return BPx, BPy, BPz, EPx, EPy, EPz
