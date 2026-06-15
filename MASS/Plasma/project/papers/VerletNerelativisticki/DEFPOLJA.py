def const_B(x, y, z, t):
	
	# Zadato je homogeno i stacionarno B polje sa komponentom Bz (duz z-ose)
	
	BPx = 0.0*1e-10
	BPy = 0.0*1e-10
	BPz = 3.0*1e-10 #U mikro Gausima [uG] se zadaje puta konverzija - 3uG u ovom primeru
	
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
