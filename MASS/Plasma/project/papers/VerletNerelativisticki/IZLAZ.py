##########################################################
# Priprema za 2D graficko predstavljanje - koordinatne ose
##########################################################

def set_plot_axes():
	
	global sp1, sp2
	sp1.set(adjustable='box', aspect='equal')
	sp2.set(adjustable='box', aspect='equal')
	sp1.set_xlabel('X')
	sp1.set_ylabel('Y')
	sp2.set_xlabel('X')
	sp2.set_ylabel('Z')
	sp1.set_xlim([-10, 10]) #korigovati po potrebi
	sp1.set_ylim([-10, 10])
	sp2.set_xlim([-10, 10])
	sp2.set_ylim([-10, 10])

##################################################################
# Priprema za 2D graficko predstavljanje - inicijalizacija grafika
##################################################################

def init_plots(plt):
	
	global sp1, sp2
	fig, (sp1, sp2) = plt.subplots(1, 2, figsize = (20, 10))
	set_plot_axes()
	return fig, sp1, sp2

###############################################################
# Priprema za 2D graficko predstavljanje - prikazivanje grafika
###############################################################

def write_plots(sp1, sp2, plt, x, y, z, t, T_sim):
	
	d = 1e8 #korigovati po potrebi, u zavisnosti od rezultata
	s = "T = " + "%7.1f" % float(t) + " of" + "%7.1f" % float(T_sim) + " [s]"
	sp1.set_title(s)
	sp1.plot(x/d, y/d, 'b.', markersize=2)
	sp2.plot(x/d, z/d, 'r.', markersize=2)
	plt.pause(0.00001)
