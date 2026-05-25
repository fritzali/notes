### preamble ###

import scipy as sci
import pandas as pd
import numpy as np
import matplotlib as mpl

import matplotlib.pyplot as plt

### theory ###

def gff(E, kT):
	return np.exp(E / (2 * kT)) * sci.special.k0(E / (2 * kT))

def eps(E, kT):
	return gff(E, kT) * np.exp(- E / kT) / kT**(1/2)

E = np.logspace(-2.0, 1.0, 333)

eps1 = eps(E, 0.1)
eps2 = eps(E, 1.0)
eps3 = eps(E, 10.0)

plt.plot([], [], label=r'$kT = \qty{e-1}{\kilo\electronvolt}$', ls=(0, (3, 1)), c='k', lw=0.8)
plt.plot([], [], label=r'$kT = \qty{e+0}{\kilo\electronvolt}$', ls=(0, (8, 2)), c='k', lw=0.8)
plt.plot([], [], label=r'$kT = \qty{e+1}{\kilo\electronvolt}$', ls='-', c='k', lw=0.8)
plt.plot(E[eps1 > 0.01], eps1[eps1 > 0.01], ls=(0.8, (3, 1)), c='k', lw=0.8)
plt.plot(E[eps2 > 0.01], eps2[eps2 > 0.01], ls=(0, (8, 2)), c='k', lw=0.8)
plt.plot(E[eps3 > 0.01], eps3[eps3 > 0.01], ls='-', c='k', lw=0.8)

plt.xlabel(r'$E \mathbin{/} \unit{\kilo\electronvolt}$')
plt.ylabel(r'$\epsilon \mathbin{/} \epsilon_0$')

plt.xscale('log')
plt.yscale('log')

plt.legend()

plt.savefig('content/theo_bremsstrahlung.pdf')

plt.close()

def nfw(r, rs = 250.0, rhos = 1e6):
	return rhos / ((r / rs) * (1 + r / rs)**2)

def cc(r, rc = 150.0, rhog = 6.5e4, b = 0.65):
	return rhog * (1 + (r / rc)**2)**(-1.5 * b)

def ncc(r, rc1 = 25.0, rhog1 = 5.0e5, b1 = 0.50, rc2 = 200.0, rhog2 = 3.0e4, b2 = 0.65):
	inner = rhog1 * (1 + (r / rc1)**2)**(-1.5 * b1)
	outer = rhog2 * (1 + (r / rc2)**2)**(-1.5 * b2)
	return inner + outer

r = np.logspace(1.0, 3.3, 333)

plt.plot(r, nfw(r), label=r'Total Density', ls='-', c='k', lw=0.8)
plt.plot(r, ncc(r), label=r'Noncool Core', ls=(0, (8, 2)), c='k', lw=0.8)
plt.plot(r, cc(r), label=r'Cool Core', ls=(0, (3, 1)), c='k', lw=0.8)

plt.xlabel(r'$r \mathbin{/} \unit{\kilo\parsec}$')
plt.ylabel(r'$\rho \mathbin{/} \unit{\solarmass\per\kilo\parsec\cubed}$')

plt.xscale('log')
plt.yscale('log')

plt.legend()

plt.savefig('content/theo_halo.pdf')

plt.close()

### spectra ###

E85, F85, Ferr85 = np.genfromtxt('catalog/Chandra/acisf15173N003_evt2-A85.txt', unpack=True)
E426, F426, Ferr426 = np.genfromtxt('catalog/Chandra/acisf04952N004_evt2-A426.txt', unpack=True)
E1644, F1644, Ferr1644 = np.genfromtxt('catalog/Chandra/acisf07922N003_evt2-A1644.txt', unpack=True)
E1656, F1656, Ferr1656 = np.genfromtxt('catalog/Chandra/acisf13996N003_evt2-A1656.txt', unpack=True)
E2029, F2029, Ferr2029 = np.genfromtxt('catalog/Chandra/acisf04977N006_evt2-A2029.txt', unpack=True)
E3158, F3158, Ferr3158 = np.genfromtxt('catalog/Chandra/acisf03712N004_evt2-A3158.txt', unpack=True)

plt.plot(E1656, F1656, label='A1656 (Coma)', lw=0.6, alpha=0.6, c='chocolate', zorder=5)
plt.plot(E426, F426, label='A426 (Perseus)', lw=0.6, alpha=0.6, c='firebrick', zorder=6)

plt.plot([], [], ' ', label=' ')
plt.plot([], [], ' ', label=' ')

plt.plot(E2029, F2029, label='A2029', lw=0.6, alpha=0.6, c='goldenrod', zorder=4)
plt.plot(E85, F85, label='A85', lw=0.6, alpha=0.6, c='olivedrab', zorder=3)
plt.plot(E3158, F3158, label='A3158', lw=0.6, alpha=0.6, c='steelblue', zorder=2)
plt.plot(E1644, F1644, label='A1644', lw=0.6, alpha=0.6, c='rebeccapurple', zorder=1)

plt.xlabel(r'$E \mathbin{/} \unit{\kilo\electronvolt}$')
plt.ylabel(r'$N \mathbin{/} \unit{\per\second\per\kilo\electronvolt}$', labelpad=-0.5)

plt.yscale('log')

plt.legend(ncol=2)

plt.savefig('content/chan_spectra.pdf')

plt.close()

### profiles ###

R85, ne85, neerr85, K85, Kerr85, P85, Perr85, T85, Terr85 = np.genfromtxt('catalog/Chandra/acisf15173N003_evt2-ABELL_0085.dat', unpack=True)
R426, ne426, neerr426, K426, Kerr426, P426, Perr426, T426, Terr426 = np.genfromtxt('catalog/Chandra/acisf04952N004_evt2-ABELL_0426.dat', unpack=True)
R1644, ne1644, neerr1644, K1644, Kerr1644, P1644, Perr1644, T1644, Terr1644 = np.genfromtxt('catalog/Chandra/acisf07922N003_evt2-ABELL_1644.dat', unpack=True)
R2029, ne2029, neerr2029, K2029, Kerr2029, P2029, Perr2029, T2029, Terr2029 = np.genfromtxt('catalog/Chandra/acisf04977N006_evt2-ABELL_2029.dat', unpack=True)
R3158, ne3158, neerr3158, K3158, Kerr3158, P3158, Perr3158, T3158, Terr3158 = np.genfromtxt('catalog/Chandra/acisf03712N004_evt2-ABELL_3158.dat', unpack=True)

plt.plot(R426[1:], ne426[1:], label='A426', lw=0.9, alpha=0.75, c='firebrick', zorder=5)
plt.plot(R2029[1:], ne2029[1:], label='A2029', lw=0.9, alpha=0.75, c='goldenrod', zorder=4)
plt.plot(R85[1:], ne85[1:], label='A85', lw=0.9, alpha=0.75, c='olivedrab', zorder=3)
plt.plot(R3158[1:], ne3158[1:], label='A3158', lw=0.9, alpha=0.75, c='steelblue', zorder=2)
plt.plot(R1644[1:], ne1644[1:], label='A1644', lw=0.9, alpha=0.75, c='rebeccapurple', zorder=1)

plt.xlabel(r'$R \mathbin{/} \unit{\mega\parsec}$')
plt.ylabel(r'$n \mathbin{/} \unit{\per\centi\meter\cubed}$')

plt.xscale('log')
plt.yscale('log')

plt.legend()

plt.savefig('content/chan_density.pdf')

plt.close()

plt.plot(R426, T426, label='A426', lw=0.9, alpha=0.75, c='firebrick', zorder=5)
plt.plot(R2029, T2029, label='A2029', lw=0.9, alpha=0.75, c='goldenrod', zorder=4)
plt.plot(R85, T85, label='A85', lw=0.9, alpha=0.75, c='olivedrab', zorder=3)
plt.plot(R3158, T3158, label='A3158', lw=0.9, alpha=0.75, c='steelblue', zorder=2)
plt.plot(R1644, T1644, label='A1644', lw=0.9, alpha=0.75, c='rebeccapurple', zorder=1)

plt.xlabel(r'$R \mathbin{/} \unit{\mega\parsec}$')
plt.ylabel(r'$kT \mathbin{/} \unit{\kilo\electronvolt}$', labelpad=5.5)

plt.ylim(0.97 * plt.ylim()[0], 1.03 * plt.ylim()[1])

plt.xscale('log')

plt.legend()

plt.savefig('content/chan_temperature.pdf')

plt.close()

plt.plot(R426[1:], K426[1:], label='A426', lw=0.9, alpha=0.75, c='firebrick', zorder=5)
plt.plot(R2029[1:], K2029[1:], label='A2029', lw=0.9, alpha=0.75, c='goldenrod', zorder=4)
plt.plot(R85[1:], K85[1:], label='A85', lw=0.9, alpha=0.75, c='olivedrab', zorder=3)
plt.plot(R3158[1:], K3158[1:], label='A3158', lw=0.9, alpha=0.75, c='steelblue', zorder=2)
plt.plot(R1644[1:], K1644[1:], label='A1644', lw=0.9, alpha=0.75, c='rebeccapurple', zorder=1)

plt.xlabel(r'$R \mathbin{/} \unit{\mega\parsec}$')
plt.ylabel(r'$K \mathbin{/} \unit{\kilo\electronvolt\centi\meter\squared}$')

plt.xscale('log')
plt.yscale('log')

plt.legend()

plt.savefig('content/chan_entropy.pdf')

plt.close()

plt.plot(R426[1:], P426[1:], label='A426', lw=0.9, alpha=0.75, c='firebrick', zorder=5)
plt.plot(R2029[1:], P2029[1:], label='A2029', lw=0.9, alpha=0.75, c='goldenrod', zorder=4)
plt.plot(R85[1:], P85[1:], label='A85', lw=0.9, alpha=0.75, c='olivedrab', zorder=3)
plt.plot(R3158[1:], P3158[1:], label='A3158', lw=0.9, alpha=0.75, c='steelblue', zorder=2)
plt.plot(R1644[1:], P1644[1:], label='A1644', lw=0.9, alpha=0.75, c='rebeccapurple', zorder=1)

plt.xlabel(r'$R \mathbin{/} \unit{\mega\parsec}$')
plt.ylabel(r'$P \mathbin{/} \unit{\erg\per\centi\meter\cubed}$')

plt.xscale('log')
plt.yscale('log')

plt.legend()

plt.savefig('content/chan_pressure.pdf')

plt.close()

fig, axs = plt.subplots(2, 2, figsize=[7.81, 5.85])

axs[0, 0].plot(R426[1:], ne426[1:], label='A426', lw=0.9, alpha=0.75, c='firebrick', zorder=5)
axs[0, 0].plot(R2029[1:], ne2029[1:], label='A2029', lw=0.9, alpha=0.75, c='goldenrod', zorder=4)
axs[0, 0].plot(R85[1:], ne85[1:], label='A85', lw=0.9, alpha=0.75, c='olivedrab', zorder=3)
axs[0, 0].plot(R3158[1:], ne3158[1:], label='A3158', lw=0.9, alpha=0.75, c='steelblue', zorder=2)
axs[0, 0].plot(R1644[1:], ne1644[1:], label='A1644', lw=0.9, alpha=0.75, c='rebeccapurple', zorder=1)

axs[0, 0].set_xlabel(r'$R \mathbin{/} \unit{\mega\parsec}$')
axs[0, 0].set_ylabel(r'$n \mathbin{/} \unit{\per\centi\meter\cubed}$')

axs[0, 0].set_xscale('log')
axs[0, 0].set_yscale('log')

axs[0, 0].legend()

axs[1, 0].plot(R426, T426, label='A426', lw=0.9, alpha=0.75, c='firebrick', zorder=5)
axs[1, 0].plot(R2029, T2029, label='A2029', lw=0.9, alpha=0.75, c='goldenrod', zorder=4)
axs[1, 0].plot(R85, T85, label='A85', lw=0.9, alpha=0.75, c='olivedrab', zorder=3)
axs[1, 0].plot(R3158, T3158, label='A3158', lw=0.9, alpha=0.75, c='steelblue', zorder=2)
axs[1, 0].plot(R1644, T1644, label='A1644', lw=0.9, alpha=0.75, c='rebeccapurple', zorder=1)

axs[1, 0].set_xlabel(r'$R \mathbin{/} \unit{\mega\parsec}$')
axs[1, 0].set_ylabel(r'$kT \mathbin{/} \unit{\kilo\electronvolt}$', labelpad=5.5)

axs[1, 0].set_ylim(0.97 * axs[1, 0].set_ylim()[0], 1.03 * axs[1, 0].set_ylim()[1])

axs[1, 0].set_xscale('log')

axs[1, 0].legend()

axs[1, 1].plot(R426[1:], K426[1:], label='A426', lw=0.9, alpha=0.75, c='firebrick', zorder=5)
axs[1, 1].plot(R2029[1:], K2029[1:], label='A2029', lw=0.9, alpha=0.75, c='goldenrod', zorder=4)
axs[1, 1].plot(R85[1:], K85[1:], label='A85', lw=0.9, alpha=0.75, c='olivedrab', zorder=3)
axs[1, 1].plot(R3158[1:], K3158[1:], label='A3158', lw=0.9, alpha=0.75, c='steelblue', zorder=2)
axs[1, 1].plot(R1644[1:], K1644[1:], label='A1644', lw=0.9, alpha=0.75, c='rebeccapurple', zorder=1)

axs[1, 1].set_xlabel(r'$R \mathbin{/} \unit{\mega\parsec}$')
axs[1, 1].set_ylabel(r'$K \mathbin{/} \unit{\kilo\electronvolt\centi\meter\squared}$')

axs[1, 1].set_xscale('log')
axs[1, 1].set_yscale('log')

axs[1, 1].legend()

axs[0, 1].plot(R426[1:], P426[1:], label='A426', lw=0.9, alpha=0.75, c='firebrick', zorder=5)
axs[0, 1].plot(R2029[1:], P2029[1:], label='A2029', lw=0.9, alpha=0.75, c='goldenrod', zorder=4)
axs[0, 1].plot(R85[1:], P85[1:], label='A85', lw=0.9, alpha=0.75, c='olivedrab', zorder=3)
axs[0, 1].plot(R3158[1:], P3158[1:], label='A3158', lw=0.9, alpha=0.75, c='steelblue', zorder=2)
axs[0, 1].plot(R1644[1:], P1644[1:], label='A1644', lw=0.9, alpha=0.75, c='rebeccapurple', zorder=1)

axs[0, 1].set_xlabel(r'$R \mathbin{/} \unit{\mega\parsec}$')
axs[0, 1].set_ylabel(r'$P \mathbin{/} \unit{\erg\per\centi\meter\cubed}$')

axs[0, 1].set_xscale('log')
axs[0, 1].set_yscale('log')

axs[0, 1].legend()

plt.savefig('content/chan_profiles.pdf')

plt.close()

R85, Z85, Zerr85 = np.genfromtxt('catalog/XMM-Newton/X-COP-A85-abun.txt', unpack=True)
R1644, Z1644, Zerr1644 = np.genfromtxt('catalog/XMM-Newton/X-COP-A1644-abun.txt', unpack=True)
R2029, Z2029, Zerr2029 = np.genfromtxt('catalog/XMM-Newton/X-COP-A2029-abun.txt', unpack=True)
R3158, Z3158, Zerr3158 = np.genfromtxt('catalog/XMM-Newton/X-COP-A3158-abun.txt', unpack=True)

plt.plot(R2029, Z2029, label='A2029', lw=0.9, alpha=0.75, c='goldenrod', zorder=4)
plt.plot(R85, Z85, label='A85', lw=0.9, alpha=0.75, c='olivedrab', zorder=3)
plt.plot(R3158, Z3158, label='A3158', lw=0.9, alpha=0.75, c='steelblue', zorder=2)
plt.plot(R1644, Z1644, label='A1644', lw=0.9, alpha=0.75, c='rebeccapurple', zorder=1)

plt.xlabel(r'$R \mathbin{/} R_{500}$')
plt.ylabel(r'$Z$')

plt.xscale('log')

plt.legend()

plt.savefig('content/xmmn_abundance.pdf')

plt.close()

R85, ne85, nelo85, nehi85 = np.genfromtxt('catalog/XMM-Newton/X-COP-A85-dens.txt', unpack=True)
R1644, ne1644, nelo1644, nehi1644 = np.genfromtxt('catalog/XMM-Newton/X-COP-A1644-dens.txt', unpack=True)
R2029, ne2029, nelo2029, nehi2029 = np.genfromtxt('catalog/XMM-Newton/X-COP-A2029-dens.txt', unpack=True)
R3158, ne3158, nelo3158, nehi3158 = np.genfromtxt('catalog/XMM-Newton/X-COP-A3158-dens.txt', unpack=True)

plt.plot(R2029, ne2029, label='A2029', lw=0.9, alpha=0.75, c='goldenrod', zorder=4)
plt.plot(R85, ne85, label='A85', lw=0.9, alpha=0.75, c='olivedrab', zorder=3)
plt.plot(R3158, ne3158, label='A3158', lw=0.9, alpha=0.75, c='steelblue', zorder=2)
plt.plot(R1644, ne1644, label='A1644', lw=0.9, alpha=0.75, c='rebeccapurple', zorder=1)

plt.xlabel(r'$R \mathbin{/} \unit{\kilo\parsec}$')
plt.ylabel(r'$n \mathbin{/} \unit{\per\centi\meter\cubed}$')

plt.xscale('log')
plt.yscale('log')

plt.legend()

plt.savefig('content/xmmn_density.pdf')

plt.close()

R85, K85, Kerr85 = np.genfromtxt('catalog/XMM-Newton/X-COP-A85-entr.txt', unpack=True)
R1644, K1644, Kerr1644 = np.genfromtxt('catalog/XMM-Newton/X-COP-A1644-entr.txt', unpack=True)
R2029, K2029, Kerr2029 = np.genfromtxt('catalog/XMM-Newton/X-COP-A2029-entr.txt', unpack=True)
R3158, K3158, Kerr3158 = np.genfromtxt('catalog/XMM-Newton/X-COP-A3158-entr.txt', unpack=True)

plt.plot(R2029, K2029, label='A2029', lw=0.9, alpha=0.75, c='goldenrod', zorder=4)
plt.plot(R85, K85, label='A85', lw=0.9, alpha=0.75, c='olivedrab', zorder=3)
plt.plot(R3158, K3158, label='A3158', lw=0.9, alpha=0.75, c='steelblue', zorder=2)
plt.plot(R1644, K1644, label='A1644', lw=0.9, alpha=0.75, c='rebeccapurple', zorder=1)

plt.xlabel(r'$R \mathbin{/} R_{500}$')
plt.ylabel(r'$K \mathbin{/} K_{500}$')

plt.xscale('log')
plt.yscale('log')

plt.legend()

plt.savefig('content/xmmn_entropy.pdf')

plt.close()

R85, fgas85, fgaslo85, fgashi85 = np.genfromtxt('catalog/XMM-Newton/X-COP-A85-fgas.txt', unpack=True)
R1644, fgas1644, fgaslo1644, fgashi1644 = np.genfromtxt('catalog/XMM-Newton/X-COP-A1644-fgas.txt', unpack=True)
R2029, fgas2029, fgaslo2029, fgashi2029 = np.genfromtxt('catalog/XMM-Newton/X-COP-A2029-fgas.txt', unpack=True)
R3158, fgas3158, fgaslo3158, fgashi3158 = np.genfromtxt('catalog/XMM-Newton/X-COP-A3158-fgas.txt', unpack=True)

plt.plot(R2029, fgas2029, label='A2029', lw=0.9, alpha=0.75, c='goldenrod', zorder=4)
plt.plot(R85, fgas85, label='A85', lw=0.9, alpha=0.75, c='olivedrab', zorder=3)
plt.plot(R3158, fgas3158, label='A3158', lw=0.9, alpha=0.75, c='steelblue', zorder=2)
plt.plot(R1644, fgas1644, label='A1644', lw=0.9, alpha=0.75, c='rebeccapurple', zorder=1)

plt.xlabel(r'$R \mathbin{/} R_{500}$')
plt.ylabel(r'$f_g$', labelpad=-0.5)

plt.xscale('log')
plt.yscale('log')

plt.legend()

plt.savefig('content/xmmn_fgas.pdf')

plt.close()

R85, P85, Perr85 = np.genfromtxt('catalog/XMM-Newton/X-COP-A85-pres.txt', unpack=True)
R1644, P1644, Perr1644 = np.genfromtxt('catalog/XMM-Newton/X-COP-A1644-pres.txt', unpack=True)
R2029, P2029, Perr2029 = np.genfromtxt('catalog/XMM-Newton/X-COP-A2029-pres.txt', unpack=True)
R3158, P3158, Perr3158 = np.genfromtxt('catalog/XMM-Newton/X-COP-A3158-pres.txt', unpack=True)

plt.plot(R2029, P2029, label='A2029', lw=0.9, alpha=0.75, c='goldenrod', zorder=4)
plt.plot(R85, P85, label='A85', lw=0.9, alpha=0.75, c='olivedrab', zorder=3)
plt.plot(R3158, P3158, label='A3158', lw=0.9, alpha=0.75, c='steelblue', zorder=2)
plt.plot(R1644, P1644, label='A1644', lw=0.9, alpha=0.75, c='rebeccapurple', zorder=1)

plt.xlabel(r'$R \mathbin{/} R_{500}$')
plt.ylabel(r'$P \mathbin{/} P_{500}$')

plt.xscale('log')
plt.yscale('log')

plt.legend()

plt.savefig('content/xmmn_pressure.pdf')

plt.close()

R85, T85, Terr85 = np.genfromtxt('catalog/XMM-Newton/X-COP-A85-temp.txt', unpack=True)
R1644, T1644, Terr1644 = np.genfromtxt('catalog/XMM-Newton/X-COP-A1644-temp.txt', unpack=True)
R2029, T2029, Terr2029 = np.genfromtxt('catalog/XMM-Newton/X-COP-A2029-temp.txt', unpack=True)
R3158, T3158, Terr3158 = np.genfromtxt('catalog/XMM-Newton/X-COP-A3158-temp.txt', unpack=True)

plt.plot(R2029, T2029, label='A2029', lw=0.9, alpha=0.75, c='goldenrod', zorder=4)
plt.plot(R85, T85, label='A85', lw=0.9, alpha=0.75, c='olivedrab', zorder=3)
plt.plot(R3158, T3158, label='A3158', lw=0.9, alpha=0.75, c='steelblue', zorder=2)
plt.plot(R1644, T1644, label='A1644', lw=0.9, alpha=0.75, c='rebeccapurple', zorder=1)

plt.xlabel(r'$R \mathbin{/} R_{500}$')
plt.ylabel(r'$T \mathbin{/} T_{500}$')

plt.xscale('log')

plt.legend()

plt.savefig('content/xmmn_temperature.pdf')

plt.close()

fig, axs = plt.subplots(2, 3, figsize=[7.81, 3.90])

R85, Z85, Zerr85 = np.genfromtxt('catalog/XMM-Newton/X-COP-A85-abun.txt', unpack=True)
R1644, Z1644, Zerr1644 = np.genfromtxt('catalog/XMM-Newton/X-COP-A1644-abun.txt', unpack=True)
R2029, Z2029, Zerr2029 = np.genfromtxt('catalog/XMM-Newton/X-COP-A2029-abun.txt', unpack=True)
R3158, Z3158, Zerr3158 = np.genfromtxt('catalog/XMM-Newton/X-COP-A3158-abun.txt', unpack=True)

axs[1, 0].plot(R2029, Z2029, label='A2029', lw=0.9, alpha=0.75, c='goldenrod', zorder=4)
axs[1, 0].plot(R85, Z85, label='A85', lw=0.9, alpha=0.75, c='olivedrab', zorder=3)
axs[1, 0].plot(R3158, Z3158, label='A3158', lw=0.9, alpha=0.75, c='steelblue', zorder=2)
axs[1, 0].plot(R1644, Z1644, label='A1644', lw=0.9, alpha=0.75, c='rebeccapurple', zorder=1)

axs[1, 0].set_xlabel(r'$R \mathbin{/} R_{500}$')
axs[1, 0].set_ylabel(r'$Z$')

axs[1, 0].set_xscale('log')

axs[1, 0].legend()

R85, ne85, nelo85, nehi85 = np.genfromtxt('catalog/XMM-Newton/X-COP-A85-dens.txt', unpack=True)
R1644, ne1644, nelo1644, nehi1644 = np.genfromtxt('catalog/XMM-Newton/X-COP-A1644-dens.txt', unpack=True)
R2029, ne2029, nelo2029, nehi2029 = np.genfromtxt('catalog/XMM-Newton/X-COP-A2029-dens.txt', unpack=True)
R3158, ne3158, nelo3158, nehi3158 = np.genfromtxt('catalog/XMM-Newton/X-COP-A3158-dens.txt', unpack=True)

axs[0, 0].plot(R2029, ne2029, label='A2029', lw=0.9, alpha=0.75, c='goldenrod', zorder=4)
axs[0, 0].plot(R85, ne85, label='A85', lw=0.9, alpha=0.75, c='olivedrab', zorder=3)
axs[0, 0].plot(R3158, ne3158, label='A3158', lw=0.9, alpha=0.75, c='steelblue', zorder=2)
axs[0, 0].plot(R1644, ne1644, label='A1644', lw=0.9, alpha=0.75, c='rebeccapurple', zorder=1)

axs[0, 0].set_xlabel(r'$R \mathbin{/} \unit{\kilo\parsec}$')
axs[0, 0].set_ylabel(r'$n \mathbin{/} \unit{\per\centi\meter\cubed}$')

axs[0, 0].set_xscale('log')
axs[0, 0].set_yscale('log')

axs[0, 0].legend()

R85, K85, Kerr85 = np.genfromtxt('catalog/XMM-Newton/X-COP-A85-entr.txt', unpack=True)
R1644, K1644, Kerr1644 = np.genfromtxt('catalog/XMM-Newton/X-COP-A1644-entr.txt', unpack=True)
R2029, K2029, Kerr2029 = np.genfromtxt('catalog/XMM-Newton/X-COP-A2029-entr.txt', unpack=True)
R3158, K3158, Kerr3158 = np.genfromtxt('catalog/XMM-Newton/X-COP-A3158-entr.txt', unpack=True)

axs[1, 2].plot(R2029, K2029, label='A2029', lw=0.9, alpha=0.75, c='goldenrod', zorder=4)
axs[1, 2].plot(R85, K85, label='A85', lw=0.9, alpha=0.75, c='olivedrab', zorder=3)
axs[1, 2].plot(R3158, K3158, label='A3158', lw=0.9, alpha=0.75, c='steelblue', zorder=2)
axs[1, 2].plot(R1644, K1644, label='A1644', lw=0.9, alpha=0.75, c='rebeccapurple', zorder=1)

axs[1, 2].set_xlabel(r'$R \mathbin{/} R_{500}$')
axs[1, 2].set_ylabel(r'$K \mathbin{/} K_{500}$')

axs[1, 2].set_xscale('log')
axs[1, 2].set_yscale('log')

axs[1, 2].legend()

R85, fgas85, fgaslo85, fgashi85 = np.genfromtxt('catalog/XMM-Newton/X-COP-A85-fgas.txt', unpack=True)
R1644, fgas1644, fgaslo1644, fgashi1644 = np.genfromtxt('catalog/XMM-Newton/X-COP-A1644-fgas.txt', unpack=True)
R2029, fgas2029, fgaslo2029, fgashi2029 = np.genfromtxt('catalog/XMM-Newton/X-COP-A2029-fgas.txt', unpack=True)
R3158, fgas3158, fgaslo3158, fgashi3158 = np.genfromtxt('catalog/XMM-Newton/X-COP-A3158-fgas.txt', unpack=True)

axs[1, 1].plot(R2029, fgas2029, label='A2029', lw=0.9, alpha=0.75, c='goldenrod', zorder=4)
axs[1, 1].plot(R85, fgas85, label='A85', lw=0.9, alpha=0.75, c='olivedrab', zorder=3)
axs[1, 1].plot(R3158, fgas3158, label='A3158', lw=0.9, alpha=0.75, c='steelblue', zorder=2)
axs[1, 1].plot(R1644, fgas1644, label='A1644', lw=0.9, alpha=0.75, c='rebeccapurple', zorder=1)

axs[1, 1].set_xlabel(r'$R \mathbin{/} R_{500}$')
axs[1, 1].set_ylabel(r'$f_g$', labelpad=-0.5)

axs[1, 1].set_xscale('log')
axs[1, 1].set_yscale('log')

axs[1, 1].legend()

R85, P85, Perr85 = np.genfromtxt('catalog/XMM-Newton/X-COP-A85-pres.txt', unpack=True)
R1644, P1644, Perr1644 = np.genfromtxt('catalog/XMM-Newton/X-COP-A1644-pres.txt', unpack=True)
R2029, P2029, Perr2029 = np.genfromtxt('catalog/XMM-Newton/X-COP-A2029-pres.txt', unpack=True)
R3158, P3158, Perr3158 = np.genfromtxt('catalog/XMM-Newton/X-COP-A3158-pres.txt', unpack=True)

axs[0, 1].plot(R2029, P2029, label='A2029', lw=0.9, alpha=0.75, c='goldenrod', zorder=4)
axs[0, 1].plot(R85, P85, label='A85', lw=0.9, alpha=0.75, c='olivedrab', zorder=3)
axs[0, 1].plot(R3158, P3158, label='A3158', lw=0.9, alpha=0.75, c='steelblue', zorder=2)
axs[0, 1].plot(R1644, P1644, label='A1644', lw=0.9, alpha=0.75, c='rebeccapurple', zorder=1)

axs[0, 1].set_xlabel(r'$R \mathbin{/} R_{500}$')
axs[0, 1].set_ylabel(r'$P \mathbin{/} P_{500}$')

axs[0, 1].set_xscale('log')
axs[0, 1].set_yscale('log')

axs[0, 1].legend()

R85, T85, Terr85 = np.genfromtxt('catalog/XMM-Newton/X-COP-A85-temp.txt', unpack=True)
R1644, T1644, Terr1644 = np.genfromtxt('catalog/XMM-Newton/X-COP-A1644-temp.txt', unpack=True)
R2029, T2029, Terr2029 = np.genfromtxt('catalog/XMM-Newton/X-COP-A2029-temp.txt', unpack=True)
R3158, T3158, Terr3158 = np.genfromtxt('catalog/XMM-Newton/X-COP-A3158-temp.txt', unpack=True)

axs[0, 2].plot(R2029, T2029, label='A2029', lw=0.9, alpha=0.75, c='goldenrod', zorder=4)
axs[0, 2].plot(R85, T85, label='A85', lw=0.9, alpha=0.75, c='olivedrab', zorder=3)
axs[0, 2].plot(R3158, T3158, label='A3158', lw=0.9, alpha=0.75, c='steelblue', zorder=2)
axs[0, 2].plot(R1644, T1644, label='A1644', lw=0.9, alpha=0.75, c='rebeccapurple', zorder=1)

axs[0, 2].set_xlabel(r'$R \mathbin{/} R_{500}$')
axs[0, 2].set_ylabel(r'$T \mathbin{/} T_{500}$')

axs[0, 2].set_xscale('log')

axs[0, 2].legend()

plt.savefig('content/xmmn_profiles.pdf')

plt.close()
