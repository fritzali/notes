### preamble ###

import scipy as sci
import pandas as pd
import numpy as np
import matplotlib as mpl

import matplotlib.pyplot as plt

### spectra ###

E85, F85, Ferr85 = np.genfromtxt('catalog/Chandra/acisf15173N003_evt2-A85.txt', unpack=True)
E426, F426, Ferr426 = np.genfromtxt('catalog/Chandra/acisf04952N004_evt2-A426.txt', unpack=True)
E1644, F1644, Ferr1644 = np.genfromtxt('catalog/Chandra/acisf07922N003_evt2-A1644.txt', unpack=True)
E1656, F1656, Ferr1656 = np.genfromtxt('catalog/Chandra/acisf13996N003_evt2-A1656.txt', unpack=True)
E2029, F2029, Ferr2029 = np.genfromtxt('catalog/Chandra/acisf04977N006_evt2-A2029.txt', unpack=True)
E3158, F3158, Ferr3158 = np.genfromtxt('catalog/Chandra/acisf03712N004_evt2-A3158.txt', unpack=True)

plt.plot(E426, F426, label='A426 (Perseus)', lw=0.6, c='k', alpha=0.6)
plt.plot(E1656, F1656, label='A1656 (Coma)', lw=0.3, c='k', alpha=0.9)

plt.plot([], [], ' ', label=' ')
plt.plot([], [], ' ', label=' ')

plt.plot(E85, F85, label='A85', lw=0.6, alpha=0.6)
plt.plot(E1644, F1644, label='A1644', lw=0.6, alpha=0.6)
plt.plot(E2029, F2029, label='A2029', lw=0.6, alpha=0.6)
plt.plot(E3158, F3158, label='A3158  ', lw=0.6, alpha=0.6)

plt.xlabel(r'$E \mathbin{/} \unit{\kilo\electronvolt}$')
plt.ylabel(r'$N \mathbin{/} \unit{\per\second\per\kilo\electronvolt}$', labelpad=-0.5)

plt.yscale('log')

plt.legend(ncol=2)

plt.savefig('content/spectra.pdf')
