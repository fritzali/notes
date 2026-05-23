### preamble ###

import scipy as sci
import pandas as pd
import numpy as np
import matplotlib as mpl

import matplotlib.pyplot as plt

### continua ###

comaE, comaF, comaFerr = np.genfromtxt('catalog/Chandra/acisf13996N003_evt2.txt', unpack=True)
persE, persF, persFerr = np.genfromtxt('catalog/Chandra/acisf04952N004_evt2.txt', unpack=True)
testE, testF, testFerr = np.genfromtxt('catalog/Chandra/acisf15173N003_evt2.txt', unpack=True)

plt.plot(persE, persF, lw=0.85, label='Perseus')
plt.plot(comaE, comaF, lw=0.85, label='Coma')
plt.plot(testE, testF, lw=0.85, label='A85')

plt.xlabel(r'$E \mathbin{/} \unit{\kilo\electronvolt}$')
plt.ylabel(r'$N \mathbin{/} \unit{\per\second\per\kilo\electronvolt}$', labelpad=-10)

plt.yscale('log')

#plt.gca().yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%g'))

plt.legend()

plt.savefig('content/continua.pdf')
