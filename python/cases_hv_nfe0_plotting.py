# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 21:57:29 2022

@author: rosha
"""

import matplotlib.pyplot as plt
import statistics
import numpy as np

eps_moea_hv = [0.65030977, 0.64177848, 0.70823343, 0.6462806,  0.67212395, 0.68623182,
 0.76396479, 0.6634377,  0.58880373, 0.67985629, 0.63096641, 0.72828856,
 0.71668831, 0.66660093, 0.71234899, 0.69977954, 0.71982867, 0.66566005,
 0.68210988, 0.68180619, 0.65146656, 0.72943225, 0.72358905, 0.76781782,
 0.66127728, 0.73619259, 0.77905526, 0.66934961, 0.65731165, 0.61904015]

aos_hv = [0.74165761, 0.77463109, 0.66654981, 0.8177392,  0.77873774, 0.6735262,
 0.69695101, 0.75146659, 0.62579524, 0.78265373, 0.7008309,  0.63573625,
 0.67244108, 0.72076158, 0.69549114, 0.69730609, 0.70497024, 0.82382158,
 0.75439242, 0.74218958, 0.70039003, 0.65372997, 0.73787065, 0.71342337,
 0.62573084, 0.83868439, 0.74036317, 0.74600215, 0.67144605, 0.70349843]

print('Eps MOEA:')
print('Median HV:')
print(statistics.median(eps_moea_hv))
print('First Quartile')
print(np.percentile(eps_moea_hv, 25))
print('Third Quartile')
print(np.percentile(eps_moea_hv, 75))

print('AOS:')
print('Median HV:')
print(statistics.median(aos_hv))
print('First Quartile')
print(np.percentile(aos_hv, 25))
print('Third Quartile')
print(np.percentile(aos_hv, 75))

plt.figure
plt.plot(eps_moea_hv, marker='*', color='#000000', label='Eps. MOEA')
plt.plot(aos_hv, marker='*', color='#E69F00', label='AOS')
plt.xlabel('Run Number')
plt.ylabel('HV at NFE=0')
plt.legend(loc='best')