# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 20:55:01 2022

@author: rosha
"""

import matplotlib.pyplot as plt
import numpy as np

eps_moea_pf = [[-3.66561221e-01,  5.69331923e+03],
 [-3.41491226e-01, 5.66545976e+03],
 [-3.29747378e-01,  5.63338054e+03],
 [-3.23012522e-01,  5.54860060e+03],
 [-2.67929062e-01,  5.53114071e+03],
 [-2.10299837e-01,  5.49409908e+03]]

aos_pf = [[-3.68182286e-01,  6.62935232e+03],
 [-3.62545844e-01,  5.71306386e+03],
 [-3.31189749e-01,  5.49908566e+03],
 [-3.08943552e-01,  5.23653744e+03],
 [-2.72644059e-01,  5.20631472e+03],
 [-2.02430127e-01,  5.19954447e+03]]

eps_moea_pf_science = [x[0] for x in eps_moea_pf]
eps_moea_pf_cost = [x[1] for x in eps_moea_pf]

aos_pf_science = [x[0] for x in aos_pf]
aos_pf_cost = [x[1] for x in aos_pf]

plt.figure()
plt.scatter(eps_moea_pf_science, eps_moea_pf_cost, marker='*', color='#000000', label='Eps. MOEA')
plt.scatter(aos_pf_science, aos_pf_cost, marker='*', color='#E69F00', label='AOS')
plt.xlabel('Science')
plt.ylabel('Cost')
plt.legend(loc='best')