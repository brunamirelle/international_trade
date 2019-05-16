#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 00:35:07 2019

@author: brunamirelle
"""

from EconomiaArtuc import EconomiaArtuc
import numpy as np
import matplotlib.pyplot as plt

economia_exercicio = EconomiaArtuc(0.5,0.97, 1, 0.31, 1, 1, 2, 1)
    
#Resolve modelo para uma liberalização efetiva que move preço para 0.7
#Sem delay
W_X_nd,W_Y_nd,V_X_nd,V_Y_nd,L_X_nd,L_Y_nd = economia_exercicio.solve_model(0.7)

#DElay de 1 período
W_X_1d,W_Y_1d,V_X_1d,V_Y_1d,L_X_1d,L_Y_1d = economia_exercicio.solve_model(0.7, effective_at = 1)

#Delay de 10 períodos
W_X_d,W_Y_d,V_X_d,V_Y_d,L_X_d,L_Y_d = economia_exercicio.solve_model(0.7, effective_at = 10)


tempo = np.arange(-1,32,1) 

plt.plot(tempo, W_X_nd, '-', color = 'blue', label = 'Without delay')
plt.plot(tempo, W_X_1d, '.-', color = 'black', label = 'Delay 1 period')
plt.plot(tempo, W_X_d, '--', color = 'red', label = 'Delay 10 periods')
plt.legend()
plt.title('Real Wages - Sector X')
plt.grid()
plt.savefig('wages_x.pdf')
plt.show()

plt.plot(tempo, W_Y_nd, '-', color = 'blue', label = 'Without delay')
plt.plot(tempo, W_Y_1d, '.-', color = 'black', label = 'Delay 1 period')
plt.plot(tempo, W_Y_d, '--', color = 'red', label = 'Delay 10 periods')
plt.legend()
plt.title('Real Wages - Sector Y')
plt.grid()
plt.savefig('wages_y.pdf')
plt.show()


plt.plot(tempo, L_X_nd, '-', color = 'blue', label = 'Without delay')
plt.plot(tempo, L_X_1d, '.-', color = 'black', label = 'Delay 1 period')
plt.plot(tempo, L_X_d, '--', color = 'red', label = 'Delay 10 periods')
plt.legend()
plt.title('Number of workers - Sector X')
plt.grid()
plt.savefig('workers_x.pdf')
plt.show()

plt.plot(tempo, L_Y_nd, '-', color = 'blue', label = 'Without delay')
plt.plot(tempo, L_Y_1d, '.-', color = 'black', label = 'Delay 1 period')
plt.plot(tempo, L_Y_d, '--', color = 'red', label = 'Delay 10 periods')
plt.legend()
plt.title('number of workers - Sector Y')
plt.grid()
plt.savefig('workers_y.pdf')
plt.show()


plt.plot(tempo, V_X_nd, '-', color = 'blue', label = 'Without delay')
plt.plot(tempo, V_X_1d, '.-', color = 'black', label = 'Delay 1 period')
plt.plot(tempo, V_X_d, '--', color = 'red', label = 'Delay 10 periods')
plt.legend()
plt.title('Value functions - Sector X')
plt.grid()
plt.savefig('value_x.pdf')
plt.show()

plt.plot(tempo, V_Y_nd, '-', color = 'blue', label = 'Without delay')
plt.plot(tempo, V_Y_1d, '.-', color = 'black', label = 'Delay 1 period')
plt.plot(tempo, V_Y_d, '--', color = 'red', label = 'Delay 10 periods')
plt.legend()
plt.title('Value functions - Sector Y')
plt.grid()
plt.savefig('value_y.pdf')
plt.show()