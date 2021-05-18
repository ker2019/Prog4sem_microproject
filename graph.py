#!/usr/bin/env python3

import matplotlib
matplotlib.use('GTK3Cairo')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys


data = pd.read_csv(sys.argv[1], sep='\t')
k = data['k'].drop_duplicates()[0]

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

x = float(sys.argv[2])
df = data[data['x'] == x]
reu = df['reu']
imu = df['imu']
y = df['y']
#ax.plot(y, reu*np.cos(k*x) - imu*np.sin(k*x))
j = 1j
ax.plot(y, abs((reu + imu*j)*np.exp(k*x*j)))

#fig.savefig('out/difraction.png', format='png')
plt.show()
