"""
Author          : Ashwin Kumar M
Date created    : 20.11.24 
Last modified   : 22.11.24
"""

#print("\033[J\033[H",end='') # Clear screen

import file_utils
import math as m
import matplotlib as mpl
import numpy as np
"""
# comparing generated data with exp distribution
[id, tau] = file_utils.readData("exp_data_linked_list.dat")

nbins = 50
binw = max(tau)/nbins
bin = np.zeros(nbins+1)
Tau = np.arange(0,nbins+1,1)

for i in tau: bin[int(i/binw)] += 1
bin = [i/(len(tau)*binw) for i in bin]

mpl.pyplot.plot(Tau, bin, lw = 1)

rate = 1*5e-4
exp_dist = [rate*m.exp(-rate*binw*i) for i in Tau]

mpl.pyplot.plot(Tau, exp_dist, "--", lw = 1)
"""
rate = 1*5e-4
x = np.arange(0, 1e4, 10)
exp_dist = [rate*m.exp(-rate*i) for i in x]

fig, ax = mpl.pyplot.subplots()
ax.plot(x, exp_dist, "-", lw = 1)