print("\033[J\033[H",end='') # Clear screen

import matplotlib as mpl
import os
import numpy as np
import file_utils

mpl.pyplot.rc('text', usetex=True)

fpath = os.getcwd().removesuffix('/tests') + "/Data1/frame0.dat"
[id,rx,ry,px,py,fx,fy,jumpx,jumpy,Pe] = file_utils.readData(fpath)

y = px

y = [i-min(y) for i in y]
leny = len(y)
bin_w = 0.1
nbin = int(max(y)/bin_w) + 1
print(nbin)
bins = [0]*nbin

for i in range(leny):
    bins[int(y[i]/bin_w)] += 1
bins = [i/leny for i in bins]

x = np.arange((-1*nbin*bin_w/2.0),(nbin*bin_w/2.0),bin_w)
mpl.pyplot.plot(x,bins)