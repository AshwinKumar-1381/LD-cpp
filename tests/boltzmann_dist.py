print("\033[J\033[H",end='') # Clear screen

import matplotlib as mpl
import os
from ast import literal_eval as liteval
import math as m
import numpy as np

mpl.pyplot.rc('text', usetex=True)

def readData(fpath):
    fobj = open(fpath, mode = 'r', encoding = "utf-8")
    nCols = len(fobj.readline().removesuffix("\n").split(sep = " "))
    Cols = [[] for i in range(nCols)]
    
    for line in fobj:
        line = line.removesuffix("\n").split(sep = " ")
        for i in range(nCols): 
            if(i==0): Cols[i].append(line[i])
            else: Cols[i].append(liteval(line[i]))
            
    fobj.close()
    return(Cols)

fpath = os.getcwd().removesuffix('/code/misc') + "/LD/LD-cpp/Data1/frame0.dat"
[id,rx,ry,px,py,fx,fy,jumpx,jumpy,Pe] = readData(fpath)

y = py
y = [i-min(y) for i in y]
leny = len(y)
bin_w = 0.1
nbin = int(m.ceil(max(y))/bin_w)
print(nbin)
bins = [0]*nbin

for i in range(leny):
    bins[int(y[i]/bin_w)] += 1

mpl.pyplot.plot(np.arange(-4,4,bin_w),bins)

bins = [i/leny for i in bins]