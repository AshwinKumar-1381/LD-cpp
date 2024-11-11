"""
Author          : Ashwin Kumar M
# Date created  : 09.11.24
# Last modified : 09.11.24
"""

print("\033[J\033[H",end='') # Clear screen

import matplotlib as mpl
import os
import math as m
import file_utils 

mpl.pyplot.rc('text', usetex=True)

def avg(x): return(sum(x)/len(x))

option = 1

if(option == 1):
    nr = 1
    runID = 2
    fpath = os.getcwd().removesuffix("/tests") + "/Data{nr}/thermo{rid}.dat" \
        .format(nr = nr,rid = runID)
    [step, pe, ke, etot, temp] = file_utils.readData(fpath)

    print("# data points = ", len(etot))
    H = etot
    H2 = [i*i for i in H]
    print("<H> = ", avg(H))
    print("<H>^2 = ", avg(H)**2)
    print("<H2> = ",avg(H2))
    print("rms delH2 = ", m.sqrt(avg(H2) - avg(H)**2))
    
    T = temp
    T2 = [i*i for i in T]
    print("<H> = ", avg(T))
    print("<H>^2 = ", avg(T)**2)
    print("<H2> = ",avg(T2))
    print("rms delH2 = ", m.sqrt(avg(T2) - avg(T)**2))

    fig, ax = mpl.pyplot.subplots()
    ax.semilogy(step, pe, lw = 0.7, label = "pe")
    ax.semilogy(step, ke, lw = 0.7, label = "ke")
    ax.semilogy(step, etot, lw = 0.7, label = "etot")
    ax.set_xlabel("Step")
    ax.set_ylabel("Per-atom energy")
    ax.legend()

else:
    
    [dt, delH2] = file_utils.readData("del_ETOT.dat")
    
    fig, ax = mpl.pyplot.subplots()
    ax.loglog(dt, delH2, "o-", lw = 0.6, ms = 2)
    ax.set_xlabel(r"$dt$", fontsize = 12)
    ax.set_ylabel(r"$\sqrt{<\delta \mathcal{H}^2>}$", fontsize = 12)
    