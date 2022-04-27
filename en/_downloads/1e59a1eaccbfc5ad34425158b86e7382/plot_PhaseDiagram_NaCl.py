# -*- coding: utf-8 -*-
"""
1. Phase diagram
==================================================
.. include:: /include.rst_
Calculate phase boundary of pure NaCl, and compare result with :cite:`Driesner2007Part1`.
"""
import os
import numpy as np
import time
import linecache
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import patches
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator
from tabulate import tabulate
import copy
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
mpl.rcParams['font.family'] = 'Arial'  # default font family
mpl.rcParams['mathtext.fontset'] = 'cm'  # font for math
dpi=100 # > 100 means use jpg format in gallery
fmt_figs = ['pdf']  # ['svg','pdf']
figpath = '.'
result_path='../../../gallery_NaCl/pT'
def savefig(figname):
    for fmt_fig in fmt_figs:
        figname_full = '%s/%s.%s' % (figpath, figname, fmt_fig)
        plt.savefig(figname_full, bbox_inches='tight',dpi=dpi)
        print('figure saved: ', figname_full)
compare = lambda a,b : float(str('%.6e'%(a)))-float(str('%.6e'%(b)))
# Import package of xThermo
from xThermo import H2O
from xThermo import NaCl
from xThermo import H2ONaCl
sw_84 = H2ONaCl.cH2ONaCl("IAPS84")
sw_95 = H2ONaCl.cH2ONaCl("IAPWS95")
salt_84 = NaCl.cNaCl("IAPS84")
salt_95 = NaCl.cNaCl("IAPWS95")

def plot_PhaseDiagram(salt,sw, scale='linear'):
    # Fig. 4 of Driesner(2007a)
    T = np.linspace(300, 1100, 100) + 273.15
    P_boiling = np.array(salt.Boiling_p(T))
    P_sublimation = np.array(salt.Sublimation_p(T))
    P_melting = np.linspace(NaCl.P_Triple, 5000E5, 100)
    T_melting = np.array(salt.Melting_T(P_melting))
    fig=plt.figure(figsize=(8,8))
    ax=plt.gca()
    # boiling curve
    ind_boil=(T>=NaCl.T_Triple)
    l,=ax.semilogy(T[ind_boil]-273.15, P_boiling[ind_boil]/1E5, color='orange',label='Boiling curve')
    # ax.semilogy(np.array(salt.Boiling_T(P_boiling[ind_boil]))-273.15, P_boiling[ind_boil]/1E5, color='red',ls='dashed',lw=0.5)
    # sublimation curve
    ind_sub=(T<=NaCl.T_Triple)
    l,=ax.semilogy(T[ind_sub]-273.15, P_sublimation[ind_sub]/1E5, color='green',label='Sublimation curve')
    # melting curve
    l,=ax.semilogy(T_melting-273.15, P_melting/1E5, color='m',label='Melting curve')
    # triple point
    ax.plot(NaCl.T_Triple-273.15, NaCl.P_Triple/1E5, 'o', mfc='red',mec='w', label='Triple point: (%.2f $^{\circ}$C, %.2f Pa)'%(NaCl.T_Triple-273.15, NaCl.P_Triple))
    # valid pressure of H2O-NaCl
    ax.axhline(sw.pmin()/1E5,color='r',label='Valid $P_{min}$ = %.1f bar\n$H2O-NaCl$ lib'%(sw.pmin()/1E5),ls='dashed')
    # phase region
    ax.fill_betweenx(np.append(P_sublimation[ind_sub], P_melting)/1E5, np.append(T[ind_sub], T_melting)-273.15, T.min()-273.15, fc='lightblue',ec='None',label='Solid phase',zorder=0)
    ax.fill(np.append(T[ind_boil][::-1], np.append(T_melting,T.max()))-273.15,np.append(P_boiling[ind_boil][::-1], np.append(P_melting,P_melting.max()))/1E5, T.min()-273.15, fc='lightgreen',ec='None',label='Liquid phase',zorder=0)
    ax.fill_betweenx(np.append(P_sublimation[ind_sub], P_boiling[ind_boil])/1E5, np.append(T[ind_sub], T[ind_boil])-273.15, T.max()-273.15, fc='gray',ec='None',label='Vapor phase',zorder=0)
    # invalid region
    ax.fill_between(T-273.15, T*0 + sw.pmin()/1E5,0,fc='w',ec='None',alpha=0.5)
    ax.legend(ncol=2,loc='lower right')
    ax.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, numticks=20))
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.set_xlim(T.min()-273.15, T.max()-273.15)
    ax.set_ylim(2E-13, 5000)
    ax.grid(lw=0.2, color='gray')
    ax.set_xlabel('Temperature ($^{\circ}$C)')
    ax.set_ylabel('P (bar)')
    savefig('PhaseDiagram_NaCl')
# %%
# Phase diagram
# --------------------------------
plot_PhaseDiagram(salt_84,sw_84,'linear')

# %%
# Benchmark comparison
# --------------------------
# Compare result of |xThermo| with :cite:`Driesner2007Part1`, and also compare result based on different EOS of water.
def benchmark_NaCl(salt,mmc1='../../H2ONaCl/Driesner2007a/1-s2.0-S0016703707002943-mmc1.txt'):
    # compare
    if(not os.path.exists(mmc1)):
        print('Please set correct mmc1 file path: %s'%(mmc1))
        exit()
    T0_melting, P0_melting = np.loadtxt(mmc1,skiprows=14,max_rows=64-14,unpack=True)
    T0_sub, P0_sub = np.loadtxt(mmc1,skiprows=78,max_rows=400-78,unpack=True)
    T0_boil, P0_boil = np.loadtxt(mmc1,skiprows=414,max_rows=496-414,unpack=True)
    # 1. calculate halite liquidus
    P_boil = np.array(salt.Boiling_p(T0_boil+273.15))/1E5
    P_sub = np.array(salt.Sublimation_p(T0_sub+273.15))/1E5
    P_melting = np.array(salt.Melting_p(T0_melting+273.15))/1E5
    # compare result dict
    Data0 = {'Melting P': P0_melting,'Boiling P': P0_boil,'Sublimation P': P0_sub}
    Data_ = {'Melting P': P_melting,'Boiling P': P_boil,'Sublimation P': P_sub}
    Err,RErr={},{}
    for key in Data0.keys(): Err[key],RErr[key] = Data0[key]-Data_[key], np.abs(Data0[key]-Data_[key])/(Data0[key])*100.0
    fig=plt.figure(figsize=(8,8))
    ax=plt.gca()
    # print to file
    for name,T0,P0,P in zip(['Melting','Boiling','Sublimation'],[T0_melting,T0_boil,T0_sub],[P0_melting,P0_boil,P0_sub],[P_melting,P_boil,P_sub]):
        fpout = open('%s/mmc1_%s_%s.csv'%(result_path,name,salt.name_backend()),'w')
        fpout.write('T[C],P(Driesner)[bar],P(xThermo),P(err)\n')
        for i in range(0,len(T0)):
            fpout.write('%f,%.6e,%.6e,%.6e\n'%(T0[i],P0[i],P[i],compare(P0[i],P[i])))
        fpout.close()
        ax.plot(T0,P,lw=4,label=name)
        ax.plot(T0,P0,'.',mfc='r',mec='w',markeredgewidth=0.5)
    ax.set_yscale('log')
    ax.legend()
    ax.set_ylim(1E-5,5000)
    ax.set_xlim(600,1000)
    ax.set_xlabel('Temperature ($^{\circ}$C)')
    ax.set_ylabel('P (bar)')
    savefig('Diff_PhaseBoundary')
    # statistics of the difference
    table=[]
    for key in list(Err.keys()):
        RErr[key] = RErr[key][~(np.isnan(RErr[key]) | np.isinf(RErr[key]))]
        table.append([key,Err[key].min(),Err[key].max(),RErr[key].min(),RErr[key].max()])
    print(tabulate(table, headers=['Pressure (bar)', 'Err. Min', 'Err. Max','RE. Min(%)','RE. Max(%)'],numalign="right",floatfmt=".6f"))

# %%
# Based on IAPS84 EOS
# ^^^^^^^^^^^^^^^^^^^^^^^
benchmark_NaCl(salt_84)

# %%
# Based on IAPWS95 EOS
# ^^^^^^^^^^^^^^^^^^^^^^^
benchmark_NaCl(salt_95)

# %%
# Result table
# -----------------------------
# .. seealso::
#
#     |mmc1| in :cite:`Driesner2007Part1` and Fig. 4 in :cite:`Driesner2007Part1`.
#
# Result data calculated by |xThermo| based on both water EOS of IAPS84 and IAPWS95.
#
# .. tab:: Melting curve
#
#     .. tab:: IAPS84
#
#         .. csv-table:: Comparison between result of :cite:`Driesner2007Part1` (|mmc1|) and result calculated by |xThermo|.
#             :file: mmc1_Melting_IAPS84.csv
#             :header-rows: 1
#
#     .. tab:: IAPWS95
#
#         .. csv-table:: Comparison between result of :cite:`Driesner2007Part1` (|mmc1|) and result calculated by |xThermo|.
#             :file: mmc1_Melting_IAPWS95.csv
#             :header-rows: 1
#
# .. tab:: Boiling curve
#
#     .. tab:: IAPS84
#
#         .. csv-table:: Comparison between result of :cite:`Driesner2007Part1` (|mmc1|) and result calculated by |xThermo|.
#             :file: mmc1_Boiling_IAPS84.csv
#             :header-rows: 1
#
#     .. tab:: IAPWS95
#
#         .. csv-table:: Comparison between result of :cite:`Driesner2007Part1` (|mmc1|) and result calculated by |xThermo|.
#             :file: mmc1_Boiling_IAPWS95.csv
#             :header-rows: 1
#
# .. tab:: Sublimation curve
#
#     .. tab:: IAPS84
#
#         .. csv-table:: Comparison between result of :cite:`Driesner2007Part1` (|mmc1|) and result calculated by |xThermo|.
#             :file: mmc1_Sublimation_IAPS84.csv
#             :header-rows: 1
#
#     .. tab:: IAPWS95
#
#         .. csv-table:: Comparison between result of :cite:`Driesner2007Part1` (|mmc1|) and result calculated by |xThermo|.
#             :file: mmc1_Sublimation_IAPWS95.csv
#             :header-rows: 1