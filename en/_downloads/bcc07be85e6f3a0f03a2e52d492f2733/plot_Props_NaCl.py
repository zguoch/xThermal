# -*- coding: utf-8 -*-
"""
2. Density and specific enthalpy
==================================================
.. include:: /include.rst_
Calculate density and specific enthalpy ofof pure NaCl.
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

def plot_Props(salt, scale='linear',enthalpy_FV='../../../gallery/H2ONaCl/Driesner2007b/Halite_enthalpy_Vehling20220112.txt'):
    # Fig. 4 of Driesner(2007a)
    T = np.linspace(0, 1100, 300) + 273.15
    P = np.linspace(1,5000,300)*1E5
    TT,PP = np.meshgrid(T,P)
    RHO,H = np.zeros_like(TT),np.zeros_like(TT)

    # calculate
    for i in range(0,TT.shape[0]):
        for j in range(0,TT.shape[1]):
            props = salt.UpdateState_TPX(TT[i][j], PP[i][j])
            RHO[i][j]=props.Rho
            H[i][j]=props.H
    fig,axes=plt.subplots(1,2,figsize=(16,7))
    for ax,prop,label,cmap in zip(axes,[RHO,H/1E6],['Density (kg/m$^{\mathregular{3}}$)','Specific enthalpy (MJ/kg)'],['rainbow','jet']):
        CS=ax.contourf(TT-273.15,PP/1E5,prop,levels=50,cmap=cmap)
        ax.set_yscale(scale)
        ax_cb=ax.inset_axes([0,1.02,1,0.05])
        plt.colorbar(CS,cax=ax_cb,label=label,orientation='horizontal',extend='both')
        ax_cb.xaxis.set_label_position('top')
        ax_cb.xaxis.set_ticks_position('top')
        # ax.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, numticks=20))
        ax.xaxis.set_major_locator(MultipleLocator(100))
        ax.set_xlim(T.min()-273.15, T.max()-273.15)
        ax.set_ylim(1, 5000)
        ax.grid(lw=0.2, color='gray')
        ax.set_xlabel('Temperature ($^{\circ}$C)')
        ax.set_ylabel('P (bar)')
    savefig('Props_NaCl_%s'%(salt.name_backend()))

    # 2. compare result with F.V. result
    # load T,P,H result from Falko Vehling table (personal communication)
    T0,P0,H0 = np.loadtxt(enthalpy_FV,skiprows=1,unpack=True)
    H_=np.zeros_like(H0)
    fpout = open('%s/Enthalpy_pT_%s.csv'%(result_path,salt.name_backend()),'w')
    fpout.write('T[C],P[bar],H(F.V.)[J/kg],H(xThermo),H(err)\n')
    for i in range(0,len(T0)):
        props = salt.UpdateState_TPX(T0[i]+273.15,P0[i]*1E5)
        H_[i]=props.H
        fpout.write('%.0f,%.0f,%.6e,%.6e,%.6e\n'%(T0[i],P0[i],H0[i],H_[i],compare(H0[i],H_[i])))
    fpout.close()
    # plot difference
    fig=plt.figure(figsize=(8,8))
    ax=plt.gca()
    # colormap
    norm = mpl.colors.CenteredNorm(vcenter=0)
    cmap=copy.deepcopy(plt.get_cmap('bwr'))
    cmap._segmentdata['red'][0,1:]=0
    cmap._segmentdata['green'][0,1:]=1
    cmap._segmentdata['blue'][0,1:]=0
    CS=ax.scatter(T0,P0, c=H0-H_,cmap=cmap,norm=norm)
    ax_cb=ax.inset_axes([0,1.02,1,0.05])
    plt.colorbar(CS,cax=ax_cb,label='Specific enthalpy (J/kg): F.V. result - xThermo',orientation='horizontal',extend='both')
    ax_cb.xaxis.set_label_position('top')
    ax_cb.xaxis.set_ticks_position('top')
    ax.patch.set_facecolor('k')
    ax.set_xlabel('Temperature ($^{\circ}$C)')
    ax.set_ylabel('P (bar)')
    savefig('Diff_NaCl_H')

# %%
# Based on IAPS84 EOS
# ^^^^^^^^^^^^^^^^^^^^^^^
plot_Props(salt_84)

# %%
# Compare specific result with F.V.'s result
# ----------------------------------------------
# Result data calculated by |xThermo| based on both water EOS of IAPS84 and IAPWS95.
#
# .. tab:: Specific enthalpy
#
#     .. csv-table:: Comparison between result of Falko Vehling(personal communication) and result calculated by |xThermo|.
#         :file: Enthalpy_pT_IAPS84.csv
#         :header-rows: 1