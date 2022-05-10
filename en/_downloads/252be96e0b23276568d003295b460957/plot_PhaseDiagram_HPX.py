# -*- coding: utf-8 -*-
"""
0. Phase diagram
==================================================
.. include:: /include.rst_
Calculate and plot phase diagram.
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
from matplotlib.patches import Patch
import copy
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
# 3d plot
import helpfunc
mpl.rcParams['font.family'] = 'Arial'  # default font family
mpl.rcParams['mathtext.fontset'] = 'cm'  # font for math
dpi=100
fmt_figs = ['pdf']  # ['svg','pdf']
figpath = '.'
result_path='../../../gallery_H2ONaCl/pT'
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

def plot_Phase(ax, XX,YY,Phase,sw, cmap="Dark2"):
    cmap = plt.get_cmap(cmap)
    # customize cmap
    phase_unique = np.sort(np.unique(Phase))
    phase_name = ['']*len(phase_unique)
    for i,phase0 in enumerate(phase_unique):
        Phase[Phase==phase0]=i+phase_unique.max()+10
        phase_name[i]=sw.phase_name(int(phase0))
    colors=list(copy.deepcopy(cmap.colors))
    colors[0:8]=['red','lightblue','lightgreen','lightgray','violet','yellow','lightcyan','k']
    cmap.colors=tuple(colors)
    CS=ax.contourf(XX,YY, Phase,cmap=cmap,vmin=Phase.min()-0.5, vmax=Phase.max()+0.5, levels=np.linspace(Phase.min()-0.5,Phase.max()+0.5,len(phase_name)+1))
    ax_cb = ax.inset_axes([0,1.01,1,0.05])
    cb=plt.colorbar(CS, cax=ax_cb, orientation='horizontal',ticklocation='top',ticks=np.arange(Phase.min(),Phase.max()+1))
    cb.ax.set_xticklabels(phase_name)


def IsobaricSection(sw, p0):
    # # H = np.linspace(0.1, 3, 100)*1E6
    # X = np.linspace(0.00000001,0.9999,100)
    # T=np.linspace(H2ONaCl.T_MIN+0.1, H2ONaCl.T_MAX, 100)
    # TT,XX = np.meshgrid(T,X)
    # PP = TT*0 + p0
    # props = sw.UpdateState_TPX(TT.reshape(-1,), PP.reshape(-1,), XX.reshape(-1,))
    # HH = np.array(props.H).reshape(TT.shape)
    # TT_=np.zeros_like(TT)
    # for i in range(0,TT.shape[0]):
    #     for j in range(0,TT.shape[1]):
    #         # print('py T=%f, H=%f\n'%(TT[i][j], HH[i][j]))
    #         TT_[i][j] = sw.T_HPX(HH[i][j], PP[i][j], XX[i][j])

    H = np.linspace(0.1, 5.5, 200)*1E6
    # X = np.linspace(0.00000001,0.9999,100)
    X = np.linspace(1E-8,1,200)
    # T=np.linspace(H2ONaCl.T_MIN+0.1, H2ONaCl.T_MAX, 100)
    HH,XX = np.meshgrid(H,X)
    PP = HH*0 + p0
    TT_=np.zeros_like(HH)
    phase = np.zeros_like(HH)
    for i in range(0,HH.shape[0]):
        for j in range(0,HH.shape[1]):
            # print('py T=%f, H=%f\n'%(TT[i][j], HH[i][j]))
            # TT_[i][j] = sw.T_HPX(HH[i][j], PP[i][j], XX[i][j])
            # phase[i][j] = sw.findPhaseRegion_HPX(HH[i][j], PP[i][j], XX[i][j])
            props = sw.UpdateState_HPX(HH[i][j], PP[i][j], XX[i][j])
            phase[i][j] = props.phase
    # calculate phase
    # props = sw.UpdateState_TPX(TT_.reshape(-1,), PP.reshape(-1,), XX.reshape(-1,))
    # phase = np.array(props.phase).reshape(HH.shape)
    # plot
    fig=plt.figure(figsize=(7,6))
    ax=plt.gca()
    ax.contourf(XX,HH,phase,levels=50)
    plot_Phase(ax, XX*100,HH,phase, sw)
    savefig('phase_HPX_p%.0fbar'%(p0/1E5))


# %%
# Isobaric section
# --------------------------
IsobaricSection(sw_84, 250E5)