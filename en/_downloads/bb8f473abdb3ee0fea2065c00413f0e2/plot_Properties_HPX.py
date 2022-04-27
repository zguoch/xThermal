# -*- coding: utf-8 -*-
"""
1. Properties
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
sw_95_CoolProp = H2ONaCl.cH2ONaCl("IAPWS95_CoolProp")

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

def plot_prop(ax, XX,YY,prop, prop_name='prop', prop_unit='unit', cmap='rainbow',levels=50):
    CS=ax.contourf(XX,YY,prop, levels=levels, cmap=cmap)
    ax_cb = ax.inset_axes([0,1.01,1,0.05])
    cb=plt.colorbar(CS, cax=ax_cb, orientation='horizontal',ticklocation='top',label='%s (%s)'%(prop_name, prop_unit))

def plot_props_P0X0(sw):
    P0,X0 = 100E5,0.1
    T=np.linspace(sw.Tmin(),sw.Tmax(),100)
    X = T*0 + X0
    P = T*0 + P0
    state = sw.UpdateState_TPX(T,P,X)
    # plot
    fig=plt.figure()
    ax=plt.gca()
    ax.plot(T, state.H,'o')
    savefig('Props_P%.0fbar_X%.1f'%(P0/1E5,X0))

def plot_props_H0(H0, sw):
    X = np.linspace(1E-6,1,400)
    P = np.linspace(1, 1000, 400)*1E5
    PP,XX = np.meshgrid(P,X)
    HH = PP*0 + H0
    phase = np.zeros_like(HH)
    Rho = np.zeros_like(HH)
    T = np.zeros_like(HH)
    Cp = np.zeros_like(HH)
    Mu_l = np.zeros_like(HH)
    Mu_v = np.zeros_like(HH)
    for i in range(0,HH.shape[0]):
        for j in range(0,HH.shape[1]):
            props = sw.UpdateState_HPX(HH[i][j], PP[i][j], XX[i][j])
            phase[i][j] = props.phase
            Rho[i][j] = props.Rho
            Cp[i][j] = props.Cp
            T[i][j] = props.T
            Mu_l[i][j] = props.Mu_l
            Mu_v[i][j] = props.Mu_v
    T[T<273.15]=np.nan
    # plot
    fig, axes=plt.subplots(2,4,figsize=(28,12),sharey=True, sharex=True,gridspec_kw={'wspace':0.1,'hspace':0.2})

    xx,yy = XX*100, PP/1E5
    # 1. phase
    ax = axes[0][0]
    plot_Phase(ax, xx,yy, phase, sw)
    # 2. Rho
    ax = axes[0][1]
    plot_prop(ax, xx,yy, Rho,prop_name="Density", prop_unit="kg/m$^\mathregular{3}$")
    # 3. H
    ax = axes[0][2]
    plot_prop(ax, xx,yy, T-273.15,prop_name="Temperature", prop_unit="$^{\circ}$C",cmap='GnBu')
    # 3. Cp
    ax = axes[0][3]
    plot_prop(ax, xx,yy, Cp,prop_name="Isobaric specific heat", prop_unit="J/kg/K",cmap='YlOrRd')
    # 1. Mu_l
    ax=axes[1][0]
    plot_prop(ax, xx,yy, (Mu_l)*1E6,prop_name="log$_{10}$ $\mu_l$", prop_unit="$\mu Pa \cdot s$",cmap='GnBu')
    # ax.grid()
    # 2. Mu_v
    ax=axes[1][1]
    plot_prop(ax, xx,yy, (Mu_v)*1E6,prop_name="log$_{10}$ $\mu_v$", prop_unit="$\mu Pa \cdot s$",cmap='BuGn')

    # set axis
    for ax in axes[1,:]: ax.set_xlabel('Bulk salinity (wt.% NaCl)')
    for ax in axes[:,0]: ax.set_ylabel('Pressure (bar)')

    savefig('Props_H%.0fMJ_kg'%(H0/1E6))

def plot_props_P0(P0, sw):
    H = np.linspace(0.1, 5.5, 400)*1E6
    X = np.linspace(1E-6,1,400)
    HH,XX = np.meshgrid(H,X)
    PP = HH*0 + P0
    phase = np.zeros_like(HH)
    Rho = np.zeros_like(HH)
    T = np.zeros_like(HH)
    Cp = np.zeros_like(HH)
    Mu_l = np.zeros_like(HH)
    Mu_v = np.zeros_like(HH)
    for i in range(0,HH.shape[0]):
        for j in range(0,HH.shape[1]):
            props = sw.UpdateState_HPX(HH[i][j], PP[i][j], XX[i][j])
            phase[i][j] = props.phase
            Rho[i][j] = props.Rho
            Cp[i][j] = props.Cp
            T[i][j] = props.T
            Mu_l[i][j] = props.Mu_l
            Mu_v[i][j] = props.Mu_v
    T[T<273.15]=np.nan
    # plot
    fig, axes=plt.subplots(2,4,figsize=(28,12),sharey=True, sharex=True,gridspec_kw={'wspace':0.1,'hspace':0.2})

    xx,yy = XX*100, HH/1E6
    # 1. phase
    ax = axes[0][0]
    plot_Phase(ax, xx,yy, phase, sw)
    # 2. Rho
    ax = axes[0][1]
    plot_prop(ax, xx,yy, Rho,prop_name="Density", prop_unit="kg/m$^\mathregular{3}$")
    # 3. H
    ax = axes[0][2]
    plot_prop(ax, xx,yy, T - 273.15,prop_name="Temperature", prop_unit="$^{\circ}$C",cmap='GnBu')
    # 3. Cp
    ax = axes[0][3]
    plot_prop(ax, xx,yy, Cp,prop_name="Isobaric specific heat", prop_unit="J/kg/K",cmap='YlOrRd')
    # 1. Mu_l
    ax=axes[1][0]
    plot_prop(ax, xx,yy, (Mu_l)*1E6,prop_name="log$_{10}$ $\mu_l$", prop_unit="$\mu Pa \cdot s$",cmap='GnBu')
    # ax.grid()
    # 2. Mu_v
    ax=axes[1][1]
    plot_prop(ax, xx,yy, (Mu_v)*1E6,prop_name="log$_{10}$ $\mu_v$", prop_unit="$\mu Pa \cdot s$",cmap='BuGn')

    # set axis
    for ax in axes[1,:]: ax.set_xlabel('Bulk salinity (wt.% NaCl)')
    for ax in axes[:,0]: ax.set_ylabel('Specific enthalpy (MJ/kg)')

    savefig('Props_P%.0fbar'%(P0/1E5))

def plot_props_X0(X0,sw):
    H = np.linspace(0.1, 5.5, 400)*1E6
    P = np.linspace(1E5, 1000E5, 400)
    HH,PP = np.meshgrid(H,P)
    XX = HH*0 + X0
    phase = np.zeros_like(HH)
    Rho = np.zeros_like(HH)
    T = np.zeros_like(HH)
    Cp = np.zeros_like(HH)
    Mu_l = np.zeros_like(HH)
    Mu_v = np.zeros_like(HH)
    for i in range(0,HH.shape[0]):
        for j in range(0,HH.shape[1]):
            props = sw.UpdateState_HPX(HH[i][j], PP[i][j], XX[i][j])
            phase[i][j] = props.phase
            Rho[i][j] = props.Rho
            Cp[i][j] = props.Cp
            T[i][j] = props.T
            Mu_l[i][j] = props.Mu_l
            Mu_v[i][j] = props.Mu_v
    T[T<273.15]=np.nan
    # plot
    fig, axes=plt.subplots(2,4,figsize=(28,12),sharey=True, sharex=True,gridspec_kw={'wspace':0.1,'hspace':0.2})

    xx,yy = HH/1E6,PP/1E5
    # 1. phase
    ax = axes[0][0]
    plot_Phase(ax, xx,yy, phase, sw)
    # 2. Rho
    ax = axes[0][1]
    plot_prop(ax, xx,yy, Rho,prop_name="Density", prop_unit="kg/m$^\mathregular{3}$")
    # 3. H
    ax = axes[0][2]
    plot_prop(ax, xx,yy, T-273.15,prop_name="Temperature", prop_unit="$^{\circ}$C",cmap='GnBu')
    # 3. Cp
    ax = axes[0][3]
    plot_prop(ax, xx,yy, Cp,prop_name="Isobaric specific heat", prop_unit="J/kg/K",cmap='YlOrRd')
    # 1. Mu_l
    ax=axes[1][0]
    plot_prop(ax, xx,yy, (Mu_l)*1E6,prop_name="log$_{10}$ $\mu_l$", prop_unit="$\mu Pa \cdot s$",cmap='GnBu')
    # ax.grid()
    # 2. Mu_v
    ax=axes[1][1]
    plot_prop(ax, xx,yy, (Mu_v)*1E6,prop_name="log$_{10}$ $\mu_v$", prop_unit="$\mu Pa \cdot s$",cmap='BuGn')

    # set axis
    for ax in axes[1,:]: ax.set_xlabel('Specific enthalpy (MJ/kg)')
    for ax in axes[:,0]: ax.set_ylabel('Pressure (bar)')

    # ax.set_xscale('log')
    savefig('Props_X%.0fwt'%(X0*100))

# plot_props_P0X0(sw_84)

# %%
# Isobaric section
# --------------------------
plot_props_P0(100E5, sw_84)

# %%
# Const enthalpy
# --------------------------
plot_props_H0(2E6, sw_84)

# %%
# Constant X
# --------------------------
plot_props_X0(0.1, sw_84)