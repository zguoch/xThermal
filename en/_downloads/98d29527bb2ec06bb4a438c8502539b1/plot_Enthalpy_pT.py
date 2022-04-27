# -*- coding: utf-8 -*-
"""
2. Specific enthalpy
=========================
"""

import numpy as np 
import time
import linecache
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import patches
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator
mpl.rcParams['font.family'] = 'Arial'  #default font family
mpl.rcParams['mathtext.fontset'] = 'cm' #font for math
fmt_figs=['pdf'] #['svg','pdf']
figpath='.'
def savefig(figname):
    for fmt_fig in fmt_figs:
        figname_full = '%s/%s.%s'%(figpath,figname,fmt_fig)
        plt.savefig(figname_full, bbox_inches='tight')
        print('figure saved: ',figname_full)

# Import package of xThermo
from xThermo import H2O
iaps84 = H2O.cIAPS84()
iapws95_CoolProp = H2O.cIAPWS95_CoolProp()
iapws95 = H2O.cIAPWS95()

# Calculate and plot result
def plot_prop(water, ax=None,name_prop='Enthalpy',unit_prop='MJ/kg',T=[],p=[],cmap='GnBu'):
    if(len(T)==0): T = np.linspace(274,1273,150)
    if(len(p)==0): p = np.linspace(1E5,2000e5,150)
    TT,PP = np.meshgrid(T,p)
    prop = np.zeros_like(TT)
    for i in range(0,TT.shape[0]):
        for j in range(0,TT.shape[1]):
            props = water.UpdateState_TPX(TT[i][j], PP[i][j])
            prop[i][j] = props.H
    # plot
    isSaveFig=False
    if(ax is None): isSaveFig = True
    if(ax==None):
        fig=plt.figure(figsize=(7,7))
        ax=plt.gca()
    CS = ax.contourf(TT-273.15,PP/1E5, prop/1E6, levels=50, cmap=cmap)
    ax_cb = ax.inset_axes([1.01, 0, 0.05, 1])
    plt.colorbar(CS, cax=ax_cb, orientation='vertical',label='%s (%s)'%(name_prop,unit_prop))
    # labels
    ax.text(0.98,0.98,water.name(),ha='right',va='top',bbox={'fc':'w','ec':'gray'}, transform=ax.transAxes)
    ax.set_xlabel('Temperature ($^{\circ}$C)')
    ax.set_ylabel('Pressure (bar)')
    if(isSaveFig): savefig('H2O_%s_%s'%(name_prop, water.name()))
    return TT-273.15,PP/1E5,prop

def plot_error(ax, XX,YY,ZZ, eosA, eosB,cmap='RdBu',unit='J/kg'):
    CS=ax.contourf(XX,YY,ZZ,levels=50,cmap=cmap, norm=mpl.colors.CenteredNorm())
    ax_cb = ax.inset_axes([1.01, 0, 0.05, 1])
    plt.colorbar(CS, cax=ax_cb, orientation='vertical',label='Diff (%s)'%(unit))
    ax.text(0.98,0.98,'%s-%s'%(eosA.name(),eosB.name()),ha='right',va='top',bbox={'fc':'w','ec':'gray'}, transform=ax.transAxes)
    ax.set_xlabel('Temperature ($^{\circ}$C)')

# %%
# IAPS84 EOS
# -------------------------
TT,PP,prop=plot_prop(iaps84)

# %%
# Comparison: IAPS84 and IAPWS95
# --------------------------------------------------

fig, axes = plt.subplots(1,5,figsize=(30,5),gridspec_kw={'wspace':0.35})
TT,PP,prop84 = plot_prop(iaps84,ax=axes[0])
TT,PP,prop95 = plot_prop(iapws95,ax=axes[1])
TT,PP,prop95_CoolProp = plot_prop(iapws95_CoolProp,ax=axes[2])
plot_error(axes[3],TT,PP,prop95 - prop84,eosA=iapws95, eosB=iaps84)
plot_error(axes[4],TT,PP,prop95 - prop95_CoolProp,eosA=iapws95, eosB=iapws95_CoolProp)
for ax in axes[1:]: ax.set_ylabel(None)
savefig('diff_84_95')

# %%
# Comparison in region close to critical point: IAPWS95 and IAPWS95
# ---------------------------------------------------------------------------
    
fig, axes = plt.subplots(1,3,figsize=(18,5),gridspec_kw={'wspace':0.35})
T = np.linspace(iaps84.T_critical()-10, iaps84.T_critical()+10, 100)
p = np.linspace(iaps84.p_critical()-10E5, iaps84.p_critical()+10E5, 100)
TT,PP,prop84 = plot_prop(iaps84,ax=axes[0], T=T, p=p)
TT,PP,prop95 = plot_prop(iapws95,ax=axes[1], T=T, p=p)
plot_error(axes[2],TT,PP,prop95 - prop84, eosA=iapws95, eosB=iaps84)
for ax in axes[1:]: ax.set_ylabel(None)
savefig('diff')
print('Critical density, IAPS84: %.4f kg/m3, IAPWS95: %.4f kg/m3, diff: %.4f kg/m3'%(iaps84.rhomass_critical(),iapws95.rhomass_critical(), iaps84.rhomass_critical()-iapws95.rhomass_critical()))