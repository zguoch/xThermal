# -*- coding: utf-8 -*-
"""
0. Phase diagram
==============================
"""

# Some python packages for data visualization
import numpy as np 
import time
import copy
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import patches
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator
mpl.rcParams['font.family'] = 'Arial'  #default font family
mpl.rcParams['mathtext.fontset'] = 'cm' #font for math
dpi=100
fmt_figs=['pdf'] #['svg','pdf']
result_path='.'
figpath=result_path
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

# Calculate
def cal_phase(water):
    H = np.linspace(0.2,4,200)*1E6
    p = np.linspace(1E5,500e5,200)
    # p = np.zeros_like(H)
    # for i in range(0,len(H)): p[i] = water.Boiling_p(T[i])
    # calculate phase index
    HH, pp = np.meshgrid(H,p)
    phase = np.zeros_like(HH)
    for i in range(0,HH.shape[0]):
        for j in range(0,HH.shape[1]):
            # props=water.UpdateState_HPX(HH[i][j], pp[i][j])
            # phase[i][j]=props.phase
            phase[i][j] = water.findPhaseRegion_HPX(HH[i][j], pp[i][j])
    phase_unique = np.sort(np.unique(phase))
    phase_name = ['']*len(phase_unique)
    for i,phase0 in enumerate(phase_unique): 
        phase[phase==phase0]=i+phase_unique.max()+10
        phase_name[i]=water.phase_name(int(phase0))
    return H,p,HH,pp,phase,phase_name

# Plot both linear and log scale
def phaseDiagram(water, axes=None):
    T,p,HH,pp,phase,phase_name = cal_phase(water)
    pp,HH = pp/1E5, HH/1E6
    if(axes==None): 
        # fig,axes=plt.subplots(1,2,figsize=(15,7),gridspec_kw={'wspace':0.02},dpi=dpi)
        fig=plt.figure(figsize=(8,7))
        ax=plt.gca()
    ax.set_ylabel('Pressure (bar)')

    # plot
    cmap = plt.get_cmap("Dark2")
    # customize cmap
    colors=list(copy.deepcopy(cmap.colors))
    colors[0:8]=['lightblue','red','lightgreen','lightgray','violet','yellow','lightcyan','m']
    cmap.colors=tuple(colors)

    ax.set_ylim(pp.min(),pp.max())
    ax.set_xlim(HH.min(),HH.max())
    ax.text(0.98,0.98,water.name(),ha='right',va='top',bbox={'fc':'w','ec':'gray'}, transform=ax.transAxes)
    ax.set_xlabel('Specific enthalpy (MJ/kg)')
    # ax.plot(T-273.15, p/1E5, label='Boiling curve')
    # ax.plot(water.T_critical()-273.15, water.p_critical()/1E5,'o',mfc='r',mec='w',label='Critical point')
    CS=ax.contourf(HH,pp,phase, cmap=cmap,vmin=phase.min()-0.5, vmax=phase.max()+0.5, levels=np.linspace(phase.min()-0.5,phase.max()+0.5,len(phase_name)+1))
    ax_cb = ax.inset_axes([0,1.03,1,0.05])
    cb=plt.colorbar(CS, cax=ax_cb, orientation='horizontal',ticklocation='top',ticks=np.arange(phase.min(),phase.max()+1))
    cb.ax.set_xticklabels(phase_name)
    # if(ax==axes[1]):
    #     ax.yaxis.set_ticks_position('right')
    #     ax.xaxis.set_minor_locator(MultipleLocator(20))
    #     ax.grid(which='major',color='gray')
    #     ax.grid(which='minor',color='lightgray')
    # ax.legend(loc='lower right')
    savefig('phase_%s'%(water.name()))


# %%
# IAPS84 Phase diagram
# -------------------------
# phaseDiagram(iaps84)

# %%
# IAPWS95 Phase diagram: build in xThermo
# --------------------------------------------------
phaseDiagram(iapws95)

# %%
# IAPWS95 Phase diagram: CoolProp
# --------------------------------------------------

# phaseDiagram(iapws95_CoolProp)