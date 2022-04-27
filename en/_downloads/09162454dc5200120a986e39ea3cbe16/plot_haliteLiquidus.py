# -*- coding: utf-8 -*-
"""
2. Halite liquidus
=========================
.. include:: /include.rst_
Calculate halite liquidus and compare result with :cite:`Driesner2007Part2`.
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
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
# 3d plot
import helpfunc
mpl.rcParams['font.family'] = 'Arial'  # default font family
mpl.rcParams['mathtext.fontset'] = 'cm'  # font for math
fmt_figs = []  # ['svg','pdf']
figpath = '.'
result_path='../../../gallery_H2ONaCl/pT'
def savefig(figname):
    for fmt_fig in fmt_figs:
        figname_full = '%s/%s.%s' % (figpath, figname, fmt_fig)
        plt.savefig(figname_full, bbox_inches='tight')
        print('figure saved: ', figname_full)
compare = lambda a,b : float(str('%.6e'%(a)))-float(str('%.6e'%(b)))
# Import package of xThermo
from xThermo import H2O
from xThermo import NaCl
from xThermo import H2ONaCl
sw_84 = H2ONaCl.cH2ONaCl("IAPS84")
sw_95 = H2ONaCl.cH2ONaCl("IAPWS95")

def plot_3d(sw,style='linear'):
    fig=plt.figure(figsize=(14,14))
    ax = fig.add_subplot(111,projection='3d',facecolor='None')
    helpfunc.set_axis_diagram_3D(ax)
    # plot
    if(style=='linear'):
        pb = sw.PhaseBoundary_HaliteLiquidus_DeformLinear()
        TT,PP,XX = np.array(pb.T), np.array(pb.p), np.array(pb.X)
        ax.plot_wireframe(XX*100,TT-273.15,PP/1E5, color='green',lw=0.5,label='Halite liquidus')
        ax.plot_wireframe(ax.get_xlim()[1],TT-273.15,PP/1E5, color='orange',lw=0.3,label='Halite liquidus: p-X plane')
        ax.plot_wireframe(XX*100,ax.get_ylim()[1],PP/1E5, color='cyan',lw=0.3,label='Halite liquidus: p-T plane')
    elif(style=='triangle'):
        # halite liquidus surface mesh
        pb = sw.PhaseBoundary_HaliteLiquidus()
        T,p,X,connection = np.array(pb.x),np.array(pb.y),np.array(pb.z), np.array(pb.connection)
        ax.plot_trisurf(X*100, T-273.15, connection,p/1E5, color='green',lw=0.1,label='Halite liquidus')
        ax.plot_trisurf(X*100, T*0 + ax.get_ylim()[1], connection,p/1E5,color='orange',lw=0.1)
        ax.plot_trisurf(X*0 + ax.get_xlim()[1], T-273.15, connection,p/1E5,color='cyan',lw=0.1)
    # ax.legend(ncol=3)
    savefig('HaliteLiquidus_3D_%s'%(style))

# %%
# General 3D view
# --------------------------------
plot_3d(sw_84)
# plot_3d('triangle')

# %%
# Benchmark comparison
# --------------------------
# Compare result of |xThermo| and :cite:`Driesner2007Part2`, and also compare result based on different EOS of water.
def plot_props(TT,PP,XX,rho_,h_,unit_H='MJ/kg'):
    wspace=0.1
    fig,axes=plt.subplots(1,4,figsize=(21,5),gridspec_kw={'wspace':wspace}, sharey=True)
    cmap_rho, cmap_h = 'YlGnBu_r','rainbow'
    # 1. Rho
    axes[0].contourf(TT-273.15,PP/1E5,rho_.reshape(TT.shape),levels=50,cmap=cmap_rho)
    CS=axes[1].contourf(XX*100,PP/1E5,rho_.reshape(TT.shape),levels=50,cmap=cmap_rho)
    ax_cb_rho = axes[0].inset_axes([0,1.03,2+wspace,0.05])
    plt.colorbar(CS,cax=ax_cb_rho,orientation='horizontal',label='Density (kg/m$^{\mathregular{3}}$)')
    # 2. H
    scale_H = 1
    if(unit_H=='MJ/kg'): scale_H=1E-6
    axes[2].contourf(TT-273.15,PP/1E5,h_.reshape(TT.shape)*scale_H,levels=50,cmap=cmap_h)
    CS=axes[3].contourf(XX*100,PP/1E5,h_.reshape(TT.shape)*scale_H,levels=50,cmap=cmap_h)
    ax_cb_h = axes[2].inset_axes([0,1.03,2+wspace,0.05])
    plt.colorbar(CS,cax=ax_cb_h,orientation='horizontal',label='Specific enthalpy (%s)'%(unit_H))
    for ax_cb in [ax_cb_rho,ax_cb_h]:
        ax_cb.xaxis.set_ticks_position('top')
        ax_cb.xaxis.set_label_position('top')
    for ax in [axes[0],axes[2]]: ax.set_xlabel('Temperature ($^{\circ}$C)')
    for ax in [axes[1],axes[3]]: ax.set_xlabel('X$_{\mathregular{NaCl}}$ (wt.% NaCl)')
    axes[0].set_ylabel('Pressure (bar)')
def plot_props_HaliteLiquidus(sw):
    # calculate
    pb = sw.PhaseBoundary_HaliteLiquidus_DeformLinear()
    TT,PP,XX = np.array(pb.T), np.array(pb.p), np.array(pb.X)
    rho_ = np.array(sw.Rho_phase(TT.reshape(-1,), PP.reshape(-1,), XX.reshape(-1,), H2ONaCl.Liquid))
    h_ = np.array(sw.H_phase(TT.reshape(-1,), PP.reshape(-1,), XX.reshape(-1,), H2ONaCl.Liquid))
    # plot
    plot_props(TT,PP,XX,rho_,h_)
    savefig('HaliteLiquidus_props_%s'%(sw.name_backend()))
    return TT,PP,XX,rho_,h_
def plot_err(ax,x,y,data,label='',cmap='rainbow'):
    # plot difference between xThermo and Driesner(2007b)
    CS=ax.scatter(x, y,c=data,cmap=cmap)
    ax_cb=ax.inset_axes([0,1.02,1,0.05])
    plt.colorbar(CS,cax=ax_cb,label=label,orientation='horizontal')
    ax_cb.xaxis.set_label_position('top')
    ax_cb.xaxis.set_ticks_position('top')
def benchmark_HaliteLiquidus(sw,mmc1='../Driesner2007b/1-s2.0-S0016703707002955-mmc1.txt'):
    # compare
    if(not os.path.exists(mmc1)):
        print('Please set correct mmc1 file path: %s'%(mmc1))
        exit()
    data=np.loadtxt(mmc1, skiprows=7)
    T0,P0,X0,rho0,h0=data[:,0]+273.15,data[:,1]*1E5,data[:,2],data[:,3],data[:,4]
    # 1. calculate halite liquidus
    X_=sw.X_HaliteLiquidus(T0,P0)
    X_mol = sw.Wt2Mol(X_)
    # 2. calculate liquid density on halite liquidus
    rho_ = np.array(sw.Rho_phase(T0, P0, sw.Mol2Wt(X0), H2ONaCl.Liquid))
    h_ = np.array(sw.H_phase(T0, P0, sw.Mol2Wt(X0), H2ONaCl.Liquid))
    # compare result dict
    Data0 = {'X':X0,'rho':rho0,'h':h0}
    Data_={'X':X_mol,'rho':rho_,'h':h_}
    Err,RErr={},{}
    for key in Data0.keys(): Err[key],RErr[key] = Data0[key]-Data_[key], np.abs(Data0[key]-Data_[key])/(Data0[key])*100.0
    # print to file
    fpout = open('%s/mmc1_%s.csv'%(result_path,sw.name_backend()),'w')
    fpout.write('T[C],P[bar],X(Driesner)[mol],X(xThermo)[mol],X(diff)[mol],Rho(Driesner)[kg/m3],Rho(xThermo),Rho(err),H(Driesner)[J/kg],H(xThermo),H(err)\n')
    for i in range(0,len(T0)):
        fpout.write('%.6e,%.6e'%(T0[i]-273.15,P0[i]))
        for key in Data0.keys():
            fpout.write(',%.6e,%.6e,%.6e'%(Data0[key][i], Data_[key][i],compare(Data0[key][i],Data_[key][i])))
        fpout.write('\n')
    fpout.close()
    # plot difference
    fig,axes=plt.subplots(1,3,figsize=(21,5),sharey=True,gridspec_kw={'wspace':0.1})
    plot_err(axes[0], T0, P0/1E5, RErr['X'],'Relative difference: Composition (%)')
    plot_err(axes[1], T0, P0/1E5, RErr['rho'],'Relative difference: Density (%)')
    plot_err(axes[2], T0, P0/1E5, RErr['h'],'Relative difference: Specific enthalpy (%)')
    for ax in axes: ax.set_xlabel('Temperature ($^{\circ}$C)')
    axes[0].set_ylabel('Pressure (bar)')
    # statistics of the difference
    table=[]
    for key,name in zip(list(Err.keys()),['Composition (mole fraction)','Density (kg/m3)','Specific enthalpy (J/kg)']):
        RErr[key] = RErr[key][~(np.isnan(RErr[key]) | np.isinf(RErr[key]))]
        table.append([name,Err[key].min(),Err[key].max(),RErr[key].min(),RErr[key].max()])
    print(tabulate(table, headers=['Critical property', 'Err. Min', 'Err. Max','RE. Min(%)','RE. Max(%)'],numalign="right",floatfmt=".6f"))

# %%
# Based on IAPS84 EOS
# ^^^^^^^^^^^^^^^^^^^^^^^
TT,PP,XX,rho_84,h_84=plot_props_HaliteLiquidus(sw_84)
benchmark_HaliteLiquidus(sw_84)

# %%
# Based on IAPWS95 EOS
# ^^^^^^^^^^^^^^^^^^^^^^^
# TT,PP,XX,rho_95,h_95=plot_props_HaliteLiquidus(sw_95)
# benchmark_HaliteLiquidus(sw_95)

# %%
# Difference between IAPS84 and IAPWS95 EOS
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# plot_props(TT,PP,XX, rho_84-rho_95, (h_84-h_95),unit_H='J/kg')

# %%
# Result table
# -----------------------------
#
# Result data calculated by |xThermo| based on both water EOS of IAPS84 and IAPWS95.
#
# .. seealso::
#
#     |mmc1| in :cite:`Driesner2007Part2` and Fig. 1,5,6 in :cite:`Driesner2007Part1`.
#
# .. tab:: IAPS84
#
#     .. csv-table:: Comparison between result of :cite:`Driesner2007Part2` and result calculated by |xThermo|.
#         :file: mmc1_IAPS84.csv
#         :header-rows: 1
#
# .. tab:: IAPWS95
#
#     .. csv-table:: Comparison between result of :cite:`Driesner2007Part2` and result calculated by |xThermo|.
#         :file: mmc1_IAPWS95.csv
#         :header-rows: 1
#
#
# .. tip::
#
#     The help function for 3D plot can be downloaded at here: :download:`helpfunc.py`
