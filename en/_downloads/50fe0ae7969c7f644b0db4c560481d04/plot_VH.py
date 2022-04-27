# -*- coding: utf-8 -*-
"""
5. Vapor + Halite coexistence surface
==================================================
.. include:: /include.rst_
Calculate VH surface and properties, and compare result with :cite:`Driesner2007Part2`.
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
fmt_figs = []  # ['svg','pdf']
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

def plot_3d(sw,scale='linear'):
    fig=plt.figure(figsize=(14,14))
    ax = fig.add_subplot(111,projection='3d',facecolor='None')
    xcenter=1 # X<1% wt.% use log scale, if X>=1% wt.% NaCl, use linear scale
    axtrans=[]
    if(scale=='loglinear'):
        axtrans=helpfunc.set_axis_diagram_3D_loglinearx(ax,xcenter=xcenter,ratio_log_lin=(1,1),zlim=(sw.pmin()/1E5,600),xMajor_loc_log=2,xMajor_loc_linear=10,xMinor_loc_linear=2,zMajor_loc=100,zMinor_loc=20,ylim=(90,820),xlim=(1E-14,100),xlabel='log$_{\mathregular{10}}$(Wt.% NaCl)')
    else:
        helpfunc.set_axis_diagram_3D(ax,zlim=(sw.pmin()/1E5,600),zMajor_loc=100,zMinor_loc=20)
    transform_X = lambda X : helpfunc.data2axis_loglin(axtrans,X) if(scale=='loglinear') else X
    # plot VLH surface as reference
    T = np.linspace(H2ONaCl.T_MIN_VLH, sw.Tmax_VLH(), 100)
    P_vlh = np.array(sw.P_VLH(T))
    # X_haliteLiquidus=np.array(sw.X_HaliteLiquidus(T,P))
    Xl_vlh, Xv_vlh = np.array(sw.X_VLH(T,P_vlh))
    # plot
    n_log,n_linear=20,40
    X = np.linspace(H2ONaCl.X_MIN, H2ONaCl.X_MAX, n_log+n_linear)
    TT,XX = np.meshgrid(T, X)
    PP = np.zeros_like(TT)
    for i in range(0,PP.shape[0]):
        PP[i,:]=P_vlh
    # 1. vapor -> liquid region
    for j in range(0,PP.shape[1]):
        XX[:,j] = np.append(10**np.linspace(np.log10(Xv_vlh[j]),np.log10(xcenter/100),n_log), np.linspace(xcenter/100,Xl_vlh[j],n_linear))
    xx_plot = transform_X(XX*100)
    ax.plot_surface(xx_plot,TT-273.15,PP/1E5,color='w',alpha=0.5)
    ax.plot_wireframe(xx_plot,TT-273.15,PP/1E5,ec='b',lw=0.8,label='VLH: V->L')
    # 2. liquid -> halite region
    for j in range(0,PP.shape[1]):
        XX[:,j] = np.linspace(Xl_vlh[j], 1,len(X))
    xx_plot = transform_X(XX*100)
    ax.plot_surface(xx_plot,TT-273.15,PP/1E5,color='w',alpha=0.5)
    ax.plot_wireframe(xx_plot,TT-273.15,PP/1E5,ec='orange',lw=0.8,label='VLH: L->H')
    # 3. Halite saturated vapor in VH region
    P = np.linspace(0,1,40)
    TT,PP = np.meshgrid(T,P)
    for i in range(0,len(T)):
        PP[:,i] = np.linspace(sw.pmin(),P_vlh[i],len(P))
    XX_VH = np.array(sw.X_VH(TT.reshape(-1,),PP.reshape(-1,)))
    xx_plot = transform_X(XX_VH.reshape(TT.shape)*100)
    ax.plot_surface(xx_plot,TT-273.15,PP/1E5,color='w',alpha=0.5)
    ax.plot_wireframe(xx_plot,TT-273.15,PP/1E5,ec='purple',lw=0.8,label='VH: vapor')
    # halite
    ax.plot(transform_X(Xl_vlh*0+100), T-273.15, P_vlh/1E5, color='k',label='VLH: halite',lw=4)
    # halite liquidus
    ax.plot(transform_X(Xl_vlh*100), T-273.15, P_vlh/1E5, color='g',label='VLH: liquid',lw=4)
    # upper pressure bound
    ax.plot(transform_X(Xv_vlh*100), T-273.15, P_vlh/1E5, color='r',label='VLH: vapor',lw=4)
    # lower pressure bound
    Xv_pmin = np.array(sw.X_VH(T,T*0 + sw.pmin()))
    ax.plot(transform_X(Xv_pmin*100), T-273.15, T*0+sw.pmin()/1E5, color='blue',label='VH: vapor (P=%.0f bar)'%(sw.pmin()/1E5),lw=4)
    # plot several isotherm profile on X-p plane: Fig. 8a in Driesner(2007a)
    for i,T0 in enumerate([300, 450, 500, 550, 650]):
        P_ = np.linspace(sw.pmin(),sw.P_VLH(T0+273.15),100)
        X_ = np.array(sw.X_VH(P_*0 + T0 + 273.15, P_))
        l,=ax.plot(transform_X(X_*100), X_*0 + ax.get_ylim()[1], P_/1E5)
        helpfunc.text3d(ax, (transform_X(10**(-13+1.2*i)),800,50+i*50), 'T=%.0f$^{\circ}$C'%(T0),zdir='y',size=0.05, fc=l.get_color(), ec="None",ha='left',va='bottom')
        # # add line of VLH
        # Xl_=sw.X_HaliteLiquidus(T0+273.15, P_[-1])
        # ax.plot(np.log10(np.array([X_[-1],Xl_])*100),np.array([ax.get_ylim()[1]]*2),np.array([P_[-1]]*2)/1E5,color='b')
        # ax.plot(np.log10(np.array([Xl_,1])*100),np.array([ax.get_ylim()[1]]*2),np.array([P_[-1]]*2)/1E5,color='orange')
    # ax.plot(np.log10(Xv_vlh*100), T*0 + ax.get_ylim()[1], P_vlh/1E5, color='r',label='VLH: vapor',lw=1,ls='dashed')
    # legend
    leg=ax.legend()
    # change legend handle of wireframe of phase boundaries to wireframe hatch
    for i in [0,1,2]:
        leg.legendHandles[i]=Patch(facecolor='white', edgecolor=leg.legendHandles[i]._color,linewidth=0.0,label=leg.texts[i]._text,hatch='++++')
    ax.legend(handles=leg.legendHandles, loc='upper left',ncol=7)
    # text
    helpfunc.text3d(ax,(transform_X(1E-12),800,sw.pmin()/1E5),"P$_{min}$=%.0f bar"%(sw.pmin()/1E5),size=0.07,angle=-90,ec='None',fc='k')
    ax.view_init(elev=25, azim=-145)
    savefig('PhaseBoundary_VH_3D_%s'%(scale))

# %%
# General 3D view: linear scale (left panel) and log-linear scale (right panel)
# ------------------------------------------------------------------------------------------------
# The VH surface connect to the Vapor side of VLH surface which can be calculated by function :code:`Xl_vlh, Xv_vlh = V_VLH(T, P)`
plot_3d(sw_84,scale='linear')
plot_3d(sw_84,scale='loglinear')

# %%
# Benchmark comparison
# --------------------------
# Compare halite saturated vapor composition of |xThermo| and :cite:`Driesner2007Part2`, and also compare result based on different EOS of water.
def benchmark_VL(sw,mmc4='../Driesner2007a/1-s2.0-S0016703707002943-mmc4.txt'):
    # compare
    if(not os.path.exists(mmc4)):
        print('Please set correct mmc4 file path: %s'%(mmc4))
        exit()
    T0,P0,X0=data=np.loadtxt(mmc4, skiprows=7,unpack=True)
    XL0_wt= np.array(sw.Mol2Wt(X0))*100
    # 1. calculate halite saturated vapor composition: XV_VH
    X_ = np.array(sw.X_VH(T0+273.15,P0*1E5))
    X_mol_ = np.array(sw.Wt2Mol(X_))
    # compare result dict
    Data0 = {'Xv':X0}
    Data_ = {'Xv':X_mol_}
    Err,RErr={},{}
    for key in Data0.keys(): Err[key],RErr[key] = Data0[key]-Data_[key], np.abs(Data0[key]-Data_[key])/(Data0[key])*100.0
    # print to file
    fpout = open('%s/mmc4a_%s.csv'%(result_path,sw.name_backend()),'w')
    fpout.write('T[C],P[bar],XV(Driesner)[mol],XV(xThermo)[mol],XV(diff)[mol]\n')
    for i in range(0,len(T0)):
        fpout.write('%.6e,%.6e'%(T0[i]-273.15,P0[i]))
        for key in Data0.keys():
            fpout.write(',%.6e,%.6e,%.6e'%(Data0[key][i], Data_[key][i],compare(Data0[key][i],Data_[key][i])))
        fpout.write('\n')
    fpout.close()
    # plot data
    fig=plt.figure(figsize=(8,8))
    ax=plt.gca()
    cmap,vmin,vmax='YlGnBu_r',sw.pmin()/1E5,H2ONaCl.P_Peak_VLH/1E5
    CS=ax.scatter(X0,T0,c=P0,cmap=cmap,vmin=vmin,vmax=vmax,s=5)
    ax_cb=ax.inset_axes([1.02,0,0.05,1])
    plt.colorbar(CS,cax=ax_cb,label='Pressure (bar)')
    # also plot VLH: vapor curve
    T = np.linspace(H2ONaCl.T_MIN_VLH, sw.Tmax_VLH(), 100)
    P_vlh = np.array(sw.P_VLH(T))
    # X_haliteLiquidus=np.array(sw.X_HaliteLiquidus(T,P))
    Xl_vlh, Xv_vlh = np.array(sw.X_VLH(T,P_vlh))
    helpfunc.plot_coloredline(ax,np.array(sw.Wt2Mol(Xv_vlh)), T-273.15, P_vlh/1E5,cmap=cmap,vmin=vmin,vmax=vmax,lw=2)
    XV_VH_pmin = np.array(sw.X_VH(T,T*0+sw.pmin()))
    helpfunc.plot_coloredline(ax,np.array(sw.Wt2Mol(XV_VH_pmin)), T-273.15, T*0+sw.pmin()/1E5 ,cmap=cmap,vmin=vmin,vmax=vmax,lw=2)
    ax.annotate("VLH: vapor",ha='center',va='center', xy=(1E-6,320),xytext=(1E-6,150), bbox={'fc':'None','ec':'r'},fontsize=14,fontweight='bold', arrowprops=dict(arrowstyle="->",connectionstyle="arc3"),)
    ax.annotate("Halite saturated vapor\nat $P_{min}$=%.0f bar"%(sw.pmin()/1E5),ha='center',va='top', xy=(2E-13,290),xytext=(2E-13,500), bbox={'fc':'None','ec':'b'},fontsize=14,fontweight='bold', arrowprops=dict(arrowstyle="->",connectionstyle="arc3"),)
    ax.grid(lw=0.2,color='gray')
    ax.set_ylim(H2ONaCl.T_MIN_VLH-273.15, H2ONaCl.T_MAX_VLH-273.15)
    ax.set_xscale('log')
    ax.set_ylabel('Temperature ($^{\circ}$C)')
    ax.set_xlabel('X$_{\mathregular{NaCl}}$ (mole fraction)')

    # savefig('mmc4a')
    # statistics of the difference
    table=[]
    for key,name in zip(list(Err.keys()),['XV (mole fraction)']):
        RErr[key] = RErr[key][~(np.isnan(RErr[key]) | np.isinf(RErr[key]))]
        table.append([name,Err[key].min(),Err[key].max(),RErr[key].min(),RErr[key].max()])
    print(tabulate(table, headers=['Critical property', 'Err. Min', 'Err. Max','RE. Min(%)','RE. Max(%)'],numalign="right",floatfmt=".6f"))

# %%
# Based on IAPS84 EOS
# ^^^^^^^^^^^^^^^^^^^^^^^
benchmark_VL(sw_84)

# %%
# Based on IAPWS95 EOS
# ^^^^^^^^^^^^^^^^^^^^^^^
benchmark_VL(sw_95)

# %%
# Result table
# -----------------------------
# .. seealso::
#
#     |mmc4| in :cite:`Driesner2007Part1` and Fig. 8 in :cite:`Driesner2007Part1`.
#
# Result data calculated by |xThermo| based on both water EOS of IAPS84 and IAPWS95.
#
# .. tab:: IAPS84
#
#     .. csv-table:: Comparison between result of :cite:`Driesner2007Part2` and result calculated by |xThermo|.
#         :file: mmc4a_IAPS84.csv
#         :header-rows: 1
#
# .. tab:: IAPWS95
#
#     .. csv-table:: Comparison between result of :cite:`Driesner2007Part2` and result calculated by |xThermo|.
#         :file: mmc4a_IAPWS95.csv
#         :header-rows: 1
#
# .. tip::
#
#     The help function for 3D plot can be downloaded at here: :download:`helpfunc.py`
