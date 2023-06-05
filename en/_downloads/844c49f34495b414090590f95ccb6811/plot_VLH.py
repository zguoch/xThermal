# -*- coding: utf-8 -*-
"""
3. Vapor + Liquid + Halite coexistence surface
==================================================
.. include:: /include.rst_
Calculate VLH surface and properties, and compare result with :cite:`Driesner2007Part2`.
"""
import os
import numpy as np
import time
import linecache
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import patches
import matplotlib.ticker as ticker
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator
from tabulate import tabulate
import copy
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

def plot_3d(sw,scale='linear'):
    fig=plt.figure(figsize=(14,14))
    ax = fig.add_subplot(111,projection='3d',facecolor='None')
    xcenter=1 # X<1% wt.% use log scale, if X>=1% wt.% NaCl, use linear scale
    axtrans=[]
    if(scale=='loglinear'):
        axtrans=helpfunc.set_axis_diagram_3D_loglinearx(ax,xcenter=xcenter,ratio_log_lin=(1,1),zlim=(sw.pmin()/1E5,600),xMajor_loc_log=1,xMajor_loc_linear=10,xMinor_loc_linear=2,zMajor_loc=100,zMinor_loc=20,yMajor_loc=200,xlim=(1E-10,100),xlabel='log$_{\mathregular{10}}$(Wt.% NaCl)')
    else:
        helpfunc.set_axis_diagram_3D(ax,zlim=(sw.pmin()/1E5,600),zMajor_loc=100,zMinor_loc=20)
    transform_X = lambda X : helpfunc.data2axis_loglin(axtrans,X) if(scale=='loglinear') else X
    # calculate boundary surface of VLH
    T = np.append(np.linspace(H2ONaCl.T_MIN_VLH, sw.Tmax_VLH()-10, 100), np.linspace(sw.Tmax_VLH()-10, sw.Tmax_VLH(), 50))
    P_vlh = np.array(sw.P_VLH(T))
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
    ax.plot_wireframe(transform_X(XX*100),TT-273.15,PP/1E5, color='darkblue',lw=1,label='V+L+H: vapor->liquid')
    # ax.plot_surface(transform_X(XX*100),TT-273.15,PP/1E5,ec='gray',linewidth=0.1)
    # 2. liquid -> halite region
    for j in range(0,PP.shape[1]):
        XX[:,j] = np.linspace(Xl_vlh[j], 1,len(X))
    ax.plot_wireframe(transform_X(XX*100),TT-273.15,PP/1E5, color='orange',lw=1,label='V+L+H: liquid->halite')
    # ax.plot_surface(transform_X((XX*100)),TT-273.15,PP/1E5,ec='gray',linewidth=0.1)
    # ax.plot(X_haliteLiquidus*100, T-273.15, P/1E5,color='orange',label='V+L+H: liquid',lw=1,zorder=10)
    # VLH: halite
    ax.plot(transform_X(Xl_vlh*0+100), T-273.15, P_vlh/1E5, color='k',label='VLH: halite',lw=2)
    # VLH: liquid
    ax.plot(transform_X(Xl_vlh*100), T-273.15, P_vlh/1E5, color='g',label='VLH: liquid',lw=2)
    # VLH: vapor
    ax.plot(transform_X(Xv_vlh*100), T-273.15, P_vlh/1E5, color='r',label='VLH: vapor',lw=2,zorder=11)
    # project halite saturated vapor on X-T plane
    ax.plot(transform_X(Xv_vlh*100), T-273.15, T*0+sw.pmin()/1E5, color='r',label='Fig.11 of Driesner & Heinrich(2007)',ls=(0,(2,1)),lw=2)
    leg=ax.legend()
    # change legend handle of wireframe of phase boundaries to wireframe hatch
    for i in [0,1]:
        leg.legendHandles[i]=Patch(facecolor='white', edgecolor=leg.legendHandles[i]._color,linewidth=0.0,label=leg.texts[i]._text,hatch='++++')
    ax.legend(handles=leg.legendHandles, loc='upper left',ncol=7)
    # text
    helpfunc.text3d(ax,(transform_X(1E-8),800,sw.pmin()/1E5),"P$_{min}$=%.0f bar"%(sw.pmin()/1E5),size=0.07,angle=-90,ec='None',fc='k')
    savefig('PhaseBoundary_VLH_3D_%s'%(scale))

# %%
# General 3D view
# --------------------------------
plot_3d(sw_84,scale='linear')
plot_3d(sw_84,scale='loglinear')

# %%
# Benchmark comparison
# --------------------------
# Compare result of |xThermo| and :cite:`Driesner2007Part2`, and also compare result based on different EOS of water.
def contourf_phase(ax,TT,pp,phase,phase_name,ax_cb=None):
    cmap = plt.get_cmap("Dark2")
    # customize cmap
    colors=list(copy.deepcopy(cmap.colors))
    colors[0:8]=['lightblue','red','lightgreen','lightgray','violet','yellow','lightcyan','lightcyan']
    cmap.colors=tuple(colors)
    CS=ax.contourf(TT,pp,phase, cmap=cmap,vmin=phase.min()-0.5, vmax=phase.max()+0.5, levels=np.linspace(phase.min()-0.5,phase.max()+0.5,len(phase_name)+1))
    if(ax_cb is None): ax_cb = ax.inset_axes([0,1.03,1,0.05])
    cb=plt.colorbar(CS, cax=ax_cb, orientation='horizontal',ticklocation='top',ticks=np.arange(phase.min(),phase.max()+1))
    cb.ax.set_xticklabels(phase_name)
    return CS,ax_cb,cb
def plot_props_VLH(sw,water,mmc5='../Driesner2007a/1-s2.0-S0016703707002943-mmc5.txt',mmc3='../Driesner2007b/1-s2.0-S0016703707002955-mmc3.txt'):
    if(not os.path.exists(mmc5)):
        print('Please set correct mmc1 file path: %s'%(mmc5))
        exit()
    data=np.loadtxt(mmc5, skiprows=7)
    T0_5,P0_5,XV0_5,XL0_5=data[:,0]+273.15,data[:,1]*1E5,data[:,2],data[:,3]
    if(not os.path.exists(mmc3)):
        print('Please set correct mmc1 file path: %s'%(mmc3))
        exit()
    data=np.loadtxt(mmc3, skiprows=5)
    T0,P0,XV0,rhoV0,hV0,XL0,rhoL0,hL0=data[:,0]+273.15,data[:,1]*1E5,data[:,2],data[:,3],data[:,4],data[:,5],data[:,6],data[:,7]
    # calculate
    T = np.linspace(H2ONaCl.T_MIN_VLH, H2ONaCl.T_MAX_VLH, 500)
    P = np.array(sw.P_VLH(T))
    # X_haliteLiquidus=np.array(sw.X_HaliteLiquidus(T,P))
    XL_,XV_ = np.array(sw.X_VLH(T,P))
    rhoV_, rhoL_ = np.array(sw.Rho_phase(T, P, XV_, H2ONaCl.Vapor)), np.array(sw.Rho_phase(T, P, XL_, H2ONaCl.Liquid))
    hV_, hL_     = np.array(sw.H_phase(T, P, XV_, H2ONaCl.Vapor)), np.array(sw.H_phase(T, P, XL_, H2ONaCl.Liquid))

    # plot
    fig,axes = plt.subplots(1,6,figsize=(30,4),gridspec_kw={'wspace':0.1},sharey=False)
    # 1. liquid salinity
    ax=axes[0]
    line=helpfunc.plot_coloredline(ax,np.array(sw.Wt2Mol(XL_)), T-273.15, P/1E5,cmap='rainbow')
    ax_cb = ax.inset_axes([0.8,0.15,0.02,0.6])
    fig.colorbar(line,cax = ax_cb,label='Pressure (bar)')
    ax.plot(XL0_5, T0_5-273.15,color='gray',ls='dashed', marker='.',markersize=5,mec='w',mfc='k',markeredgewidth=0.3,label='Driesner & Heinrich(2007)')
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(MultipleLocator(0.04))
    ax.set_xlabel('Liquid composition X$_{\mathregular{NaCl}}$ (mole fraction)')
    # 2. vapor salinity
    ax=axes[1]
    helpfunc.plot_coloredline(ax,np.array(sw.Wt2Mol(XV_)), T-273.15, P/1E5,cmap='rainbow')
    ax.set_xscale('log')
    ax.plot(XV0_5, T0_5-273.15,color='gray',ls='dashed', marker='.',markersize=5,mec='w',mfc='k',markeredgewidth=0.3,label='Driesner & Heinrich(2007)')
    ax.set_xlabel('Vapor composition X$_{\mathregular{NaCl}}$ (mole fraction)')
    # 3. liquid density
    ax=axes[2]
    helpfunc.plot_coloredline(ax,rhoL_, T-273.15, P/1E5,cmap='rainbow')
    ax.plot(rhoL0, T0-273.15,color='gray',ls='dashed', marker='.',markersize=5,mec='w',mfc='k',markeredgewidth=0.3,label='Driesner(2007)')
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.xaxis.set_minor_locator(MultipleLocator(20))
    ax.set_xlabel('Liquid density (kg/m$^{\mathregular{3}}$)')
    # 4. Vapor density
    ax=axes[3]
    helpfunc.plot_coloredline(ax,rhoV_, T-273.15, P/1E5,cmap='rainbow')
    ax.plot(rhoV0, T0-273.15,color='gray',ls='dashed', marker='.',markersize=5,mec='w',mfc='k',markeredgewidth=0.3,label='Driesner(2007)')
    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.xaxis.set_minor_locator(MultipleLocator(4))
    ax.set_xlabel('Vapor density (kg/m$^{\mathregular{3}}$)')
    # 5. liquid enthalpy
    ax=axes[4]
    helpfunc.plot_coloredline(ax,hL_/1E6, T-273.15, P/1E5,cmap='rainbow')
    ax.plot(hL_/1E6, T-273.15,'.')
    ax.plot(hL0/1E6, T0-273.15,color='gray',ls='dashed', marker='.',markersize=5,mec='w',mfc='k',markeredgewidth=0.3,label='Driesner(2007)')
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.set_xlabel('Liquid enthalpy (MJ/kg)')
    # 6. vapor enthalpy
    ax=axes[5]
    helpfunc.plot_coloredline(ax,hV_/1E6, T-273.15, P/1E5,cmap='rainbow')
    ax.plot(hV0/1E6, T0-273.15,color='gray',ls='dashed', marker='.',markersize=5,mec='w',mfc='k',markeredgewidth=0.3,label='Driesner(2007)')
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.set_xlabel('Vapor enthalpy (MJ/kg)')

    for ax in axes:
        ax.yaxis.set_major_locator(MultipleLocator(100))
        ax.yaxis.set_minor_locator(MultipleLocator(20))
        ax.grid(which='major',lw=0.04,color='k')
        ax.grid(which='minor',lw=0.04,color='gray')
        ax.axhline(H2ONaCl.T_MIN_VLH-273.15,label='T=%.2f $^{\circ}$C, p=1 bar'%(H2ONaCl.T_MIN_VLH-273.15),color='r',ls='dotted')
        ax.legend()
        # ax.set_ylim(300,T.max()-273.15 + 20)
    axes[0].set_ylabel('Temperature ($^{\circ}$C)')
    savefig('VLH_props_%s'%(sw.name_backend()))

    # let's checkout what happens in the (low T, low p) and (high T, low p) region for the liquid enthalpy
    q1,q2 = np.array(sw.q1q2_Tstar_H(P, XL_))
    q1_v,q2_v = np.array(sw.q1q2_Tstar_H(P, XV_))
    Tstar_H,Tstar_H_v = q1 + q2*(T-273.15) + 273.15, q1_v + q2_v*(T-273.15) + 273.15
    n1,n2 = np.array(sw.n1n2_Tstar_V(P, XL_))
    n1_v,n2_v = np.array(sw.n1n2_Tstar_V(P, XV_))
    Tstar_V,Tstar_V_v = n1 + n2*(T-273.15) + 273.15, n1_v + n2_v*(T-273.15) + 273.15
    T_water,p_water=np.linspace(np.array([T.min(),Tstar_H.min(),Tstar_H_v.min(),Tstar_V.min(),Tstar_V_v.min()]).min(),np.array([T.max(),Tstar_H.max(),Tstar_H_v.max(),Tstar_V.max(),Tstar_V_v.max()]).max()+100,1000), np.linspace(np.log10(P.min()),np.log10(P.max())+0.1,500)
    TT,pp=np.meshgrid(T_water,10**p_water)
    phase,rho,h = np.zeros_like(TT),np.zeros_like(TT),np.zeros_like(TT)
    for i in range(0,TT.shape[0]):
        for j in range(0,TT.shape[1]):
            props = water.UpdateState_TPX(TT[i][j], pp[i][j])
            phase[i][j]=props.phase
            rho[i][j]=props.Rho
            h[i][j]=props.H
    # get phase names
    phase_unique = np.sort(np.unique(phase))
    phase_name = ['']*len(phase_unique)
    for i,phase0 in enumerate(phase_unique):
        phase[phase==phase0]=i+phase_unique.max()+10
        phase_name[i]=water.phase_name(int(phase0))
    fig,axes=plt.subplots(1,3,figsize=(21,5),gridspec_kw={'wspace':0.05},sharey=True)
    axes[0].set_ylabel('Pressure (bar)')
    l_L,l_V=[],[]
    for ax,prop, cmap,label in zip(axes, [phase, rho, h/1E6],['Paired','YlGnBu_r','RdBu'],['Phase','Density (kg/m$^{\mathregular{3}}$)','Specific enthalpy (MJ/kg)']):
        CS=[]
        if(ax==axes[0]):
            CS,ax_cb,cb=contourf_phase(ax,TT-273.15,pp/1E5,prop,phase_name,ax.inset_axes([0,1.03,1,0.03]))
        else:
            CS=ax.contourf(TT-273.15,pp/1E5,prop,levels=50,cmap=cmap)
            ax_cb=ax.inset_axes([0,1.02,1,0.03])
            plt.colorbar(CS,cax=ax_cb,label=label,orientation='horizontal')
            ax_cb.xaxis.set_label_position('top')
            ax_cb.xaxis.set_ticks_position('top')
        l_L,=ax.plot(Tstar_H-273.15, P/1E5,lw=2,label='$T^*_h$: liquid')
        l_V,=ax.plot(Tstar_H_v-273.15, P/1E5,lw=2,label='$T^*_h$: vapor')
        ax.plot(Tstar_V-273.15, P/1E5,lw=2,ls='dashed',label='$T^*_V$: liquid')
        ax.plot(Tstar_V_v-273.15, P/1E5,lw=2,ls='dashed',label='$T^*_V$: vapor')
        ax.plot(T-273.15, P/1E5,lw=0.8,marker='.',markevery=20,mec='w',mew=0.5,ms=10,label='$T_{VLH}$')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlabel('Temperature ($^{\circ}$C)')
    axes[0].legend(loc='center left') #,bbox_to_anchor=[1.01,0]
    # x0,y0=Tstar_H[int(len(Tstar_H)/4)]-273.15,P[int(len(Tstar_H)/4)]/1E5
    # axes[0].annotate("$T^{*}_h - p$ path of\nsaturated liquid phase\non VLH coexistence",
    #                     xy=(x0,y0),xytext=(x0-100,y0), ha='right',va='center',bbox={'fc':'None','ec':l_L.get_color()},fontsize=14,fontweight='bold',
    #                     arrowprops=dict(arrowstyle="->",connectionstyle="arc3"),)
    # x0,y0=Tstar_H_v[int(len(Tstar_H)/5)]-273.15,P[int(len(Tstar_H)/5)]/1E5
    # axes[0].annotate("$T^{*}_h - p$ path of\nsaturated vapor phase\non VLH coexistence",
    #                  xy=(x0,y0),xytext=(x0,y0/100), ha='center',va='center',bbox={'fc':'None','ec':l_V.get_color()},fontsize=14,fontweight='bold',
    #                  arrowprops=dict(arrowstyle="->",connectionstyle="arc3"),)
    savefig('Tstar_VLH_%s'%(sw.name_backend()))
def mscatter(ax,x,y, m=None, **kw):
    import matplotlib.markers as mmarkers
    if not ax: ax=plt.gca()
    sc = ax.scatter(x,y,**kw)
    if (m is not None) and (len(m)==len(x)):
        paths = []
        for marker in m:
            if isinstance(marker, mmarkers.MarkerStyle):
                marker_obj = marker
            else:
                marker_obj = mmarkers.MarkerStyle(marker)
            path = marker_obj.get_path().transformed(
                marker_obj.get_transform())
            paths.append(path)
        sc.set_paths(paths)
    return sc
def plot_err(ax,x,y,data,label='',cmap='rainbow',markers_VL=None,scale_data='log',vmin=1E-6,vmax=1):
    # plot difference between xThermo and Driesner(2007b)
    norm = None
    if(scale_data=='log'):
        norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = mpl.colors.CenteredNorm(vcenter=1)
        cmap = 'seismic'
    CS=mscatter(ax, x, y,c=data,m=markers_VL,cmap=cmap,norm=norm,zorder=3)
    ax_cb=ax.inset_axes([0,1.02,1,0.05])
    plt.colorbar(CS,cax=ax_cb,label=label,orientation='horizontal',extend='both')
    ax_cb.xaxis.set_label_position('top')
    ax_cb.xaxis.set_ticks_position('top')
    if(markers_VL is not None):
        ax.set_xlim(ax.get_xlim())
        ax.set_ylim(ax.get_ylim())
        ax.plot(-1,-1,markers_VL[0],label='Vapor')
        ax.plot(-1,-1,markers_VL[-1],label='Liquid')
        ax.legend(ncol=2)
def benchmark_VLH(sw,mmc3='../Driesner2007b/1-s2.0-S0016703707002955-mmc3.txt'):
    # compare
    if(not os.path.exists(mmc3)):
        print('Please set correct mmc1 file path: %s'%(mmc3))
        exit()
    data=np.loadtxt(mmc3, skiprows=5)
    T0,P0,XV0,rhoV0,hV0,XL0,rhoL0,hL0=data[:,0]+273.15,data[:,1]*1E5,data[:,2],data[:,3],data[:,4],data[:,5],data[:,6],data[:,7]
    ind=(T0>sw.Tmin_VLH())
    # only compare the result in valid range of pressure: >1bar
    T0,P0,XV0,rhoV0,hV0,XL0,rhoL0,hL0 = T0[ind],P0[ind],XV0[ind],rhoV0[ind],hV0[ind],XL0[ind],rhoL0[ind],hL0[ind]
    # 1. calculate halite liquidus
    P_=np.array(sw.P_VLH(T0))
    XL_,XV_ = np.array(sw.X_VLH(T0,P_))
    XL_mol_,XV_mol_ = np.array(sw.Wt2Mol(XL_)), np.array(sw.Wt2Mol(XV_))
    # 2. calculate saturated liquid density and vapor density
    rhoV_, rhoL_ = np.array(sw.Rho_phase(T0, P_, XV_, H2ONaCl.Vapor)), np.array(sw.Rho_phase(T0, P_, XL_, H2ONaCl.Liquid))
    hV_, hL_     = np.array(sw.H_phase(T0, P_, XV_, H2ONaCl.Vapor)), np.array(sw.H_phase(T0, P_, XL_, H2ONaCl.Liquid))
    # compare result dict
    Data0 = {'p':P0/1E5,'XV':XV0,'rhoV':rhoV0,'hV':hV0,'XL':XL0,'rhoL':rhoL0,'hL':hL0}
    Data_ = {'p':P_/1E5,'XV':XV_mol_,'rhoV':rhoV_,'hV':hV_,'XL':XL_mol_,'rhoL':rhoL_,'hL':hL_}
    Err,RErr={},{}
    for key in Data0.keys(): Err[key],RErr[key] = Data0[key]-Data_[key], np.abs(Data0[key]-Data_[key])/(Data0[key])*100.0
    # print to file
    fpout = open('%s/mmc3_%s.csv'%(result_path,sw.name_backend()),'w')
    fpout.write('T[C],P(Driesner)[bar],P(xThermo),P(diff)[bar],XV(Driesner)[mol],XV(xThermo)[mol],XV(diff)[mol],RhoV(Driesner)[kg/m3],RhoV(xThermo),RhoV(err),HV(Driesner)[J/kg],HV(xThermo),HV(err),XL(Driesner)[mol],XL(xThermo)[mol],XL(diff)[mol],RhoL(Driesner)[kg/m3],RhoL(xThermo),RhoL(err),HL(Driesner)[J/kg],HL(xThermo),HL(err)\n')
    for i in range(0,len(T0)):
        fpout.write('%.6e'%(T0[i]-273.15))
        for key in Data0.keys():
            fpout.write(',%.6e,%.6e,%.6e'%(Data0[key][i], Data_[key][i],compare(Data0[key][i],Data_[key][i])))
        fpout.write('\n')
    fpout.close()
    # plot difference
    offset_p = -100
    fig,axes=plt.subplots(1,4,figsize=(21,5),sharey=True,gridspec_kw={'wspace':0.1})
    plot_err(axes[0], T0, P0/1E5, RErr['p'],'Relative difference: Pressure (%)')
    markers=np.repeat(["o", "*"], len(T0))
    plot_err(axes[1], np.append(T0,T0), np.append(P0/1E5,P0/1E5+offset_p), np.append(RErr['XV'],RErr['XL']),'Relative difference: Composition (%)',markers_VL=markers)
    plot_err(axes[2], np.append(T0,T0), np.append(P0/1E5,P0/1E5+offset_p), np.append(RErr['rhoV'],RErr['rhoL']),'Relative difference: Density (%)',markers_VL=markers)
    plot_err(axes[3], np.append(T0,T0), np.append(P0/1E5,P0/1E5+offset_p), np.append(RErr['hV'],RErr['hL']),'Relative difference: Specific enthalpy (%)',markers_VL=markers)
    for ax in axes:
        ax.yaxis.set_major_locator(MultipleLocator(100))
        ax.yaxis.set_minor_locator(MultipleLocator(20))
        ax.xaxis.set_major_locator(MultipleLocator(100))
        ax.xaxis.set_minor_locator(MultipleLocator(20))
        ax.grid(which='major',lw=0.04,color='k')
        ax.grid(which='minor',lw=0.04,color='gray')
        ax.set_xlabel('Temperature ($^{\circ}$C)')
    axes[0].set_ylabel('Pressure (bar)')
    savefig('diff_VLH')
    # statistics of the difference
    table=[]
    for key,name in zip(list(Err.keys()),['Pressure (bar)','XV (mole fraction)','RhoV (kg/m3)','HV (J/kg)','XL (mole fraction)','RhoL (kg/m3)','HL (J/kg)']):
        RErr[key] = RErr[key][~(np.isnan(RErr[key]) | np.isinf(RErr[key]))]
        table.append([name,Err[key].min(),Err[key].max(),RErr[key].min(),RErr[key].max()])
    print(tabulate(table, headers=['Critical property', 'Err. Min', 'Err. Max','RE. Min(%)','RE. Max(%)'],numalign="right",floatfmt=".6f"))

# %%
# Based on IAPS84 EOS
# ^^^^^^^^^^^^^^^^^^^^^^^
plot_props_VLH(sw_84,H2O.cIAPS84())
benchmark_VLH(sw_84)

# %%
# Based on IAPWS95 EOS
# ^^^^^^^^^^^^^^^^^^^^^^^
# plot_props_VLH(sw_95,H2O.cIAPWS95_CoolProp())
# benchmark_VLH(sw_95)

# %%
# .. warning::
#
#     The :math:`T^*_h` equation (22) in :cite:`Driesner2007Part2` seems not valid for saturated liquid phase on VLH surface when in
#     (low T, low p) and (high T, low p) regions, because the corrected :math:`T^*_h-p` path cross the boiling curve of water (see above),
#     therefore there are some value jupms on liquid enthalpy curve.
#     However, :cite:`Driesner2007Part2` gives some smooth values in low-p regions (see dashed lines with dots in above figures), how to get such values?
#     Well, for the hydrothermal modeling, this issue is not that critical because the pressure always high (>=50bar).

# %%
# Result table
# -----------------------------
#
# Result data calculated by |xThermo| based on both water EOS of IAPS84 and IAPWS95.
#
# .. seealso::
#
#     |mmc3| in :cite:`Driesner2007Part2` and Fig. 11 in :cite:`Driesner2007Part1`.
#
# .. tab:: IAPS84
#
#     .. csv-table:: Comparison between result of :cite:`Driesner2007Part2` and result calculated by |xThermo|.
#         :file: mmc3_IAPS84.csv
#         :header-rows: 1
#
# .. tab:: IAPWS95
#
#     .. csv-table:: Comparison between result of :cite:`Driesner2007Part2` and result calculated by |xThermo|.
#         :file: mmc3_IAPWS95.csv
#         :header-rows: 1
#
#
# .. tip::
#
#     The help function for 3D plot can be downloaded at here: :download:`helpfunc.py`