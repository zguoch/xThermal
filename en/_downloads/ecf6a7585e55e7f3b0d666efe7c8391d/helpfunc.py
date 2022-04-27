# -*- coding: utf-8 -*-
"""
Help function definition
===========================
.. include:: /include.rst_
"""
import numpy as np
import time
import linecache
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import patches
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
mpl.rcParams['font.family'] = 'Arial'  #default font family
mpl.rcParams['mathtext.fontset'] = 'cm' #font for math
import os
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

import copy
from matplotlib.path import Path
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
from matplotlib.transforms import Affine2D
from mpl_toolkits.mplot3d import Axes3D, art3d
from matplotlib.patches import PathPatch

def text3d(ax, xyz, s, zdir="z", size=0.1, angle=0,font='Arial',weight='normal',ha='left',va='center', **kwargs):
    x, y, z = xyz
    xlim,ylim,zlim=ax.get_xlim(),ax.get_ylim(),ax.get_zlim()
    xmin,xmax,ymin,ymax,zmin,zmax=xlim[0],xlim[1],ylim[0],ylim[1],zlim[0],zlim[1]
    xlen,ylen,zlen=xmax-xmin, ymax-ymin, zmax-zmin
    minLen_axis=np.min([xlen,ylen,zlen])
    aspect=ax.get_box_aspect()
    fontscale=1
    if   zdir == "y": xy, z, fontscale = (x,z), y, (xlen/zlen)/(aspect[0]/aspect[2])
    elif zdir == "x": xy, z, fontscale = (y,z), x, (ylen/zlen)/(aspect[1]/aspect[2])
    else:             xy, z, fontscale = (x,y), z, (xlen/ylen)/(aspect[0]/aspect[1])
    path = TextPath((0, 0), s, size=size, prop = FontProperties(family=font,weight=weight))
    V = copy.deepcopy(path.vertices)
    if(ha=='center'):
        V[:,0] -= (V[:,0].max() - V[:,0].min())/2 #
    elif(ha=='right'):
        V[:,0] -= V[:,0].max()
    if(va=='center'):
        V[:,1] -= (V[:,1].max() - V[:,1].min())/2 # 居中
    elif(va=='top'):
        V[:,1] -= V[:,1].max()
    trans = Affine2D().rotate(angle/180*np.pi).scale(1,1/fontscale).translate(xy[0], xy[1])
    # path = PathPatch(trans.transform_path(path), clip_on=False, **kwargs)
    path = PathPatch(trans.transform_path(Path(V, path._codes)), clip_on=False, **kwargs)
    ax.add_patch(path)
    art3d.pathpatch_2d_to_3d(path, z=z, zdir=zdir)

def niceAxis_3D(ax, lw_major=0.5,lw_minor=0.25,color_major="0.50",color_minor="0.75",
             alpha_pane=0.2,fill_pane=True,ec_pane='None',label3D=True,fs_label=0.1,
             length_major=0.05,length_minor=0.02,offset_axislabel=3,frame_on=True,scaled=True):
    # 坐标轴自定义
    xlim,ylim,zlim=ax.get_xlim(),ax.get_ylim(),ax.get_zlim()
    xmin,xmax,ymin,ymax,zmin,zmax=xlim[0],xlim[1],ylim[0],ylim[1],zlim[0],zlim[1]
    # 1. 根据坐标轴范围设置坐标轴比例：假设xyz三轴以自然比例显示
    xlen,ylen,zlen=xmax-xmin, ymax-ymin, zmax-zmin
    minLen_axis=np.min([xlen,ylen,zlen])
    xratio,yratio,zratio=[xlen,ylen,zlen]/minLen_axis
    if(scaled):
        ax.set_box_aspect((4*xratio,4*yratio,3*zratio))  # 默认比例是4：4：3
    # 3. 重新定义网格线：注意mpl的三维绘图目前不支持grid()函数分别定义主刻度网格和副刻度网格属性
    ax.grid(False) # 首先关闭默认网格
    zorder_grid=-50
    try:
        ax.get_figure().draw_without_rendering() # new in 3.5: https://matplotlib.org/stable/users/prev_whats_new/whats_new_3.5.0.html#figure-now-has-draw-without-rendering-method
    except:
        try:
            ax.get_figure().canvas.draw() # Old version
        except:
            pass
    if(ax.get_xticklabels()[0]._text==''): return
    # 主刻度线
    linewidth, color, length = lw_major, color_major, length_major
    if(type(length_major)==type([1,1,1])):
        length=length_major[0]
    # x轴
    tickmin,tickmax,tickaxis,axislim=xmin,xmax,ax.xaxis,ax.get_xlim()
    tickpos1,tickpos2=[ymin,ymax,ymax],[zmin,zmin,zmax]
    for tick,label in zip(ax.get_xticks(),ax.get_xticklabels()):
        if((tick<axislim[0]) | (tick>axislim[1])):
            continue
        ax.plot([tick,tick,tick],tickpos1,tickpos2, lw=linewidth,color=color,zorder=zorder_grid)
        if(label3D):
            text3d(ax, (tick, ymin-length*1.2*ylen, zmin), label.get_text(),size=fs_label*minLen_axis, fc=label.get_color(), ec="None", ha='center',va='top')
            label.set_alpha(0)
            ax.plot([tick,tick],[ymin,ymin-length*ylen],[zmin,zmin],lw=linewidth,color=label.get_color())
    if(label3D):
        ax.plot([xmin,xmax],[ymin,ymin],[zmin,zmin],lw=linewidth,color=label.get_color())
        text3d(ax, ((xmin+xmax)/2, ymin-length*offset_axislabel*ylen, zmin), ax.get_xlabel(),
               size=fs_label*minLen_axis, fc=label.get_color(), ec="None",ha='center',va='top')
    # y轴
    tickmin,tickmax,tickaxis,axislim=ymin,ymax,ax.yaxis,ax.get_ylim()
    tickpos1,tickpos2=[xmax,xmax,xmin],[zmax,zmin,zmin]
    if(type(length_major)==type([1,1,1])):
        length=length_major[1]
    for tick,label in zip(ax.get_yticks(),ax.get_yticklabels()):
        if((tick<axislim[0]) | (tick>axislim[1])):
            continue
        ax.plot(tickpos1,[tick,tick,tick],tickpos2,
                lw=linewidth,color=color,zorder=zorder_grid)
        if(label3D):
            text3d(ax, (xmin-length*1.1*xlen, tick, zmin), label.get_text(),size=fs_label*minLen_axis,
                   fc=label.get_color(), ec="None",ha='center',va='top',angle=-90)
            label.set_alpha(0)
            ax.plot([xmin,xmin-length*xlen],[tick,tick],[zmin,zmin],lw=linewidth,color=label.get_color())
    if(label3D):
        ax.plot([xmax,xmax],[ymin,ymax],[zmin,zmin],lw=linewidth,color=label.get_color(),ls='dashed')
        text3d(ax, (xmin-length*offset_axislabel*xlen,(ymin+ymax)/2, zmin), ax.get_ylabel(),
               size=fs_label*minLen_axis, fc=label.get_color(), ec="None",ha='center',va='top',angle=-90)

    # z轴
    tickmin,tickmax,tickaxis,axislim=zmin,zmax,ax.zaxis,ax.get_zlim()
    tickpos1,tickpos2=[xmax,xmax,xmin],[ymin,ymax,ymax]
    if(type(length_major)==type([1,1,1])):
        length=length_major[2]
    for tick,label in zip(ax.get_zticks(),ax.get_zticklabels()):
        if((tick<axislim[0]) | (tick>axislim[1])):
            continue
        ax.plot(tickpos1,tickpos2,[tick,tick,tick],
                lw=linewidth,color=color,zorder=zorder_grid)
        if(label3D):
            text3d(ax, (xmin-1.1*length*xlen,ymax,tick), label.get_text(),zdir='y',size=fs_label*minLen_axis,
                   fc=label.get_color(), ec="None",ha='center',va='bottom',angle=90)
            label.set_alpha(0)
            ax.plot([xmin,xmin-length*minLen_axis],[ymax,ymax],[tick,tick],lw=linewidth,color=label.get_color(),clip_on=False)
    if(label3D):
        ax.plot([xmin,xmin],[ymax,ymax],[zmin,zmax],lw=linewidth,color=label.get_color())
        text3d(ax, (xmin-length*offset_axislabel*xlen,ymax, (zmin+zmax)/2), ax.get_zlabel(),zdir='y',
               size=fs_label*minLen_axis, fc=label.get_color(), ec="None",ha='center',va='bottom',angle=90)
    if(frame_on==True):
        ax.plot([xmax,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymax],[zmin,zmax,zmax,zmax,zmin],color=label.get_color(),lw=linewidth)
        ax.plot([xmin,xmin,xmin],[ymin,ymax,ymax],[zmin,zmin,zmax],color=label.get_color(),lw=linewidth)
        ax.plot([xmin,xmax],[ymax,ymax],[zmin,zmin],color=label.get_color(),lw=linewidth,ls='dashed')
    # 副刻度线
    linewidth, color, length = lw_minor, color_minor,length_minor
    # x轴
    tickmin,tickmax,tickaxis,axislim=xmin,xmax,ax.xaxis,ax.get_xlim()
    tickpos1,tickpos2=[ymin,ymax,ymax],[zmin,zmin,zmax]
    for tick in tickaxis.get_minor_locator().tick_values(tickmin,tickmax):
        if((tick<axislim[0]) | (tick>axislim[1])):
            continue
        ax.plot([tick,tick,tick],tickpos1,tickpos2,
                lw=linewidth,color=color,zorder=zorder_grid)
        if(label3D):
            ax.plot([tick,tick],[ymin,ymin-length*ylen],[zmin,zmin],lw=linewidth,color=label.get_color())
    # y轴
    tickmin,tickmax,tickaxis,axislim=ymin,ymax,ax.yaxis,ax.get_ylim()
    tickpos1,tickpos2=[xmax,xmax,xmin],[zmax,zmin,zmin]
    for tick in tickaxis.get_minor_locator().tick_values(tickmin,tickmax):
        if((tick<axislim[0]) | (tick>axislim[1])):
            continue
        ax.plot(tickpos1,[tick,tick,tick],tickpos2,
                lw=linewidth,color=color,zorder=zorder_grid)
        if(label3D):
            ax.plot([xmin,xmin-length*xlen],[tick,tick],[zmin,zmin],lw=linewidth,color=label.get_color())
    # z轴
    tickmin,tickmax,tickaxis,axislim=zmin,zmax,ax.zaxis,ax.get_zlim()
    tickpos1,tickpos2=[xmax,xmax,xmin],[ymin,ymax,ymax]
    for tick in tickaxis.get_minor_locator().tick_values(tickmin,tickmax):
        if((tick<axislim[0]) | (tick>axislim[1])):
            continue
        ax.plot(tickpos1,tickpos2,[tick,tick,tick],
                lw=linewidth,color=color,zorder=zorder_grid)
        if(label3D):
            ax.plot([xmin,xmin-length*xlen],[ymax,ymax],[tick,tick],lw=linewidth,color=label.get_color(),clip_on=False)
    # 3. 设置ticks和labels以及spines和panes
    for axis,inward_factor in zip([ax.xaxis,ax.yaxis,ax.zaxis],[0.3,0.2,0.2]):
        axis._axinfo['tick']['inward_factor'] = inward_factor
        axis._axinfo['tick']['outward_factor'] = 0
    ax.tick_params(which='minor', color=(0,0,1,0) ) # 取消副刻度
    # xaxis
    color,alpha,pad,axis,waxis=(1,0,0),alpha_pane,-3,ax.xaxis,ax.w_xaxis
    ax.tick_params(axis='x', which='major', pad=pad)
    axis.set_pane_color(color+(alpha,))
    axis.pane.set_edgecolor(ec_pane)
    axis.pane.fill = fill_pane
    waxis.line.set_color(ax.get_xticklines()[0].get_color())
    # yaxis
    color,alpha,pad,axis,waxis=(0,1,0),alpha_pane,-3,ax.yaxis,ax.w_yaxis
    ax.tick_params(axis='y', which='major', pad=pad)
    axis.set_pane_color(color+(alpha,))
    axis.pane.set_edgecolor(ec_pane)
    axis.pane.fill = fill_pane
    waxis.line.set_color(ax.get_yticklines()[0].get_color())
    # zaxis
    color,alpha,pad,axis,waxis=(0,0,1),alpha_pane,0,ax.zaxis,ax.w_zaxis
    ax.tick_params(axis='z', which='major', pad=pad)
    axis.set_pane_color(color+(alpha,))
    axis.pane.set_edgecolor(ec_pane)
    axis.pane.fill = fill_pane
    waxis.line.set_color(ax.get_zticklines()[0].get_color())
    if(label3D):
        ax.axis('off')

def niceAxis_3D_loglinearx(ax, data_lim,xcenter, lw_major=0.5,lw_minor=0.25,color_major="0.50",color_minor="0.75",
             alpha_pane=0.2,fill_pane=True,ec_pane='None',label3D=True,fs_label=0.1,
             length_major=0.05,length_minor=0.02,offset_axislabel=3,frame_on=True,scaled=True,xMajor_loc_log=2,xMinor_loc_log=0.4,xMajor_loc_linear=10,xMinor_loc_linear=2,process_x=False):
    color_log,color_linear = 'k','k'
    # 坐标轴自定义
    xlim,ylim,zlim=ax.get_xlim(),ax.get_ylim(),ax.get_zlim()
    xmin,xmax,ymin,ymax,zmin,zmax=xlim[0],xlim[1],ylim[0],ylim[1],zlim[0],zlim[1]
    # 1. 根据坐标轴范围设置坐标轴比例：假设xyz三轴以自然比例显示
    xlen,ylen,zlen=xmax-xmin, ymax-ymin, zmax-zmin
    minLen_axis=np.min([xlen,ylen,zlen])
    xratio,yratio,zratio=[xlen,ylen,zlen]/minLen_axis
    if(scaled):
        ax.set_box_aspect((8*xratio,4*yratio,3*zratio))  # 默认比例是4：4：3
    # 3. 重新定义网格线：注意mpl的三维绘图目前不支持grid()函数分别定义主刻度网格和副刻度网格属性
    ax.grid(False) # 首先关闭默认网格
    zorder_grid=-50
    try:
        ax.get_figure().draw_without_rendering() # new in 3.5: https://matplotlib.org/stable/users/prev_whats_new/whats_new_3.5.0.html#figure-now-has-draw-without-rendering-method
    except:
        try:
            ax.get_figure().canvas.draw() # Old version
        except:
            pass
    if(ax.get_xticklabels()[0]._text==''): return
    # 主刻度线
    linewidth, color, length = lw_major, color_major, length_major
    if(type(length_major)==type([1,1,1])):
        length=length_major[0]
    # x轴
    if(process_x):
        tickmin,tickmax,tickaxis,axislim=xmin,xmax,ax.xaxis,ax.get_xlim()
        # log段
        xmin_data,xmax_data = np.log10(data_lim[0]),np.log10(xcenter)
        xMajorLocator = mpl.ticker.MultipleLocator(xMajor_loc_log).tick_values(xmin_data,xmax_data)
        # 计算tick在坐标轴上的真正位置，将数据范围投影到坐标轴范围中
        xmin_axis,xmax_axis=ax.get_xlim()[0],0
        length_ = xmax_axis- xmin_axis
        scale_data2axes = length_/(xmax_data - xmin_data)
        tickpos1,tickpos2=[ymin,ymax,ymax],[zmin,zmin,zmax]
        for tick_data in xMajorLocator:
            if((tick_data<xmin_data) | (tick_data>xmax_data)):
                continue
            tick = (tick_data-xmin_data)*scale_data2axes + xmin_axis
            ticklabel = '10$^{\mathregular{%d}}$'%(tick_data)
            ax.plot([tick,tick,tick],tickpos1,tickpos2, lw=linewidth,color=color,zorder=zorder_grid)
            if(label3D):
                text3d(ax, (tick, ymin-length*1.2*ylen, zmin), ticklabel,size=fs_label*minLen_axis, fc=color_log, ec="None", ha='center',va='top')
                ax.plot([tick,tick],[ymin,ymin-length*ylen],[zmin,zmin],lw=linewidth,color=color_log)
        # linear和log的分界线
        # ax.plot([0,0],[ymin,ymax],[zmin,zmin],lw=linewidth,color='r')
        if(label3D):
            ax.plot([xmin_axis,xmax_axis],[ymin,ymin],[zmin,zmin],lw=linewidth,color=color_log)
            text3d(ax, ((xmin_axis+xmax_axis)/2, ymin-length*offset_axislabel*1.3*ylen, zmin), 'log scale',
                   size=fs_label*minLen_axis, fc=color_log, ec="None",ha='center',va='top')
        # linear段
        xmin_data,xmax_data = xcenter,data_lim[1]
        xMajorLocator = mpl.ticker.MultipleLocator(xMajor_loc_linear).tick_values(xmin_data,xmax_data)
        # 计算tick在坐标轴上的真正位置，将数据范围投影到坐标轴范围中
        xmin_axis,xmax_axis=0,ax.get_xlim()[1]
        length_ = xmax_axis- xmin_axis
        scale_data2axes = length_/(xmax_data - xmin_data)
        tickpos1,tickpos2=[ymin,ymax,ymax],[zmin,zmin,zmax]
        for tick_data in xMajorLocator:
            if((tick_data<xmin_data) | (tick_data>xmax_data)):
                continue
            tick = (tick_data-xmin_data)*scale_data2axes + xmin_axis
            ticklabel = '%.0f'%(tick_data)
            ax.plot([tick,tick,tick],tickpos1,tickpos2, lw=linewidth,color=color,zorder=zorder_grid)
            if(label3D):
                text3d(ax, (tick, ymin-length*1.2*ylen, zmin), ticklabel,size=fs_label*minLen_axis, fc=color_linear, ec="None", ha='center',va='top')
                ax.plot([tick,tick],[ymin,ymin-length*ylen],[zmin,zmin],lw=linewidth,color=color_linear)
        # linear和log的分界线
        ax.plot([0,0],[ymin,ymax],[zmin,zmin],lw=linewidth,color='k')
        if(label3D):
            ax.plot([xmin_axis,xmax_axis],[ymin,ymin],[zmin,zmin],lw=linewidth,color=color_linear)
            text3d(ax, ((xmin_axis+xmax_axis)/2, ymin-length*offset_axislabel*1.3*ylen, zmin), 'linear scale',
                   size=fs_label*minLen_axis, fc=color_linear, ec="None",ha='center',va='top')
            # 总的xlabel
            text3d(ax, ((ax.get_xlim()[0]+ax.get_xlim()[1])/2, ymin-length*offset_axislabel*1.5*ylen, zmin), ax.get_xlabel(),
                   size=fs_label*minLen_axis, fc='k', ec="None",ha='center',va='top')
    # y轴
    tickmin,tickmax,tickaxis,axislim=ymin,ymax,ax.yaxis,ax.get_ylim()
    tickpos1,tickpos2=[xmax,xmax,xmin],[zmax,zmin,zmin]
    if(type(length_major)==type([1,1,1])):
        length=length_major[1]
    for tick,label in zip(ax.get_yticks(),ax.get_yticklabels()):
        if((tick<axislim[0]) | (tick>axislim[1])):
            continue
        ax.plot(tickpos1,[tick,tick,tick],tickpos2,
                lw=linewidth,color=color,zorder=zorder_grid)
        if(label3D):
            text3d(ax, (xmin-length*1.1*xlen, tick, zmin), label.get_text(),size=fs_label*minLen_axis,
                   fc=label.get_color(), ec="None",ha='center',va='top',angle=-90)
            label.set_alpha(0)
            ax.plot([xmin,xmin-length*xlen],[tick,tick],[zmin,zmin],lw=linewidth,color=label.get_color())
    if(label3D):
        ax.plot([xmax,xmax],[ymin,ymax],[zmin,zmin],lw=linewidth,color=label.get_color(),ls='dashed')
        text3d(ax, (xmin-length*offset_axislabel*xlen,(ymin+ymax)/2, zmin), ax.get_ylabel(),
               size=fs_label*minLen_axis, fc=label.get_color(), ec="None",ha='center',va='top',angle=-90)

    # z轴
    tickmin,tickmax,tickaxis,axislim=zmin,zmax,ax.zaxis,ax.get_zlim()
    tickpos1,tickpos2=[xmax,xmax,xmin],[ymin,ymax,ymax]
    if(type(length_major)==type([1,1,1])):
        length=length_major[2]
    for tick,label in zip(ax.get_zticks(),ax.get_zticklabels()):
        if((tick<axislim[0]) | (tick>axislim[1])):
            continue
        ax.plot(tickpos1,tickpos2,[tick,tick,tick],
                lw=linewidth,color=color,zorder=zorder_grid)
        if(label3D):
            text3d(ax, (xmin-1.1*length*xlen,ymax,tick), label.get_text(),zdir='y',size=fs_label*minLen_axis,
                   fc=label.get_color(), ec="None",ha='center',va='bottom',angle=90)
            label.set_alpha(0)
            ax.plot([xmin,xmin-length*minLen_axis],[ymax,ymax],[tick,tick],lw=linewidth,color=label.get_color(),clip_on=False)
    if(label3D):
        ax.plot([xmin,xmin],[ymax,ymax],[zmin,zmax],lw=linewidth,color=label.get_color())
        text3d(ax, (xmin-length*offset_axislabel*xlen,ymax, (zmin+zmax)/2), ax.get_zlabel(),zdir='y',
               size=fs_label*minLen_axis, fc=label.get_color(), ec="None",ha='center',va='bottom',angle=90)
    if(frame_on==True):
        ax.plot([xmax,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymax],[zmin,zmax,zmax,zmax,zmin],color=label.get_color(),lw=linewidth)
        ax.plot([xmin,xmin,xmin],[ymin,ymax,ymax],[zmin,zmin,zmax],color=label.get_color(),lw=linewidth)
        ax.plot([xmin,xmax],[ymax,ymax],[zmin,zmin],color=label.get_color(),lw=linewidth,ls='dashed')
    # 副刻度线
    linewidth, color, length = lw_minor, color_minor,length_minor
    # x轴
    if(process_x):
        tickmin,tickmax,tickaxis,axislim=xmin,xmax,ax.xaxis,ax.get_xlim()
        # log段
        if(xMajor_loc_log==1):
            xmin_data,xmax_data = np.log10(data_lim[0]),np.log10(xcenter)
            xMajorLocator = mpl.ticker.MultipleLocator(xMajor_loc_log).tick_values(xmin_data,xmax_data)
            # 计算tick在坐标轴上的真正位置，将数据范围投影到坐标轴范围中
            xmin_axis,xmax_axis=ax.get_xlim()[0],0
            length_ = xmax_axis- xmin_axis
            scale_data2axes = length_/(xmax_data - xmin_data)
            tickpos1,tickpos2=[ymin,ymax,ymax],[zmin,zmin,zmax]
            for tick_data in xMajorLocator:
                for ii in np.log10(range(1,10)):
                    minortick_data = tick_data+ii
                    if((minortick_data<xmin_data) | (minortick_data>xmax_data)):
                        continue
                    tick = (minortick_data-xmin_data)*scale_data2axes + xmin_axis
                    # ax.plot([tick,tick,tick],tickpos1,tickpos2, lw=linewidth,color=color,zorder=zorder_grid)
                    ax.plot([tick,tick],[ymin,ymin-length*ylen],[zmin,zmin],lw=linewidth,color='k')
        # linear段
        xmin_data,xmax_data = xcenter,data_lim[1]
        xMajorLocator = mpl.ticker.MultipleLocator(xMajor_loc_linear).tick_values(xmin_data,xmax_data)
        xMinorLocator = mpl.ticker.MultipleLocator(xMinor_loc_linear).tick_values(xmin_data,xmax_data)
        # 计算tick在坐标轴上的真正位置，将数据范围投影到坐标轴范围中
        xmin_axis,xmax_axis=0,ax.get_xlim()[1]
        length_ = xmax_axis- xmin_axis
        scale_data2axes = length_/(xmax_data - xmin_data)
        tickpos1,tickpos2=[ymin,ymax,ymax],[zmin,zmin,zmax]
        for tick_data in xMinorLocator:
            if((tick_data<xmin_data) | (tick_data>xmax_data) | (tick_data in xMajorLocator)):
                continue
            tick = (tick_data-xmin_data)*scale_data2axes + xmin_axis
            ax.plot([tick,tick,tick],tickpos1,tickpos2, lw=linewidth,color=color,zorder=zorder_grid)
            if(label3D):
                ax.plot([tick,tick],[ymin,ymin-length*ylen],[zmin,zmin],lw=linewidth,color='k')

    # y轴
    tickmin,tickmax,tickaxis,axislim=ymin,ymax,ax.yaxis,ax.get_ylim()
    tickpos1,tickpos2=[xmax,xmax,xmin],[zmax,zmin,zmin]
    for tick in tickaxis.get_minor_locator().tick_values(tickmin,tickmax):
        if((tick<axislim[0]) | (tick>axislim[1])):
            continue
        ax.plot(tickpos1,[tick,tick,tick],tickpos2,
                lw=linewidth,color=color,zorder=zorder_grid)
        if(label3D):
            ax.plot([xmin,xmin-length*xlen],[tick,tick],[zmin,zmin],lw=linewidth,color=label.get_color())
    # z轴
    tickmin,tickmax,tickaxis,axislim=zmin,zmax,ax.zaxis,ax.get_zlim()
    tickpos1,tickpos2=[xmax,xmax,xmin],[ymin,ymax,ymax]
    for tick in tickaxis.get_minor_locator().tick_values(tickmin,tickmax):
        if((tick<axislim[0]) | (tick>axislim[1])):
            continue
        ax.plot(tickpos1,tickpos2,[tick,tick,tick],
                lw=linewidth,color=color,zorder=zorder_grid)
        if(label3D):
            ax.plot([xmin,xmin-length*xlen],[ymax,ymax],[tick,tick],lw=linewidth,color=label.get_color(),clip_on=False)
    # 3. 设置ticks和labels以及spines和panes
    for axis,inward_factor in zip([ax.xaxis,ax.yaxis,ax.zaxis],[0.3,0.2,0.2]):
        axis._axinfo['tick']['inward_factor'] = inward_factor
        axis._axinfo['tick']['outward_factor'] = 0
    ax.tick_params(which='minor', color=(0,0,1,0) ) # 取消副刻度
    # xaxis
    if(process_x):
        color,alpha,pad,axis,waxis=(1,0,0),alpha_pane,-3,ax.xaxis,ax.w_xaxis
        ax.tick_params(axis='x', which='major', pad=pad)
        axis.set_pane_color(color+(alpha,))
        axis.pane.set_edgecolor(ec_pane)
        axis.pane.fill = fill_pane
        waxis.line.set_color(ax.get_xticklines()[0].get_color())
    # yaxis
    color,alpha,pad,axis,waxis=(0,1,0),alpha_pane,-3,ax.yaxis,ax.w_yaxis
    ax.tick_params(axis='y', which='major', pad=pad)
    axis.set_pane_color(color+(alpha,))
    axis.pane.set_edgecolor(ec_pane)
    axis.pane.fill = fill_pane
    waxis.line.set_color(ax.get_yticklines()[0].get_color())
    # zaxis
    color,alpha,pad,axis,waxis=(0,0,1),alpha_pane,0,ax.zaxis,ax.w_zaxis
    ax.tick_params(axis='z', which='major', pad=pad)
    axis.set_pane_color(color+(alpha,))
    axis.pane.set_edgecolor(ec_pane)
    axis.pane.fill = fill_pane
    waxis.line.set_color(ax.get_zticklines()[0].get_color())
    if(label3D):
        ax.axis('off')

def set_axis_diagram_3D(ax,xlim=(0,100),ylim=(0,1000),zlim=(0,2500),xMajor_loc=10,xMinor_loc=2,yMajor_loc=100,yMinor_loc=20,zMajor_loc=500,zMinor_loc=100,xlabel="Wt.% NaCl",ylabel="Temperature ($^{\circ}$C)",zlabel="Pressure (bar)"):
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_zlabel(zlabel)

    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    ax.set_zlim(zlim)
    # 设置坐标轴刻度间隔
    ax.xaxis.set_major_locator(MultipleLocator(xMajor_loc))
    ax.xaxis.set_minor_locator(MultipleLocator(xMinor_loc))
    ax.yaxis.set_major_locator(MultipleLocator(yMajor_loc))
    ax.yaxis.set_minor_locator(MultipleLocator(yMinor_loc))
    ax.zaxis.set_major_locator(MultipleLocator(zMajor_loc))
    ax.zaxis.set_minor_locator(MultipleLocator(zMinor_loc))

    # 重新自定义坐标轴属性
    niceAxis_3D(ax,fill_pane=False,label3D=True,fs_label=0.024,length_major=0.02,length_minor=0.01, scaled=False)
    # ax.legend()
    ax.view_init(elev=25, azim=-145)

def set_axis_diagram_3D_loglinearx(ax,xcenter=1,ratio_log_lin=(1,1),xlim=(1E-14,100),ylim=(0,1000),zlim=(0,2500),xMajor_loc_log=10,xMinor_loc_log=2,xMajor_loc_linear=10,xMinor_loc_linear=2,yMajor_loc=100,yMinor_loc=20,zMajor_loc=500,zMinor_loc=100,xlabel="Wt.% NaCl",ylabel="Temperature ($^{\circ}$C)",zlabel="Pressure (bar)"):
    ax.set_box_aspect((8,4,3))
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_zlabel(zlabel)
    ax.set_ylim(ylim)
    # 把x轴的范围投影到-1，1之间
    data_lim = xlim
    ax.set_xlim((-ratio_log_lin[0]/ratio_log_lin[1],1))
    ax.set_zlim(zlim)
    # 设置坐标轴刻度间隔
    # ax.xaxis.set_major_locator(MultipleLocator(xMajor_loc))
    # ax.xaxis.set_minor_locator(MultipleLocator(xMinor_loc))
    ax.yaxis.set_major_locator(MultipleLocator(yMajor_loc))
    ax.yaxis.set_minor_locator(MultipleLocator(yMinor_loc))
    ax.zaxis.set_major_locator(MultipleLocator(zMajor_loc))
    ax.zaxis.set_minor_locator(MultipleLocator(zMinor_loc))

    # 重新自定义坐标轴属性
    # 设置yz轴
    niceAxis_3D_loglinearx(ax,data_lim,xcenter=xcenter,fill_pane=False,label3D=True,fs_label=0.024,length_major=0.02,length_minor=0.01, scaled=False, process_x=False)
    # 专门处理x轴
    ax.set_xlabel(xlabel)
    niceAxis_3D_loglinearx(ax,data_lim,xcenter=xcenter,xMajor_loc_log=xMajor_loc_log,xMinor_loc_log=xMinor_loc_log,xMajor_loc_linear=xMajor_loc_linear,xMinor_loc_linear=xMinor_loc_linear,
                           fill_pane=False,label3D=True,fs_label=0.024,length_major=0.02,length_minor=0.01, scaled=False, process_x=True)
    # ax.legend()
    ax.view_init(elev=25, azim=-145)

    return {'xlim':xlim,'xcenter':xcenter,'ratio_log_lin':ratio_log_lin}

# transform original data to log-linear axis data
def data2axis_loglin(axtrans,x):
    x_new = x*0
    if(type(x)==type(np.linspace(0,1,2))):   # array
        ind_log=(x<=axtrans['xcenter'])
        ind_linear=(x>axtrans['xcenter'])
        # linear part
        xmin=axtrans['xcenter']
        xmax=axtrans['xlim'][1]
        x_new[ind_linear] = (x[ind_linear]-xmin)/(xmax-xmin)
        # log part
        xmin=np.log10(axtrans['xlim'][0])
        xmax=np.log10(axtrans['xcenter'])
        x_log=np.log10(x)
        x_new[ind_log] = (x_log[ind_log]-xmax)/(xmax-xmin)
    else:                                    # number
        if(x<=axtrans['xcenter']):
            # log part
            xmin=np.log10(axtrans['xlim'][0])
            xmax=np.log10(axtrans['xcenter'])
            x_log=np.log10(x)
            x_new = (x_log-xmax)/(xmax-xmin)
        elif(x>axtrans['xcenter']):
            # linear part
            xmin=axtrans['xcenter']
            xmax=axtrans['xlim'][1]
            x_new = (x-xmin)/(xmax-xmin)
    return x_new

def plot_coloredline(ax,x,y,data,cmap='YlGnBu',lw=4):
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    norm = plt.Normalize(data.min(), data.max())
    lc = LineCollection(segments, cmap=cmap, norm=norm)
    lc.set_array(data)
    lc.set_linewidth(lw)
    line = ax.add_collection(lc)
    return line

def plot_coloredline(ax,x,y,data,vmin=None,vmax=None,cmap='YlGnBu',lw=4):
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    if(vmin is None): vmin = data.min()
    if(vmax is None): vmax = data.max()
    norm = plt.Normalize(vmin,vmax)
    lc = LineCollection(segments, cmap=cmap, norm=norm)
    lc.set_array(data)
    lc.set_linewidth(lw)
    line = ax.add_collection(lc)
    return line

