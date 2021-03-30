import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import re
import astropy.units as u
from textwrap import wrap

import matplotlib._color_data as mcd

#Routines adapted from the DRC plot routines from Gordon Richard's group DCR RNAAS paper at: https://github.com/RichardsGroup/LSST_DCR/blob/main/notebooks/02_RNAAS.ipynb

# return median metric data
def get_Opsim_dist_median(mb, data_func=None):
    
    mask = mb.metricValues.mask
    data = mb.metricValues.data[~mask]
    data = data[~(np.isnan(data) | np.isinf(data))]
    if data_func is not None:
        data = data_func(data)
    
    return np.median(data)

# get the median values from all opsims 
# for normaliation in plotting
def get_metric_medians(key, bd, data_func=None):
    
    mds = []
    for run in bd:
        keys = [*bd[run].keys()]
        run_key = [elem for elem in keys if elem[1] == key[1]][0]
        mds.append(get_Opsim_dist_median(bd[run][run_key], data_func))

    return mds

def plot_OpSims_hist(Key, bundleDicts_input, order_func=get_metric_medians, data_func=None, 
                     color_map=mpl.cm.summer, xlabel=None, healpix_pixarea=6.391586616190171e-05*u.sr, 
                     figsize=(10, 15), dpi=200, FBS=None, datamin=None, datamax=None, title=None,
                     bins=60):
    
    #First, select the runs to use by FBS version if requested.
    if FBS is None:
        bundleDicts = bundleDicts_input
    else:
        bundleDicts = dict()
        all_runs = list(bundleDicts_input.keys())
        for run in all_runs:
            if not re.search(FBS,run):
                continue
            bundleDicts[run] = bundleDicts_input[run]
    
    #If no metrics, do not continue. 
    if len(list(bundleDicts.keys()))==0:
        return
    
    # get plotting order
    unsort_mds = order_func(Key, bundleDicts, data_func)
    runs = list(bundleDicts.keys())
    sort_order = np.argsort(unsort_mds)
    mds = np.sort(unsort_mds)
    mds_offset = mds[0] - 1e-3*(mds[-1]-mds[0])
    mds -= mds_offset

    #Print the names of the extreme metrics according to the sorting function. 
    print(runs[sort_order[ 0]], unsort_mds[sort_order[ 0]])
    print(runs[sort_order[-1]], unsort_mds[sort_order[-1]])
    
    # create normalization object. If more than one color map is being used, set the normalization as a list.
    if type(color_map) is list:
        Norm = list()
        n_mds = len(mds)
        n_chunks = int(np.ceil(len(mds)/len(color_map)))
        for i in range(len(color_map)):
            j1 = i*n_chunks
            j2 = (i+1)*n_chunks
            if j2>=len(mds):
                j2 = len(mds)-1
            Norm.append(mpl.colors.LogNorm(vmin=mds[j1], vmax=mds[j2]))
    else:
        Norm = mpl.colors.LogNorm(vmin=mds[0], vmax=mds[-1])

    # other plot setting
    density = False
    #bins = 60

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    for k, order in enumerate(sort_order):
    
        run = runs[order]

        # need to mask the pixels that have no available data
        mask = bundleDicts[run][Key].metricValues.mask
        data = bundleDicts[run][Key].metricValues.data[~mask]
        data = data[~(np.isnan(data) | np.isinf(data))]
        if data_func is not None:
            data = data_func(data)
        if datamin is not None:
            data = data[data>=datamin]
        if datamax is not None:
            data = data[data<=datamax]

        if type(color_map) is list:
            j = int(k/n_chunks)
            c = color_map[j](Norm[j](unsort_mds[order]-mds_offset))
        else:
            c = color_map(Norm(unsort_mds[order]-mds_offset))
        _ = ax.hist(data, bins=bins, histtype='step', color=c, \
                 density=density, label=run)
    
    # label & legend
    ax.set_xlabel(xlabel, fontsize=12)
    #ncol_legend = 1
    #if len(mds)>40:
    #    ncol_legend = 2
    ncol_legend = 1 + int(len(mds)/60.)
    ax.legend(fontsize=7.5, bbox_to_anchor=(1.0, 1.0), edgecolor='k', loc=2, labelspacing=0.45, ncol=ncol_legend)
    
    #ax.yaxis.set_major_locator(plt.FixedLocator(np.array([500, 1000, 1500, 2000])/(healpix_pixarea/60)**2))
    y_vals = ax.get_yticks()
    ax.set_yticklabels(['{:.0f}'.format(x * healpix_pixarea.to(u.deg**2).value) for x in y_vals], rotation=90)
    ax.set_ylabel('Area ($\mathrm{degree^{2}}$)', labelpad=7)
        
    #Set the xlabel range.
    xmin, xmax = ax.get_xlim()
    if datamin is not None:
        xmin = datamin
    if datamax is not None:
        xmax = datamax
    ax.set_xlim([xmin,xmax])
    
    #Add the title if provided.
    if title is not None:
        ax.set_title(title)
    
####

def get_data(bundleDicts, run, Key, data_func=None, datamin=None, datamax=None):
    
    # need to mask the pixels that have no available data
    mask = bundleDicts[run][Key].metricValues.mask
    data = bundleDicts[run][Key].metricValues.data[~mask]
    data = data[~(np.isnan(data) | np.isinf(data))]
    if data_func is not None:
        data = data_func(data)
    if datamin is not None:
        data = data[data>=datamin]
    if datamax is not None:
        data = data[data<=datamax]
    return data
    

def plot_OpSims_hist_extremes(Key, bundleDicts_input, order_func=get_metric_medians, 
                              data_func=None, color_map_top=mpl.cm.summer_r, 
                              color_map_bottom=mpl.cm.cool_r, xlabel=None, 
                              healpix_pixarea=6.391586616190171e-05*u.sr, 
                              figsize=(10, 15), dpi=200, FBS=None, datamin=None, 
                              datamax=None, title=None, bins=60, ymin_use=None, ymax_use=None, 
                              percentile=10, top_axis=False, data_func_top=None, 
                              top_xlabel=None, survey_label=None, legend_box_in_plot=False):
    
    #First, select the runs to use by FBS version if requested.
    if FBS is None:
        bundleDicts = bundleDicts_input
    else:
        bundleDicts = dict()
        all_runs = list(bundleDicts_input.keys())
        for run in all_runs:
            if not re.search(FBS,run):
                continue
            bundleDicts[run] = bundleDicts_input[run]
    
    #If no metrics, do not continue. 
    if len(list(bundleDicts.keys()))==0:
        return
    
    # get plotting order
    unsort_mds = order_func(Key, bundleDicts, data_func)
    runs = list(bundleDicts.keys())
    sort_order = np.argsort(unsort_mds)
    mds = np.sort(unsort_mds)
    mds_offset = mds[0] - 1e-3*(mds[-1]-mds[0])
    mds -= mds_offset

    #Print the names of the extreme metrics according to the sorting function. 
    print(runs[sort_order[ 0]], unsort_mds[sort_order[ 0]])
    print(runs[sort_order[-1]], unsort_mds[sort_order[-1]])
    
    #Create normalization objects. We use the top and bottom percentile requested.
    pcs = [0.0, percentile, 100.-percentile, 100.]
    mds_lims = np.percentile(mds,pcs)
    Norm = [None]*2
    #Norm[0] = mpl.colors.LogNorm(vmin=mds_lims[0], vmax=mds_lims[1])
    #Norm[1] = mpl.colors.LogNorm(vmin=mds_lims[2], vmax=mds_lims[3])
    Norm[0] = mpl.colors.Normalize(vmin=mds_lims[0], vmax=mds_lims[1])
    Norm[1] = mpl.colors.Normalize(vmin=mds_lims[2], vmax=mds_lims[3])
    
    # other plot setting
    density = False

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    #First draw only the ones in gray, not in the top and bottom percentiles.
    for k, order in enumerate(sort_order):  
        run = runs[order]
        data = get_data(bundleDicts,run,Key,data_func,datamin,datamax)
        mds_aux = unsort_mds[order]-mds_offset
        if  mds_aux>mds_lims[1] and mds_aux<mds_lims[2]:
            c = mcd.XKCD_COLORS["xkcd:light grey"].upper()
        else:
            continue
        _ = ax.hist(data, bins=bins, histtype='step', color=c, density=density, label=None)

    #Now draw the bottom percentile. Do it backwards so that the bottom run is drawn last.
    for k, order in enumerate(sort_order[::-1]):
        run = runs[order]
        data = get_data(bundleDicts,run,Key,data_func,datamin,datamax)
        mds_aux = unsort_mds[order]-mds_offset
        if mds_aux <= mds_lims[1]:
            c = color_map_bottom(Norm[0](unsort_mds[order]-mds_offset))
        else:
            continue
        _ = ax.hist(data, bins=bins, histtype='step', color=c, density=density, label=run)
        
        
    #Finally, draw the top percentile.
    for k, order in enumerate(sort_order):  
        run = runs[order]
        data = get_data(bundleDicts,run,Key,data_func,datamin,datamax)
        mds_aux = unsort_mds[order]-mds_offset
        if mds_aux >= mds_lims[2]: 
            c = color_map_top(Norm[1](unsort_mds[order]-mds_offset))
        else:
            continue
        _ = ax.hist(data, bins=bins, histtype='step', color=c, density=density, label=run)

    #Get all the legends and then sort them. 
    handles_raw, labels_raw = ax.get_legend_handles_labels()
    labels_raw = np.array(labels_raw)
    handles = list()
    labels = list()
    for order in sort_order:
        run = runs[order]
        if run in labels_raw:
            k = np.argwhere(labels_raw==run)[0][0]
            labels.append('\n'.join(wrap(labels_raw[k],35)))
            handles.append(handles_raw[k])
    
    # label & legend
    ax.set_xlabel(xlabel, fontsize=12)
    ncol_legend = 1
    #ncol_legend = 1 + int(len(mds[])/60.)
    if legend_box_in_plot:
        bbox_to_anchor=(0.03, 0.90)
        fontsize=5.0
    else:
        bbox_to_anchor=(1.0, 1.0)
        fontsize=7.5
    ax.legend(handles, labels, fontsize=fontsize, bbox_to_anchor=bbox_to_anchor, 
              edgecolor='k', loc=2, labelspacing=0.45, ncol=ncol_legend)
    
    #Set the y-limit range and then the y ticks. 
    ymin, ymax = ax.get_ylim()
    if ymin_use is not None:
        ymin = ymin_use/healpix_pixarea.to(u.deg**2).value
    if ymax_use is not None:
        ymax = ymax_use/healpix_pixarea.to(u.deg**2).value
    ax.set_ylim([ymin, ymax])
    
    y_vals = ax.get_yticks()
    pixarea = healpix_pixarea.to(u.deg**2).value
    ax.set_yticklabels(['{:.0f}'.format(x * pixarea) for x in y_vals], rotation=90)
    ax.set_ylabel('Area ($\mathrm{degree^{2}}$)', fontsize=12, labelpad=7)
        
    #Set the xlabel range.
    xmin, xmax = ax.get_xlim()
    if datamin is not None:
        xmin = datamin
    if datamax is not None:
        xmax = datamax
    ax.set_xlim([xmin,xmax])
    
    #Add the title if provided.
    if title is not None:
        ax.set_title(title)
    
    if top_axis:
        ax2 = ax.twiny()
        xlim = ax.get_xlim()
        if data_func_top is not None:
            xlim = data_func_top(xlim)
        ax2.set_xlim(xlim)
        #ax2.set_xticks(new_tick_locations)
        ax2.set_xlabel(top_xlabel, fontsize=12, labelpad=7)

    if survey_label is not None:
        if legend_box_in_plot:
            xloc = 0.17
            yloc = 0.85
        else:
            xloc = 0.125
            yloc = 0.85
        ax.annotate(survey_label,
            xy=(xloc, yloc), xycoords='figure fraction',
            horizontalalignment='left', verticalalignment='top',
            fontsize=14)
  
####

def plot_OpSims_color_excess_redshift_extremes(Key, bundleDicts_input, zs, quasar_colors, 
                                               order_func=get_metric_medians, 
                                               color_map_top=mpl.cm.summer_r,
                                               color_map_bottom=mpl.cm.cool_r,
                                               ylabel=None, figsize=(10, 15), dpi=200, 
                                               FBS=None, datamin=None, datamax=None,
                                               percentile=10):
    
    #First, select the runs to use by FBS version if requested.
    if FBS is None:
        bundleDicts = bundleDicts_input
    else:
        bundleDicts = dict()
        all_runs = list(bundleDicts_input.keys())
        for run in all_runs:
            if not re.search(FBS,run):
                continue
            bundleDicts[run] = bundleDicts_input[run]
    
    # get plotting order
    unsort_mds = order_func(Key, bundleDicts)
    runs = list(bundleDicts.keys())
    sort_order = np.argsort(unsort_mds)
    mds = np.sort(unsort_mds)
    mds_offset = mds[0] - 1e-3*(mds[-1]-mds[0])
    mds -= mds_offset

    #Print the names of the extreme metrics according to the sorting function. 
    print(runs[sort_order[ 0]], unsort_mds[sort_order[ 0]])
    print(runs[sort_order[-1]], unsort_mds[sort_order[-1]])
    
    
    #Create normalization objects. We use the top and bottom percentile requested.
    pcs = [0.0, percentile, 100.-percentile, 100.]
    mds_lims = np.percentile(mds,pcs)
    Norm = [None]*2
    #Norm[0] = mpl.colors.LogNorm(vmin=mds_lims[0], vmax=mds_lims[1])
    #Norm[1] = mpl.colors.LogNorm(vmin=mds_lims[2], vmax=mds_lims[3])
    Norm[0] = mpl.colors.Normalize(vmin=mds_lims[0], vmax=mds_lims[1])
    Norm[1] = mpl.colors.Normalize(vmin=mds_lims[2], vmax=mds_lims[3])

    # other plot setting
    density = False

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    #First, plot all the grey ones.
    for k, order in enumerate(sort_order):  
        run = runs[order]
        data = get_data(bundleDicts, run, Key)
        mds_aux = unsort_mds[order]-mds_offset
        if mds_aux > mds_lims[1] and mds_aux < mds_lims[2]:
            c = mcd.XKCD_COLORS["xkcd:light grey"].upper()
        else:
            continue
        color_cond = (~np.isnan(quasar_colors)) & (~np.isinf(quasar_colors))
        color_excess_z = unsort_mds[order] - quasar_colors[color_cond]
        zs_use = zs[color_cond]
        ax.plot(zs_use, color_excess_z, color=c, label=None)

    #Now, plot the ones in the top and bottom percentiles.
    for k, order in enumerate(sort_order):  
        run = runs[order]
        data = get_data(bundleDicts, run, Key)
        mds_aux = unsort_mds[order]-mds_offset
        if mds_aux <= mds_lims[1]:
            c = color_map_bottom(Norm[0](mds_aux))
        elif mds_aux >= mds_lims[2]: 
            c = color_map_top(Norm[1](mds_aux))
        else:
            continue
        color_cond = (~np.isnan(quasar_colors)) & (~np.isinf(quasar_colors))
        color_excess_z = unsort_mds[order] - quasar_colors[color_cond]
        zs_use = zs[color_cond]
        ax.plot(zs_use, color_excess_z, color=c, label=run)

        
    # label & legend
    ax.set_xlabel("Redshift", fontsize=12)
    ncol_legend = 1
    ax.legend(fontsize=7.5, bbox_to_anchor=(1.0, 1.0), edgecolor='k', 
              loc=2, labelspacing=0.45, ncol=ncol_legend)
    ax.set_ylabel(ylabel, fontsize=12, labelpad=7)
        
    #Set the xlabel range.
    xmin, xmax = ax.get_xlim()
    if datamin is not None:
        xmin = datamin
    if datamax is not None:
        xmax = datamax
    ax.set_xlim([xmin,xmax])
    
    return

####

def plot_OpSims_color_excess_redshift(Key, bundleDicts_input, zs, quasar_colors, order_func=get_metric_medians, 
                               color_map=mpl.cm.summer, ylabel=None, figsize=(10, 15), dpi=200, 
                               FBS=None, datamin=None, datamax=None):
    
    #First, select the runs to use by FBS version if requested.
    if FBS is None:
        bundleDicts = bundleDicts_input
    else:
        bundleDicts = dict()
        all_runs = list(bundleDicts_input.keys())
        for run in all_runs:
            if not re.search(FBS,run):
                continue
            bundleDicts[run] = bundleDicts_input[run]
    
    # get plotting order
    unsort_mds = order_func(Key, bundleDicts)
    runs = list(bundleDicts.keys())
    sort_order = np.argsort(unsort_mds)
    mds = np.sort(unsort_mds)
    mds_offset = mds[0] - 1e-3*(mds[-1]-mds[0])
    mds -= mds_offset

    #Print the names of the extreme metrics according to the sorting function. 
    print(runs[sort_order[ 0]], unsort_mds[sort_order[ 0]])
    print(runs[sort_order[-1]], unsort_mds[sort_order[-1]])
    
    # create normalization object. If more than one color map is being used, set the normalization as a list.
    if type(color_map) is list:
        Norm = list()
        n_mds = len(mds)
        n_chunks = int(np.ceil(len(mds)/len(color_map)))
        for i in range(len(color_map)):
            j1 = i*n_chunks
            j2 = (i+1)*n_chunks
            if j2>=len(mds):
                j2 = len(mds)-1
            Norm.append(mpl.colors.LogNorm(vmin=mds[j1], vmax=mds[j2]))
    else:
        Norm = mpl.colors.LogNorm(vmin=mds[0], vmax=mds[-1])

    # other plot setting
    density = False
    bins = 60

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    for k, order in enumerate(sort_order):
    
        run = runs[order]

        # need to mask the pixels that have no available data
        mask = bundleDicts[run][Key].metricValues.mask
        data = bundleDicts[run][Key].metricValues.data[~mask]
        data = data[~(np.isnan(data) | np.isinf(data))]
        
        if type(color_map) is list:
            j = int(k/n_chunks)
            c = color_map[j](Norm[j](unsort_mds[order]-mds_offset))
        else:
            c = color_map(Norm(unsort_mds[order]-mds_offset))
            
        color_cond = (~np.isnan(quasar_colors)) & (~np.isinf(quasar_colors))
        color_excess_z = unsort_mds[order] - quasar_colors[color_cond]
        zs_use = zs[color_cond]
        
        ax.plot(zs_use, color_excess_z, color=c, label=run)
    
    # label & legend
    ax.set_xlabel("Redshift", fontsize=12)
    ncol_legend = 1 + int(len(mds)/60.)
    ax.legend(fontsize=7.5, bbox_to_anchor=(1.0, 1.0), edgecolor='k', loc=2, labelspacing=0.45, ncol=ncol_legend)
    
    ax.set_ylabel(ylabel, fontsize=12, labelpad=7)
        
    #Set the xlabel range.
    xmin, xmax = ax.get_xlim()
    if datamin is not None:
        xmin = datamin
    if datamax is not None:
        xmax = datamax
    ax.set_xlim([xmin,xmax])

####

def plot_OpSims_Nqso_hist(Key, bundleDicts_input, xlabel=None, figsize=(10, 15), dpi=200, FBS=None, datamin=None, datamax=None, bins=60, title=None):
    
    #First, select the runs to use by FBS version if requested.
    if FBS is None:
        bundleDicts = bundleDicts_input
    else:
        bundleDicts = dict()
        all_runs = list(bundleDicts_input.keys())
        for run in all_runs:
            if not re.search(FBS,run):
                continue
            bundleDicts[run] = bundleDicts_input[run]
    
    #Determine the number of quasars found.
    dbRuns = list(bundleDicts.keys())
    Nqso   = np.zeros(len(dbRuns))
    for k, run in enumerate(dbRuns):
        mb  = bundleDicts[run][Key]
        mask = mb.metricValues.mask
        area_sq_deg = (mb.slicer.pixArea * u.sr).to(u.deg**2).value
        Nqso[k] = np.sum(mb.metricValues[~mask]*area_sq_deg)
    
    #Make a histogram of Nqso. 
    density = False
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    
    _ = ax.hist(Nqso, bins=bins, histtype='step', density=density)
    
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(r'Number of OpSim Runs', labelpad=7)
    
    #Set the xlabel range.
    xmin, xmax = ax.get_xlim()
    if datamin is not None:
        xmin = datamin
    if datamax is not None:
        xmax = datamax
    ax.set_xlim([xmin,xmax])
    
    #Add the title if provided.
    if title is not None:
        ax.set_title(title)

###

def plot_OpSims_Nqso_hist_v2(Key, Nqso, xlabel=None, figsize=(5, 10), dpi=200,
                             datamin=None, datamax=None, bins=60, title=None):
      
    #Make a histogram of Nqso. 
    density = False
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    
    _ = ax.hist(Nqso[Key], bins=bins, histtype='step', density=density)
    
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(r'Number of OpSim Runs', labelpad=7)
    
    #Set the xlabel range.
    xmin, xmax = ax.get_xlim()
    if datamin is not None:
        xmin = datamin
    if datamax is not None:
        xmax = datamax
    ax.set_xlim([xmin,xmax])
    
    #Add the title if provided.
    if title is not None:
        ax.set_title(title)

    return