import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import re
import astropy.units as u

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
                     figsize=(10, 15), dpi=200, FBS=None, datamin=None, datamax=None, mds_offset_cm=0):
    
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
    sort_order = np.argsort(np.abs(unsort_mds))
    mds = np.sort(np.abs(unsort_mds)) + mds_offset_cm

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
        if data_func is not None:
            data = data_func(data)
        if datamin is not None:
            data = data[data>=datamin]
        if datamax is not None:
            data = data[data<=datamax]

        if type(color_map) is list:
            j = int(k/n_chunks)
            c = color_map[j](Norm[j](np.abs(unsort_mds[order])))
        else:
            c = color_map(Norm(np.abs(unsort_mds[order])))
        _ = ax.hist(data, bins=bins, histtype='step', color=c, \
                 density=density, label=run)
    
    # label & legend
    ax.set_xlabel(xlabel, fontsize=12)
    ncol_legend = 1
    if len(mds)>40:
        ncol_legend = 2
    ax.legend(fontsize=7.5, bbox_to_anchor=(1.0, 1.02), edgecolor='k', loc=2, labelspacing=0.45, ncol=ncol_legend)
    
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
    ax.legend(fontsize=7.5, bbox_to_anchor=(1.0, 1.02), edgecolor='k', loc=2, labelspacing=0.45, ncol=2)
    
    ax.set_ylabel(ylabel, fontsize=12, labelpad=7)
        
    #Set the xlabel range.
    xmin, xmax = ax.get_xlim()
    if datamin is not None:
        xmin = datamin
    if datamax is not None:
        xmax = datamax
    ax.set_xlim([xmin,xmax])
