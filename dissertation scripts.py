

def xHII_splatter(snapshot):
    import matplotlib as mpl
    from matplotlib import ticker
    import matplotlib.gridspec as gridspec
    import matplotlib.pylab as plt
    import numpy as np
    from seren3.analysis.visualization import engines, operators
    from seren3.utils import plot_utils

    op = operators.MassWeightedOperator("xHII", snapshot.C.none)
    eng = engines.CustomSplatterEngine(snapshot.g, 'xHII', op, extra_fields=['rho'])
    cam = snapshot.camera()
    cam.map_max_size = 2048
    proj = eng.process(cam)
    im = proj.save_plot()
    im.show()

import seren3
snapshot=seren3.load_snapshot('/lustre/scratch/astro/dn99/bc03',00050)


def spherical_filter(snapshot, center=None, radius=None):
 
    if (center is None):
        center = [0.5, 0.5, 0.5]  # center in usual box units of [0., 1.)]

    if (radius is None):
        radius = 0.005

    sphere = snapshot.get_sphere(center, radius)

    sub_snapshot = snapshot[sphere]

    return sub_snapshot


sub_snapshot = spherical_filter(snapshot, center=None, radius=None)
dset = sub_snapshot.g[['dx','rho','xHII']].flatten()
#can alter snap to be much smaller to allow for error testing much more rapidly









def load_xHII_table(path='./'):
    import pickle, os
    fname = "{PATH}/xHII_reion_history.p".format(PATH=path)
    if os.path.isfile(fname):
        data = pickle.load( open(fname, "rb") )
        table = {}
        for item in data:
            table[item.idx] = item.result
        return table
    else:
        raise IOError("No such file: {FNAME}".format(FNAME=fname))

path = '/its/home/dn99/scratch/bc03/pickle/'

a = load_xHII_table(path)

i = 0
z=[]
m=[]
v=[]
while i < len(a):
    i = i + 1
    z.append(a[i]['z'])
    m.append(a[i]['mass_weighted'])
    v.append(a[i]['volume_weighted'])

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


plt.title('reionization history for bc03 simulation')
plt.xlabel('redshift')
plt.ylabel('Ionized fraction, xHII')
axes = plt.gca()
axes.set_ylim([0,1.1])
axes.set_xlim([5.9,14])
plt.plot(z,v,'r')
plt.plot(z,m,'b')
red = mpatches.Patch(color='red', label = 'Volume weighted ionized fraction')
blue = mpatches.Patch(color='blue', label = 'Mass weighted ionized fraction')
plt.legend(handles=[red,blue])

plt.show()





def load(snap):
    import pickle

    data = pickle.load( open("%s/pickle/ConsistentTrees/time_int_fesc_all_halos_%05i.p" % (snap.path, snap.ioutput), "rb") )
    print len(data)
    return data

def load_halo(halo):
    data = load(halo.base)

    for i in range(len(data)):
        if (int(data[i].idx) == int(halo["id"])):
            return data[i].result
    return None

import seren3
import numpy as np
snapshot = seren3.load_snapshot('/research/prace/david/bpass/bc03',00100)

a = load(snapshot)
i = 0

#while loop helps to find halo with most data associated to it
while i < len(a):
    num = len(a[i].result['fesc'])
    m = a[i].result['Mvir']
    mass = np.log10(m)
    halo_id = a[i].idx
    if num > 18:
        print 'halo id is %s, number of entries is %s, Mvir of halo is %s, Output_ID is %s' %(halo_id,num,mass,i)
        i=i+1
    else:
        i=i+1
        
        
        
        
def plot(snap, idata, data=None):
    from matplotlib import rcParams
    rcParams['figure.figsize'] = 16, 6
    rcParams['axes.linewidth'] = 1.5
    rcParams['xtick.labelsize'] = 14
    rcParams['ytick.labelsize'] = 14
    rcParams['axes.labelsize'] = 20
    rcParams['xtick.major.pad'] = 10
    rcParams['ytick.major.pad'] = 10

    import matplotlib.pylab as plt

    if (data is None):
        data = load(snap)
    item = data[idata]

    halos = snap.halos(finder="ctrees")
    hid = int(item.idx)
    h = halos.with_id(hid)

    print "Mvir = %1.2e [Msol/h]" % h["Mvir"]

    tint_fesc = item.result["tint_fesc_hist"]
    fesc = item.result["fesc"]
    lbtime = item.result["lbtime"]

    print (fesc*100.).min(), (fesc*100.).max()

    from seren3.analysis.stars import sfr as sfr_fn

    nbins=75
    agerange = [0, lbtime.max()]
    sSFR, SFR, sSFR_lookback_time, sSFR_binsize = sfr_fn(h, ret_sSFR=True, nbins=nbins, agerange=agerange)

    ax1 = plt.gca()
    ax2 = ax1.twinx()
    lbtime_units = "Myr"
    ax1.semilogy(lbtime.in_units(lbtime_units), tint_fesc*100., color='red', linewidth=2.5, linestyle='--')
    ax1.semilogy(lbtime.in_units(lbtime_units), fesc*100., color='blue', linewidth=2.5)

    # sSFR
    sSFR_col = '#0099cc'
    ax2.fill_between(sSFR_lookback_time.in_units(lbtime_units), sSFR, color=sSFR_col, alpha=0.2)
    ax2.set_ylim(1e-2, 1.1e2)    #changed axis to 1e-2 from 1e-1
    ax2.set_yscale("log")

    for tl in ax2.get_yticklabels():
        tl.set_color(sSFR_col)
    ax2.set_ylabel(r"sSFR [Gyr$^{-1}$]")
    ax2.yaxis.label.set_color(sSFR_col)

    for ax in [ax1, ax2]:
        ax.tick_params('both', length=20, width=2, which='major')
        ax.tick_params('both', length=10, width=1, which='minor')

    ax1.set_ylabel(r"f$_{\mathrm{esc}}$ [%]")
    ax1.set_xlabel(r"Lookback-time [%s]" % lbtime_units)

    ax1.set_ylim(1e-1, 1.1e2)    

    plt.xlim(0., 300.)
    plt.tight_layout()
    plt.show()
    
    
    
    
    
def _plot(fesc, tot_mass, label, alpha, ax, c, nbins):
    import numpy as np
    import random
    # from seren3.analysis.plots import fit_scatter
    from seren3.analysis import plots

    reload(plots)

    print 'length of fesc is %s' % (len(fesc))

    # Deal with fesc>1. following Kimm & Cen 2014
    bad = np.where( fesc > 1. )
    nbad = float(len(bad[0]))
    print "%f %% of points above 1" % ((nbad/float(len(fesc)))*100.)
    for ix in bad:
        fesc[ix] =  fesc[ix].flatten()
        fesc[ix] = random.uniform(0.9,1.0)
    #    print fesc[ix]

    keep = np.where( np.logical_and(fesc >= 0., fesc <=1.) )
    keep_mass = np.where(np.logical_and(tot_mass>=10.**(7.5), tot_mass>=0)) #added 
    fesc = fesc[keep]
    fesc = fesc[keep_mass]
    tot_mass = tot_mass[keep]
    tot_mass = tot_mass[keep_mass]

    log_tot_mass = np.log10(tot_mass)
    log_fesc = np.log10(100. * fesc)

    fesc_percent = fesc * 100.

    remove = np.where( np.logical_or( np.isinf(log_fesc), np.isnan(log_fesc) ) )

    log_fesc = np.delete( log_fesc, remove )
    fesc_percent = np.delete( fesc_percent, remove )
    log_tot_mass = np.delete( log_tot_mass, remove )

    remove = np.where( log_fesc < -1. ) 

    log_fesc = np.delete( log_fesc, remove )
    fesc_percent = np.delete( fesc_percent, remove )
    log_tot_mass = np.delete( log_tot_mass, remove )

    print 'length of log_fesc is %s' % (len(log_fesc))

    bin_centres, mean, std_dev, std_err = plots.fit_scatter(log_tot_mass, fesc_percent, nbins=nbins, ret_sterr=True)

    # Plot
    ax.scatter(log_tot_mass, fesc_percent, marker='+', color=c, alpha=alpha)

    e = ax.errorbar(bin_centres, mean, yerr=std_dev, color=c, label=label,\
         fmt="o", markerfacecolor=c, mec='k', capsize=2, capthick=2, elinewidth=2, linestyle='-')

    # e = ax.plot(bin_centres, median, color=c, label=label, linestyle='-', linewidth=2.)

    ax.set_yscale("log")
    ax.set_xlim(7.5, 10.) 






def plot(paths, ioutputs, labels, nbins=4, alpha=0.75, ax=None, cols=None):
    import seren3
    import matplotlib.pylab as plt
    import numpy as np
    import pickle

    if ax is None:
        fig, ax = plt.subplots()

    if cols is None:
        from seren3.utils import plot_utils
        cols = plot_utils.ncols(len(paths))

    for path, iout, label, c in zip(paths, ioutputs, labels, cols):
        
        # fname = "%s/pickle/ConsistentTrees/fesc_filt_%05i.p" % (path, iout)
      
        fname = "%s/pickle/fesc_do_multigroup_filt_no_age_limit%05i.p" % (path, iout)   # added by Dylan
        with open(fname, "rb") as f:
            dset = pickle.load( f )

            fesc = np.zeros(len(dset)); tot_mass = np.zeros(len(dset))
            for i in range(len(dset)):
                res = dset[i].result
                fesc[i] = res["fesc"]
                tot_mass[i] = res["tot_mass"]

            _plot(fesc, tot_mass, label, alpha, ax, c, nbins)


    ax.set_xlabel(r"log$_{10}$(M$_{\mathrm{vir}}$ [M$_{\odot}$])")
    ax.set_ylabel(r"log$_{10}$(f$_{\mathrm{esc}}$ [%])")
    plt.legend()
    plt.show()
    return ax

path = '/its/home/dn99/scratch/bc03' , '/its/home/dn99/scratch/bc03', '/its/home/dn99/scratch/bc03', '/its/home/dn99/scratch/bc03'
iout = 40, 60, 80, 100 
labels = 'z = 12.5', 'z = 9.5', 'z = 7.5', 'z = 6.5'
cols = 'red','blue','orange','purple'    








import numpy as np
import seren3

def load_dbs(halo):
    from seren3.scripts.mpi import time_int_fesc_all_halos, history_mass_flux_all_halos

    fesc_res = time_int_fesc_all_halos.load_halo(halo)
    mass_flux_res = history_mass_flux_all_halos.load_halo(halo)

    return fesc_res, mass_flux_res

def plot_nion_esc(snapshot, ax=None, **kwargs):

#   Plot instantaneous photons escaped as a function of halo mass; adapated from David code which calculated integrated rate

    import matplotlib.pylab as plt
    from seren3.analysis import plots
    from seren3.scripts.mpi import time_int_fesc_all_halos
    from seren3.array import SimArray

    if (ax is None):
        ax = plt.gca()

    color = kwargs.pop("color", "k")
    nbins = kwargs.pop("nbins", 5)
    label = kwargs.pop("label", None)
    ls = kwargs.pop("ls", "-")
    lw = kwargs.pop("lw", 1.)
    
    #edited from here...
    
    fesc_db = time_int_fesc_all_halos.load(snapshot)
    nphotons_escaped = np.zeros(len(fesc_db))
    mvir = np.zeros(len(fesc_db))
    
    for i in range(len(fesc_db)): #collation of instantaneous fesc data vs halo mass data for each halo 
        hid = int(fesc_db[i].idx)
        fesc_res = fesc_db[i].result
        nphotons_escaped[i] = fesc_res["I1"][0]
        mvir[i] = fesc_res["Mvir"]

	# to here. The above collects the instantaneous fesc data attributed 
    #  to each halo. The original function used a cumulative calculation (integral)
 

    log_mvir = np.log10(mvir)
    log_nphotons = np.log10(nphotons_escaped)

    ix = np.where( np.logical_and(np.isfinite(log_nphotons), log_mvir >= 7.5) ) #filters out low mass halos below significant resolution
    log_nphotons = log_nphotons[ix]
    log_mvir = log_mvir[ix]

    bc, mean, std, sterr = plots.fit_scatter(log_mvir, log_nphotons, nbins=nbins, ret_sterr=True)

    ax.scatter(log_mvir, log_nphotons, s=10, color='blue', alpha=0.233)     #alpha determines colour intensity
    e = ax.errorbar(bc, mean, yerr=std, color=color, label=label,\
         fmt="o", markerfacecolor=color, mec='k', capsize=2, capthick=2, elinewidth=2, linestyle=ls, linewidth=lw)

    if (kwargs.pop("legend", False)):
        ax.legend(loc="lower right", frameon=False, prop={"size":16})

    ax.set_xlabel(r'log$_{10}$ M$_{\mathrm{vir}}$ [M$_{\odot}$/h]', fontsize=20)
    ax.set_ylabel(r'log$_{10}$ $\dot{\mathrm{N}}_{\mathrm{ion}}(t)$ f$_{\mathrm{esc}}$ ($t$) [#$s^{-1}$]', fontsize=20)
    plt.show()

snapshot=seren3.load_snapshot('/research/prace/david/bpass/bc03',00100)




def main(path, iout, pickle_path):
    import seren3
    from seren3.core.serensource import DerivedDataset
    from seren3.utils import derived_utils
    from seren3.analysis.escape_fraction import fesc
    from seren3.analysis.parallel import mpi
    from seren3.exceptions import NoParticlesException
    import pickle, os

    mpi.msg("Loading snapshot")
    snap = seren3.load_snapshot(path, iout)
    snap.set_nproc(1)
    halos = None

    star_Nion_d_fn = derived_utils.get_derived_field(snap.s, "Nion_d")
    nIons = snap.info_rt["nIons"]

    halos = snap.halos(finder="ctrees")
    halo_ix = None
    if mpi.host:
        halo_ix = halos.halo_ix(shuffle=True)

    dest = {}
    for i, sto in mpi.piter(halo_ix, storage=dest, print_stats=True):
        h = halos[i]
        star_dset = h.s[["age", "metal", "mass"]].flatten()
        if (len(star_dset["mass"]) > 0):
            dt = h.dt
            ix_young = np.where((star_dset["age"].in_units("Gyr") - dt.in_units("Gyr")) >= 0.)
            if len(ix_young[0] > 0):
                # Work on this halo
                sto.idx = int(h["id"])
                try:
                    dm_dset = h.d["mass"].flatten()
                    gas_dset = h.g["mass"].flatten()

                    star_mass = star_dset["mass"]
                    star_age = star_dset["age"]
                    star_metal = star_dset["metal"]

                    dict_stars = {"age" : star_age, "metal" : star_metal, "mass" : star_mass}
                    dset_stars = DerivedDataset(snap.s, dict_stars)

                    Nion_d_all_groups = snap.array(np.zeros(len(dset_stars["age"])), "s**-1 Msol**-1")

                    for ii in range(nIons):
                        Nion_d_now = star_Nion_d_fn(snap, dset_stars, group=ii+1, dt=0.)
                        Nion_d_all_groups += Nion_d_now

                    tot_mass = dm_dset["mass"].in_units("Msol h**-1").sum() + star_dset["mass"].in_units("Msol h**-1").sum()\
                                     + gas_dset["mass"].in_units("Msol h**-1").sum()

                    h_fesc = fesc(h.subsnap, nside=2**3, filt=False, do_multigroup=True, denoise=True)
                    mpi.msg("%1.2e \t %1.2e" % (h["Mvir"], h_fesc))
                    sto.result = {"fesc" : h_fesc, "tot_mass" : tot_mass, \
                        "Nion_d_now" : Nion_d_all_groups, "star_mass" : star_mass,\
                        "hprops" : h.properties}

                except NoParticlesException as e:
                    mpi.msg(e.message)

    if mpi.host:
        if pickle_path is None:
            pickle_path = "%s/pickle/%s/" % (path, halos.finder.lower())
        # pickle_path = "%s/" % path
        if os.path.isdir(pickle_path) is False:
            os.mkdir(pickle_path)
        unpacked_dest = mpi.unpack(dest)
        fesc_dict = {}
        for i in range(len(unpacked_dest)):
            fesc_dict[int(unpacked_dest[i].idx)] = unpacked_dest[i].result
       #pickle.dump( fesc_dict, open( "%s/fesc_database_no_filt_denoise_%05i.p" % (pickle_path, iout), "wb" ) )
        pickle.dump( fesc_dict, open("%s/time_int_fesc_all_halos_%05i.p" % (pickle_path, iout), "wb") )