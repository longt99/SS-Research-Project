# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.12.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
from radiospectra.spectrogram2 import Spectrogram

from sunpy.net import Fido, attrs as a
from stixpy.net.client import STIXClient
import sunpy.timeseries
from sunpy.timeseries import TimeSeries
from stixpy import timeseries
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from matplotlib.colors import LogNorm

import datetime
from datetime import datetime
import numpy as np

from stixpy import science
from stixpy.science import ScienceData
from stixpy.science import Spectrogram
from sunpy import timeseries as ts

from astropy.time import Time
import astropy.units as u
from astropy.io import fits

import radiospectra
import radiospectra.net
from radiospectra.net import sources


from radiospectra import net #let Fido know about the radio clients
from radiospectra.spectrogram2 import Spectrogram # in the process of updating old spectrogram
from radiospectra.spectrogram2.sources import WAVESSpectrogram
from radiospectra.spectrogram2.sources import EOVSASpectrogram

import matplotlib.style
import matplotlib as mpl
import matplotlib.ticker as ticker
plt.rcParams.update(plt.rcParamsDefault)
from cycler import cycler
mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['image.cmap'] = 'rainbow'
mpl.rcParams['mathtext.fontset'] = 'cm'
from matplotlib import rc
rc('text', usetex=True)
rc('text.latex', preamble=r'\usepackage{cmbright}')

from scipy.optimize import curve_fit


# -

# ## RADIO SPECTROGRAM, TIMESERIES, VELOCITY APPROXIMATION

# +
# USING THIS CODE IN ORDER TO DETERMINE THE DRIFT VELOCITY FROM THE RADIO DATA 
# AND TO ALSO GET TIMESERIES PLOTS FROM RADIO DATA

# Fido search for data...

radio_query = Fido.search(a.Time('2021-06-28T08:00', '2021-06-28T10:00'), 
                    a.Instrument.ecallisto)

# checking to see which observatories have data...

print(np.unique(radio_query['CALLISTO']['Observatory']).data)
# -

str(radio_query[0][0][0])[:]

# +
# fetching the data using the observatory location...

index = ((radio_query['CALLISTO']['Observatory'] == 'GLASGOW') )

callisto_files = Fido.fetch(radio_query['CALLISTO'][index], max_conn=1)

# processing files retrieved into spectrogram object...

callisto_specs = Spectrogram(sorted(callisto_files), silence_errors=True)

# finding the limits for clipping...

lims = [np.percentile(s.data, (1, 97)) for s in callisto_specs]

#callisto_specs

# +
# plotting the radio data

fig, axes = plt.subplots(figsize=(10,8))

[sp.plot(axes=axes, vmin=l[0], vmax=l[1]) for sp,l in zip(callisto_specs, lims)]

# plotting points used to get slope...

axes.plot(datetime(2021, 6, 28, 9, 31, 19, 250000), 60, color='black', marker='X', ms='10')

axes.plot(datetime(2021, 6, 28, 9, 31, 21, 0), 45.5, color='black', marker='X', ms='10')

# order not guarnteed so need to set plot range

axes.set_xlim(datetime(2021, 6, 28, 9, 31, 10), datetime(2021, 6, 28, 9, 31, 40))
axes.set_ylim(45, 65)
#axes.set_ylim(16,70)

#axes.set_yscale('log')
axes.invert_yaxis()

axes.set_xlabel('Time (UT)', fontsize='12')
axes.set_ylabel('Frequency (MHz)', fontsize='12')

date_format = mdates.DateFormatter('%H:%M:%S')
axes.xaxis.set_major_formatter(date_format)
axes.xaxis.set_tick_params(rotation=0, labelbottom=True)
plt.xticks(ha='center')

fig.savefig("/Users/thomas/SS Research Project/week 7/events/"+str(radio_query[0][0][0])[:10]+"_"+str(radio_query[0][0][0])[11:16]+".png")

plt.show()
# -
callisto_specs

# +
# creating objects to be used for plotting radio timeseries

index = -3

time_s = callisto_specs[index].times.datetime
#time_s

# three frequency bands, summed "row-wise" for timeseries plotting

set1 = callisto_specs[index].data[40:60]
set2 = callisto_specs[index].data[90:110]
set3 = callisto_specs[index].data[140:160]

time_series1 = set1.sum(axis=0)
time_series2 = set2.sum(axis=0)
time_series3 = set3.sum(axis=0)

freq_labels = ['$69.75-73.31\;MHz$', '$60.5-64.0\;MHz$', '$51.06-54.69\;MHz$']
# -


print(callisto_specs[index].frequencies[40:60][0], callisto_specs[index].frequencies[40:60][-1])
print(callisto_specs[index].frequencies[90:110][0], callisto_specs[index].frequencies[90:110][-1])
print(callisto_specs[index].frequencies[140:160][0], callisto_specs[index].frequencies[140:160][-1])

# +
# creating timeseries plot

fig, axes = plt.subplots(4,1, sharex=True, figsize=(9,12), gridspec_kw={'height_ratios': [2, 1, 1, 1]})

callisto_specs[index].plot(axes=axes[0], vmin=lims[index][0], vmax=lims[index][1])

axes[0].set_ylabel('', fontsize='12')
axes[0].set_title('')
axes[0].tick_params(which='both', bottom=True, top=True, left=True, right=True, direction='in', labelsize='12', length=5)


spectra = [time_series1,time_series2,time_series3]
colors = ['red','orange','blue']

for ax,spec,c,label in zip(axes[1:],spectra,colors,freq_labels):
    ax.plot(time_s, spec/(np.max(spec)), color=c, label=label)
#     axes.plot(time_s, spec, color=c, label=label)
    ax.set_yscale('linear')
    ax.set_xlim(time_s[100], time_s[1000])
    ax.legend(frameon=False, fontsize='14', loc='upper right')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True, direction='in', labelsize='12', length=5)


date_format = mdates.DateFormatter('%H:%M')
axes[3].xaxis.set_major_formatter(date_format)
axes[3].xaxis.set_tick_params(rotation=0, labelbottom=True)
axes[3].set_xlabel('Time (UT)', fontsize='12')
#axes[2].set_ylabel('Normalised Intensity (arb. units)', fontsize='12')
plt.xticks(ha='center')

fig.text(0.072,0.41, "Normalised Intensity (arb. units)", ha="center", va="center", rotation=90, fontsize='12')
fig.text(0.072,0.75, "Frequency (MHz)", ha="center", va="center", rotation=90, fontsize='12')

plt.subplots_adjust(wspace=0, hspace=0)

fig.savefig("/Users/thomas/SS Research Project/week 7/events/radio_timeseries_"+str(radio_query[0][0][0])[:10]+"_"+str(radio_query[0][0][0])[11:16]+".png")
    
plt.show()
# -

# ## ENERGY SPECTRUM

# +
# querying data

stix_query = Fido.search(a.Time('2021-06-28 00:00:00', '2021-06-28 23:59:59'), a.Instrument.stix,
                    a.stix.DataProduct.sci_xray_cpd)

stix_query[0][6:]

# +
# retrieving file and creating ScienceData object

cpd_file = Fido.fetch(stix_query[0][3])

stix_cpd = ScienceData.from_fits(cpd_file[0])

# Quick look at the spectrogram
stix_cpd.plot_timeseries()
plt.show()

# Check the requested energy channels
stix_cpd.control['energy_bin_mask']

# get indices
_, einds = np.where(stix_cpd.control['energy_bin_mask'].data > 0)

# get peak time (there are better ways)
tind_max = stix_cpd.data['counts'].sum(axis=(1,2,3)).argmax()

# get "lowest" time (there are better ways)
tind_min = stix_cpd.data['counts'].sum(axis=(1,2,3)).argmin()

stix_cpd.energy_masks
# -

einds

# +
# get background data (need to chose this away from any flare)
bg_data = stix_cpd.get_data(time_indices=[tind_min],
                            #time_indices=[[0,4]],# this will sum counts 0-4 time
                            energy_indices=list(range(1,13)), # use all requested channels
#                             energy_indices=einds,
                            detector_indices=[[0,31]],  # sum all detectors
                            pixel_indices=[[0,11]])  # sum all pixels

# get flare data
flare_data = stix_cpd.get_data(time_indices=[tind_max], # this will sum counts 0-4 time
                               energy_indices=list(range(1,13)),
#                                energy_indices=einds,
                               detector_indices=[[0,31]],  # sum all detectors
                               pixel_indices=[[0,11]])  # sum all pixels


# get_data returns a tuple of the counts, count errors, integrated time and energies

# get background subtracted counts
counts = (flare_data[0] - bg_data[0]).reshape(-1)

# calculate errors
count_errs = np.sqrt(flare_data[1]**2 + bg_data[1]**2).reshape(-1)
#count_errs[:] = 0

# defining energy bin widths and points for plotting
e_width = (flare_data[-1]['e_high'] - flare_data[-1]['e_low']).reshape(-1)
e_center = (flare_data[-1]['e_low'] + 0.5 * e_width).reshape(-1)

# fit the count spectra from ~ 6 to 16 keV

def power_law(x, c, alpha):
    return c * x **(-alpha)

popt, pcov = curve_fit(power_law, e_center[2:10].value, counts[2:10].value, 
                       sigma=count_errs[2:10].value)

#plt.plot(e_center[2:14], power_law(e_center[2:14].value, *popt), color='red')

#print(*popt)

# plt.show()

fit_err = np.sqrt(np.diag(pcov))

# +
# plotting the flare count spectrum 

fig, axes = plt.subplots(figsize=(7,5))

plt.errorbar(e_center, counts,  yerr=count_errs, xerr=0.5*e_width, marker='o', ms='3', linestyle='', color='black', capsize=2, label='STIX Science Data')

plt.plot(e_center[2:10], power_law(e_center[2:10].value, *popt), color='red', linestyle='--', linewidth=2, label='Fit to STIX Science Data')

axes.set_ylabel('Counts ($keV^{-1}\;s^{-1}$)', fontsize='14')
axes.set_xlabel('Energy ($keV$)', fontsize='14')

axes.set_xscale('log')
axes.set_yscale('log')

axes.tick_params(which='both', bottom=True, top=True, left=True, right=True, direction='in', labelsize='14', length=5)

# axes.set_xlim(3.75,26)
# axes.set_ylim(0.1,2000)

plt.legend(fontsize=14, loc='lower left', frameon=False)

fig.savefig("/Users/thomas/SS Research Project/week 7/events/spectrum_"+str(stix_query[0][0][0])[:10]+"_"+str(stix_query[0][0][0])[11:16]+".png")

plt.show()
# -

print(popt[1],"+/-",fit_err[1])


