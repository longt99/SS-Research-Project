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


from radiospectra import net 
from radiospectra.spectrogram2 import Spectrogram
from radiospectra.spectrogram2.sources import WAVESSpectrogram
from radiospectra.spectrogram2.sources import EOVSASpectrogram

import matplotlib.style
import matplotlib as mpl
import matplotlib.ticker as ticker
plt.rcParams.update(plt.rcParamsDefault)
from cycler import cycler
mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['image.cmap'] = 'viridis'
mpl.rcParams['mathtext.fontset'] = 'cm'
from matplotlib import rc
rc('text', usetex=True)
rc('text.latex', preamble=r'\usepackage{cmbright}')

# +
# USING THIS CODE IN ORDER TO DETERMINE THE DRIFT VELOCITY FROM THE RADIO DATA 
# AND TO ALSO GET TIMESERIES PLOTS FROM RADIO DATA

# Fido search for data...

query = Fido.search(a.Time('2021-05-23T10:00', '2021-05-23T12:59'), 
                    a.Instrument.ecallisto)

query

# +
# checking to see which observatories have data...

np.unique(query['CALLISTO']['Observatory']).data

# +
# fetching the data using the observatory location...

index = ((query['CALLISTO']['Observatory'] == 'GLASGOW') )

callisto_files = Fido.fetch(query['CALLISTO'][index], max_conn=1)

# processing files retrieved into spectrogram object...

callisto_specs = Spectrogram(sorted(callisto_files), silence_errors=True)

# finding the limits for clipping...

lims = [np.percentile(s.data, (2, 99)) for s in callisto_specs]

callisto_specs

# +
# plotting the radio data

fig, axes = plt.subplots(figsize=(10,8))

[sp.plot(axes=axes, vmin=l[0], vmax=l[1]) for sp,l in zip(callisto_specs, lims)]
# [sp.plot(axes=axes) for sp in callisto_specs]

# callisto_specs.plot(axes=axes)


# plotting points used to get slope...

axes.plot(datetime(2021, 5, 23, 11, 4, 23, 500000), 80, 'bx', ms='8')
axes.plot(datetime(2021, 5, 23, 11, 4, 25, 500000), 55.5, 'bx', ms='8')

# order not guarnteed so need to set plot range

axes.set_xlim(datetime(2021, 5, 23, 11, 4, 10), datetime(2021, 5, 23, 11, 4, 40))

#axes.set_yscale('log')
axes.invert_yaxis()

axes.set_xlabel('Time (UT)', fontsize='12')
axes.set_ylabel('Frequency (MHz)', fontsize='12')

date_format = mdates.DateFormatter('%H:%M:%S')
axes.xaxis.set_major_formatter(date_format)
axes.xaxis.set_tick_params(rotation=0, labelbottom=True)
plt.xticks(ha='center')

#fig.savefig("/Users/thomas/SS Research Project/week 6/speeds/05-23_2.png")

plt.show()

# +
# checking the time cadence for the radio data...

(callisto_specs[0].times - callisto_specs[0].times[0]).to('s')

# +
# creating objects to be used for plotting radio timeseries

time_s = callisto_specs[4].times.datetime
#time_s

# three frequency bands, summed "row-wise" for timeseries plotting

set1 = callisto_specs[4].data[40:60]
set2 = callisto_specs[4].data[90:110]
set3 = callisto_specs[4].data[140:160]

time_series1 = set1.sum(axis=0)
time_series2 = set2.sum(axis=0)
time_series3 = set3.sum(axis=0)

freq_labels = ['$69.75-73.31\;MHz$', '$60.5-64.0\;MHz$', '$51.06-54.69\;MHz$']

time_s

# +
# creating timeseries plot

fig, axes = plt.subplots(3,1, sharex=True, figsize=(10,6))

spectra = [time_series1,time_series2,time_series3]
colors = ['red','blue','orange']


for ax,spec,c,label in zip(axes,spectra,colors,freq_labels):
    ax.plot(time_s, spec/(np.max(spec)), color=c, label=label)
#     axes.plot(time_s, spec, color=c, label=label)
    ax.set_yscale('linear')
    ax.set_xlim(datetime(2021, 5, 23, 11, 4), datetime(2021, 5, 23, 11, 6))
    ax.legend(frameon=False, fontsize='14', loc='upper right')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True, direction='in', labelsize='12', length=5)


date_format = mdates.DateFormatter('%H:%M:%S')
axes[2].xaxis.set_major_formatter(date_format)
axes[2].xaxis.set_tick_params(rotation=0, labelbottom=True)
axes[2].set_xlabel('Time (UT)', fontsize='12')
#axes[2].set_ylabel('Normalised Intensity (arb. units)', fontsize='12')


fig.text(0.059,0.5, "Normalised Intensity (arb. units)", ha="center", va="center", rotation=90, fontsize='12')
plt.subplots_adjust(wspace=0, hspace=0)

# fig.savefig("/Users/thomas/SS Research Project/week 6/norm2_radio_timeseries.png")
    
plt.show()

# +
# TO-DO #

# add this with STIX science data time series
# -














