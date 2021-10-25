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

from scipy.optimize import curve_fit

from astropy.time import Time
import astropy.units as u
from astropy.io import fits
from astropy import constants as const

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
mpl.rcParams['image.cmap'] = 'viridis'
mpl.rcParams['mathtext.fontset'] = 'cm'

from matplotlib import rc
rc('text', usetex=True)

rc('text.latex', preamble=r'\usepackage{cmbright}')

# +
# querying data

query = Fido.search(a.Time('2021-05-23 00:00:00', '2021-05-23 23:00:00'), a.Instrument.stix,
                    a.stix.DataProduct.sci_xray_cpd)
# query
# -

query[0][:]


# +
# retrieving file and creating ScienceData object

cpd_file = Fido.fetch(query[0][-7])

stix_cpd = ScienceData.from_fits(cpd_file[0])
# -

cpd_file

# +
# Quick look at the spectrogram
# stix_cpd.plot_spectrogram()
# plt.show()

# Check the requested energy channels
stix_cpd.control['energy_bin_mask']

# get indices
_, einds = np.where(stix_cpd.control['energy_bin_mask'].data > 0)

# get peak time (there are better ways)
tind_max = stix_cpd.data['counts'].sum(axis=(1,2,3)).argmax()

# get "lowest" time (there are better ways)
tind_min = stix_cpd.data['counts'].sum(axis=(1,2,3)).argmin()
# -

stix_cpd.energy_masks

# +

# get background data (need to chose this away from any flare)
bg_data = stix_cpd.get_data(time_indices=[tind_min],
                            #time_indices=[[0,4]],# this will sum counts 0-4 time
                            energy_indices=list(range(1,18)), # use all requested channels
#                             energy_indices=einds,
                            detector_indices=[[0,31]],  # sum all detectors
                            pixel_indices=[[0,11]])  # sum all pixels

# get flare data
flare_data = stix_cpd.get_data(time_indices=[tind_max], # this will sum counts 0-4 time
                               energy_indices=list(range(1,18)),
#                                energy_indices=einds,
                               detector_indices=[[0,31]],  # sum all detectors
                               pixel_indices=[[0,11]])  # sum all pixels


# get_data returns a tuple of the counts, count errors, integrated time and energies

# get background subtracted counts
counts = (flare_data[0] - bg_data[0]).reshape(-1)

# calculate errors
count_errs = np.sqrt(flare_data[1]**2 + bg_data[1]**2).reshape(-1)

# defining energy bin widths and points for plotting
e_width = (flare_data[-1]['e_high'] - flare_data[-1]['e_low']).reshape(-1)
e_center = (flare_data[-1]['e_low'] + 0.5 * e_width).reshape(-1)

# fit the count spectra from ~ 6 to 16 keV

def power_law(x, c, alpha):
    return c * x **(-alpha)

popt, pcov = curve_fit(power_law, e_center[2:-3].value, counts[2:-3].value, 
                       sigma=count_errs[2:-3].value)

plt.plot(e_center[2:-3], power_law(e_center[2:-3].value, *popt), color='red')

print(*popt)

# plt.show()
# -

flare_data

# +
# plotting the flare count spectrum 

fig, axes = plt.subplots(figsize=(7,5))

plt.errorbar(e_center, counts,  yerr=count_errs, xerr=0.5*e_width, marker='o', ms='3', linestyle='', color='black', capsize=2, label='STIX Science Data')

plt.plot(e_center[2:-3], power_law(e_center[2:-3].value, *popt), color='red', linestyle='--', linewidth=2, label='Fit to STIX Science Data')

axes.set_ylabel('Counts ($keV^{-1}\;s^{-1}$)', fontsize='14')
axes.set_xlabel('Energy ($keV$)', fontsize='14')

axes.set_xscale('log')
axes.set_yscale('log')

axes.tick_params(which='both', bottom=True, top=True, left=True, right=True, direction='in', labelsize='14', length=5)

axes.set_xlim(3.75,26)
axes.set_ylim(0.3,2000)

plt.legend(fontsize=14, loc='lower left', frameon=False)

#fig.savefig("/Users/thomas/SS Research Project/week 6/cts_vs_energy.png")

plt.show()
# -


