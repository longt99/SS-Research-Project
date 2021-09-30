# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
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
from sunpy.net import Fido, attrs as a
from stixpy.net.client import STIXClient
import sunpy.timeseries
from sunpy.timeseries import TimeSeries
from stixpy import timeseries
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import LogNorm
import datetime
import numpy as np

from astropy.time import Time

import radiospectra
import radiospectra.net
from radiospectra.net import sources
from radiospectra.spectrogram2 import Spectrogram

from radiospectra.spectrogram2.spectrogram import GenericSpectrogram

import matplotlib.style
import matplotlib as mpl
plt.rcParams.update(plt.rcParamsDefault)
from cycler import cycler
mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')


# -

def summary_plot(start, end):
    
    '''
    Returns a figure containing data from GOES XRS, STIX, STEREO WAVES (HFR and LFR),
    WIND WAVES (R1 and R2) and I-LOFAR instruments using a start and end date.
    ''' 
    
    # Obtaining data...
    
    query = Fido.search(a.Time(start,end), a.Instrument.goes | a.Instrument.stix |a.Instrument.swaves |a.Instrument.waves |a.Instrument.ilofar)
    print(query)
    
    # GOES
    
    goes_data = Fido.fetch(query['xrs'])
    goes = TimeSeries(goes_data[0]) #
    
    # STIX 
    
    stix_data = Fido.fetch(query['stix'][0])
    stix = TimeSeries(stix_data[0])
    
    # SWAVES (both HFR and LFR)
    
    swaves_data1 = Fido.fetch(query['swaves'][0])
    swaves_spec_l = Spectrogram(sorted(swaves_data1))
    
    lims_sw_l = np.percentile(swaves_spec_l.data, (4, 96))
    
    
    swaves_data2 = Fido.fetch(query['swaves'][2])
    swaves_spec_h = Spectrogram(sorted(swaves_data2))
    
    lims_sw_h = np.percentile(swaves_spec_h.data, (4, 96))
    
    
    # WIND (both Rad 1 and Rad 2)
    
    wind_data = Fido.fetch(query['waves'])
    wind_spec = Spectrogram(sorted(wind_data))
    lims_w = [np.percentile(s.data, (2, 98)) for s in wind_spec]
    
    
    # I-LOFAR
    
    ilofar_data = Fido.fetch(query['ilofarmode357'][0])
    ilofar_spec = Spectrogram(ilofar_data[0])
    lims_i = [np.percentile(s.data, (2, 98)) for s in ilofar_spec]
    
    
    # Creating the summary plot
    
    fig, axes = plt.subplots(7, 1, sharex=True, figsize = (12,18),gridspec_kw={'height_ratios': [1, 1, 1, 1, 2, 1, 1]})
    
    box = dict(facecolor='white', alpha=0.8)
    
    [spec.plot(axes=axes[4], norm=LogNorm(vmin=l[0], vmax=l[1])) for spec, l in zip(ilofar_spec, lims_i)]
    axes[4].text(ilofar_spec[0].times.datetime[150], 30, 'I-LOFAR',color='black', size='10', bbox=box)
    
    wind_spec[0].plot(axes=axes[1], norm=LogNorm(vmin=lims_w[0][0], vmax=lims_w[0][1]))
    axes[1].text(ilofar_spec[0].times.datetime[150], 120, 'WIND R1',color='black', size='10', bbox=box)
    wind_spec[1].plot(axes=axes[2], norm=LogNorm(vmin=lims_w[1][0], vmax=lims_w[1][1]))
    axes[2].text(ilofar_spec[0].times.datetime[150], 2, 'WIND R2',color='black', size='10', bbox=box)
    
 
    swaves_spec_h.plot(axes=axes[3], vmin=lims_sw_h[0], vmax=lims_sw_h[1])
    axes[3].text(ilofar_spec[0].times.datetime[150], 2000, 'STEREO WAVES HFR',color='black', size='10', bbox=box)
    
    swaves_spec_l.plot(axes=axes[0], vmin=lims_sw_l[0], vmax=lims_sw_l[1])
    axes[0].text(ilofar_spec[0].times.datetime[150], 20, 'STEREO WAVES LFR',color='black', size='10', bbox=box)
        
    goes.plot(axes=axes[5], color=['black','red'], markersize=0.1)
    stix.plot(axes=axes[6], markersize=0.1)
    
    # Inverting the axes of the spectrograms so that they are ordered from lower to higher frequency
    
    for ax in axes[:5]:
        ax.invert_yaxis()
    
    # Setting limits of x-axis to the start and end times of the I-LOFAR data
    # and removing the title from each subplot

    for ax in axes:
        ax.set_xlim(ilofar_spec[0].times.datetime[0], ilofar_spec[0].times.datetime[-1])
        ax.set_title('')
    
    axes[5].set_yscale('log')
    axes[5].set_ylabel("Flux (Watts m$^{-2}$)")
    axes[5].legend([r'$1-8\:A$', r'$0.4-5\:A$'],fontsize=10, loc='upper right', frameon=True)
    axes[5].grid(axis='y')
    axes[5].text(ilofar_spec[0].times.datetime[150], 7*10**(-7), 'GOES XRS',color='black', size='10', bbox=box)
    axes[5].set_ylim(10**(-10.5),3*10**(-6))
    
    date_format = mdates.DateFormatter('%H:%M:%S')
    axes[6].xaxis.set_major_formatter(date_format)
    axes[6].legend(fontsize=10, loc='upper right', frameon=True)
    axes[6].set_ylim(5,2*10**(3))
    axes[6].text(ilofar_spec[0].times.datetime[150], 10**(3), 'STIX',color='black', size='10', bbox=box)
    axes[6].set_xlabel('Time (UT)', fontsize='12')
    
    fig.suptitle('Summary Plot for '+ start_day_sample, fontsize='18')
    plt.subplots_adjust(wspace=0, hspace=0)
    fig.subplots_adjust(top=0.95)
    fig.savefig("/Users/thomas/SS Research Project/week 2/sum_plt1.png")
    plt.show()


# +
start_day_sample = "2021/09/07"
end_day_sample = "2021/09/08"

summary_plot(start_day_sample, end_day_sample)
# -


