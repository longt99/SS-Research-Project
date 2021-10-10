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
from datetime import datetime
from datetime import timedelta
from datetime import date
import numpy as np

from astropy.time import Time

import radiospectra
import radiospectra.net
from radiospectra.net import sources
from radiospectra.spectrogram2 import Spectrogram

from radiospectra.spectrogram2.spectrogram import GenericSpectrogram

import matplotlib.ticker as ticker
import matplotlib.style
import matplotlib as mpl
plt.rcParams.update(plt.rcParamsDefault)
from cycler import cycler
mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['image.cmap'] = 'inferno'
mpl.rcParams['mathtext.fontset'] = 'cm'
from matplotlib import rc
rc('text', usetex=True)
rc('text.latex', preamble=r'\usepackage{cmbright}')


# +
def query_data(start, end):
    
    '''
    Scapes GOES XRS, STIX, SWAVES, WAVES and I-LOFAR data
    using a start and end time.
    
    ''' 
    
    # Initiating Fido query...
    
    queries = []
    
    for inst in [a.Instrument.goes, a.Instrument.stix, a.Instrument.swaves, a.Instrument.waves, a.Instrument.ilofar]:
        try:
            query = Fido.search(a.Time(start, end), inst)
            
            # if statement which removes any Fido search with no results
        
#             if len(query[0]) == 0:
#                 continue
            
            queries.append(query)
            
        except Exception:
            print('There was a prboblem with %s', inst)

    
#     query = Fido.search(a.Time(start,end), a.Instrument.goes | a.Instrument.stix |a.Instrument.swaves |a.Instrument.waves |a.Instrument.ilofar)
    
    
    return queries


def fetch_data(start_time, end_time, query):
    
    '''
    Fetches and assigns data to variables from query
    '''
     
    date = query[0][0][0][0]
    
    start = str(date)[:10] + ' ' + start_time
    end = str(date)[:10] + ' ' + end_time
    
    start_time_obj = datetime.strptime(start, '%Y-%m-%d %H:%M:%S')
    end_time_obj = datetime.strptime(end, '%Y-%m-%d %H:%M:%S')
    
    
    
    #fig, axes = plt.subplots(4,1, sharex=True, figsize = (10,14), gridspec_kw={'height_ratios': [2, 3, 2, 2]}) #gridspec_kw={'height_ratios': [1, 1, 1, 1, 2, 1, 1]}
    
    # Defining the text box features
    
    box = dict(facecolor='white', alpha=1.0, pad=0.1, edgecolor = 'black')
    
    
    # If both I-LOFAR and STIX return no results
    
    if ((len(query[4][0]) == 0) & (len(query[1][0]) == 0)):
        
        
        fig, axes = plt.subplots(4,1, sharex=True, figsize = (10,14), gridspec_kw={'height_ratios': [3, 3, 2, 2]})
        
        axes[1].set_yticks([])
        axes[1].tick_params(which='both', bottom=False, top=False, left=False, right=False)
        
        axes[3].set_yticks([])
        axes[3].tick_params(which='both', bottom=False, top=False, left=False, right=False)
    
    
        goes_data = Fido.fetch(query[0][0][0])
        goes = TimeSeries(goes_data)
        goes.plot(axes=axes[2], color=['black','red'], markersize=0.1)
        
        axes[2].set_xlim(start_time_obj, end_time_obj)
        
        axes[2].set_yscale('log')
        axes[2].set_ylabel("Flux ($W\;m^{-2}$)", fontsize='12')
        axes[2].legend([r'$1-8\:\AA$', r'$0.4-5\:\AA$'],fontsize=10, loc='upper right', frameon=True)
        axes[2].grid(axis='y')
        axes[2].set_ylim(10**(-8),10**(-3))
        axes[2].yaxis.set_major_locator(ticker.FixedLocator([10**-7, 10**-6, 10**-5, 10**-4]))
        axes[2].tick_params(which='major', bottom=False, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)
        axes[2].tick_params(which='minor', bottom=False, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)

        flare_class = axes[2].twinx()    
        flare_class.set_yticks(range(0,6))
        flare_class.yaxis.set_major_formatter(ticker.NullFormatter())
        flare_class.yaxis.set_minor_locator(ticker.FixedLocator([0.5,1.5,2.5,3.5,4.5]))

        flare_class.yaxis.set_minor_formatter(ticker.FixedFormatter(['A', 'B', 'C', 'M', 'X']))
        flare_class.tick_params(which='major', bottom=True, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)
        #flare_class.tick_params(which='minor', bottom=False, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)

        
        
        axes[3].set_xlim(start_time_obj, end_time_obj)
        

        
        swaves_data1 = Fido.fetch(query[2][1][0])
        swaves_data2 = Fido.fetch(query[2][1][1])
        
        swaves_spec1 = Spectrogram(swaves_data1)
        swaves_spec2 = Spectrogram(swaves_data2)
        
        sw_data = np.concatenate((swaves_spec1.data, swaves_spec2.data), axis=0)
        sw_freq = np.concatenate((swaves_spec1.frequencies, swaves_spec2.frequencies), axis=0)
        lims_sw = np.percentile(sw_data, (2, 98))
        
        #axes[0].pcolormesh(swaves_spec1.times.datetime, sw_freq/1000, sw_data, vmin=lims_sw[0], vmax=lims_sw[1])
        
        wind_data1 = Fido.fetch(query[3][0][0])
        wind_data2 = Fido.fetch(query[3][0][1])
        wind_spec1 = Spectrogram(wind_data1)
        wind_spec2 = Spectrogram(wind_data2)
        
        wind_data = np.concatenate((wind_spec1.data, wind_spec2.data), axis=0)
        wind_freq = np.concatenate((wind_spec1.frequencies, wind_spec2.frequencies), axis=0)
        lims_w = np.percentile(wind_data.data, (2, 98))
        
        axes[0].pcolormesh(wind_spec1.times.datetime, wind_freq/1000, wind_data, vmin=lims_w[0], vmax=lims_w[1])
        
        axes[0].invert_yaxis()
        axes[0].set_xlim(start_time_obj, end_time_obj)
        axes[0].set_ylabel(' ')
        axes[0].set_yscale('log')
        axes[0].tick_params(which='major', bottom=False, top=False, left=True, right=False, direction='inout', labelsize='12', length=3)
      
        
     
        
        
    # If STIX returns no results
        
    elif len(query[1][0]) == 0:
        
        fig, axes = plt.subplots(4,1, sharex=True, figsize = (10,14), gridspec_kw={'height_ratios': [3, 3, 2, 2]})
        
        axes[3].set_yticks([])
        axes[3].tick_params(which='both', bottom=False, top=False, left=False, right=False, direction='inout', labelsize='12', length=3)
    
        
        goes_data = Fido.fetch(query[0][0][0])
        goes = TimeSeries(goes_data)
        goes.plot(axes=axes[2], color=['black','red'], markersize=0.1)
        
        axes[2].set_xlim(start_time_obj, end_time_obj)
        
        axes[2].set_yscale('log')
        axes[2].set_ylabel("Flux ($W\;m^{-2}$)", fontsize='12')
        axes[2].legend([r'$1-8\:\AA$', r'$0.4-5\:\AA$'],fontsize=10, loc='upper right', frameon=True)
        axes[2].grid(axis='y')
        axes[2].set_ylim(10**(-8),10**(-3))
        axes[2].yaxis.set_major_locator(ticker.FixedLocator([10**-7, 10**-6, 10**-5, 10**-4]))
        axes[2].tick_params(which='major', bottom=False, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)
        axes[2].tick_params(which='minor', bottom=False, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)

        flare_class = axes[2].twinx()    
        flare_class.set_yticks(range(0,6))
        flare_class.yaxis.set_major_formatter(ticker.NullFormatter())
        flare_class.yaxis.set_minor_locator(ticker.FixedLocator([0.5,1.5,2.5,3.5,4.5]))

        flare_class.yaxis.set_minor_formatter(ticker.FixedFormatter(['A', 'B', 'C', 'M', 'X']))
        flare_class.tick_params(which='major', bottom=True, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)
        #flare_class.tick_params(which='minor', bottom=False, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)



        
        axes[3].set_xlim(start_time_obj, end_time_obj)
        
        
        swaves_data1 = Fido.fetch(query[2][1][0])
        swaves_data2 = Fido.fetch(query[2][1][1])

        swaves_spec1 = Spectrogram(swaves_data1)
        swaves_spec2 = Spectrogram(swaves_data2)
        
        sw_data = np.concatenate((swaves_spec1.data, swaves_spec2.data), axis=0)
        sw_freq = np.concatenate((swaves_spec1.frequencies, swaves_spec2.frequencies), axis=0)
        lims_sw = np.percentile(sw_data, (2, 98))
        
        #axes[0].pcolormesh(swaves_spec1.times.datetime, sw_freq/1000, sw_data, vmin=lims_sw[0], vmax=lims_sw[1])
        
      
        wind_data1 = Fido.fetch(query[3][0][0])
        wind_data2 = Fido.fetch(query[3][0][1])
        
        wind_spec1 = Spectrogram(wind_data1)
        wind_spec2 = Spectrogram(wind_data2)
        
        wind_data = np.concatenate((wind_spec1.data, wind_spec2.data), axis=0)
        wind_freq = np.concatenate((wind_spec1.frequencies, wind_spec2.frequencies), axis=0)
        lims_w = np.percentile(wind_data.data, (2, 98))
        
        axes[0].pcolormesh(wind_spec1.times.datetime, wind_freq/1000, wind_data, vmin=lims_w[0], vmax=lims_w[1])
        
        axes[0].invert_yaxis()
        axes[0].set_xlim(start_time_obj, end_time_obj)
        axes[0].set_ylabel(' ')
        axes[0].set_yscale('log')
        axes[0].tick_params(which='major', bottom=False, top=False, left=True, right=False, direction='inout', labelsize='12', length=3)
      
        
        ilofar_data = [Fido.fetch(query[4][0][i]) for i in np.arange(len(query[4][0])) if i % 2 == 0]
        
        ilofar_specs = []
    
        for i, file in enumerate(ilofar_data):
            ilofar_file = Spectrogram(ilofar_data[i])
            ilofar_specs.append(ilofar_file)

        lims_i = []

        for ilofar_spec in ilofar_specs:
            lims = [np.percentile(s.data, (2, 98)) for s in ilofar_spec]
            lims_i.append(lims)
        
        for ilofar_spec, lims in zip(ilofar_specs, lims_i):
            [spec.plot(axes=axes[1], norm=LogNorm(vmin=l[0], vmax=l[1])) for spec, l in zip(ilofar_spec, lims)]
    
        axes[1].set_yscale('log')
        axes[1].invert_yaxis()
        axes[1].set_xlim(start_time_obj, end_time_obj)
        axes[1].set_ylabel(' ', fontsize='12')
        axes[1].tick_params(which='major', bottom=False, top=False, left=True, right=False, direction='inout', labelsize='12', length=3)



        
    # If S/WAVES returns no results
    
    elif len(query[2][1]) == 0:
        
        fig, axes = plt.subplots(4,1, sharex=True, figsize = (10,14), gridspec_kw={'height_ratios': [3, 3, 2, 2]})
        
        axes[0].set_yticks([])
        axes[0].tick_params(which='both', bottom=False, top=False, left=False, right=False, direction='out', labelsize='8', length=5)
    
        
        goes_data = Fido.fetch(query[0][0][0])
        goes = TimeSeries(goes_data)
        goes.plot(axes=axes[2], color=['black','red'], markersize=0.1)
        
        axes[2].set_xlim(start_time_obj, end_time_obj)
        
        axes[2].set_yscale('log')
        axes[2].set_ylabel("Flux ($W\;m^{-2}$)", fontsize='12')
        axes[2].legend([r'$1-8\:\AA$', r'$0.4-5\:\AA$'],fontsize=10, loc='upper right', frameon=True)
        axes[2].grid(axis='y')
        axes[2].set_ylim(10**(-8),10**(-3))
        axes[2].yaxis.set_major_locator(ticker.FixedLocator([10**-7, 10**-6, 10**-5, 10**-4]))
        axes[2].tick_params(which='major', bottom=False, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)
        axes[2].tick_params(which='minor', bottom=False, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)

        flare_class = axes[2].twinx()    
        flare_class.set_yticks(range(0,6))
        flare_class.yaxis.set_major_formatter(ticker.NullFormatter())
        flare_class.yaxis.set_minor_locator(ticker.FixedLocator([0.5,1.5,2.5,3.5,4.5]))

        flare_class.yaxis.set_minor_formatter(ticker.FixedFormatter(['A', 'B', 'C', 'M', 'X']))
        flare_class.tick_params(which='major', bottom=True, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)
        #flare_class.tick_params(which='minor', bottom=False, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)


        
        stix_data = Fido.fetch(query[1][0][0])
        stix = TimeSeries(stix_data)
        stix.plot(axes=axes[3], markersize=0.1)
        
        axes[3].set_xlim(start_time_obj, end_time_obj)
        axes[3].legend(fontsize=10, loc='upper right', frameon=True)
        axes[3].set_ylim(2*10**(0),3.162*10**(4))
        axes[3].set_ylabel('Counts ($s^{-1}\;keV^{-1}$)', fontsize='12')
        
        
        
        wind_data1 = Fido.fetch(query[3][0][0])
        wind_data2 = Fido.fetch(query[3][0][1])

        wind_spec1 = Spectrogram(wind_data1)
        wind_spec2 = Spectrogram(wind_data2)
        
        wind_data = np.concatenate((wind_spec1.data, wind_spec2.data), axis=0)
        wind_freq = np.concatenate((wind_spec1.frequencies, wind_spec2.frequencies), axis=0)
        lims_w = np.percentile(wind_data.data, (2, 98))
        
        axes[0].pcolormesh(wind_spec1.times.datetime, wind_freq/1000, wind_data, vmin=lims_w[0], vmax=lims_w[1])
        
        axes[0].invert_yaxis()
        axes[0].set_xlim(start_time_obj, end_time_obj)
        axes[0].set_ylabel(' ')
        axes[0].set_yscale('log')
        axes[0].tick_params(which='major', bottom=False, top=False, left=True, right=False, direction='inout', labelsize='12', length=3)
      
        
#         ilofar_data = Fido.fetch(query[4][0])
#         ilofar_spec = Spectrogram(ilofar_data[0])
#         lims_i = [np.percentile(s.data, (2, 98)) for s in ilofar_spec]
        
#         [spec.plot(axes=axes[1], norm=LogNorm(vmin=l[0], vmax=l[1])) for spec, l in zip(ilofar_spec, lims_i)]
        
        ilofar_data = [Fido.fetch(query[4][0][i]) for i in np.arange(len(query[4][0])) if i % 2 == 0]
        
        ilofar_specs = []
    
        for i, file in enumerate(ilofar_data):
            ilofar_file = Spectrogram(ilofar_data[i])
            ilofar_specs.append(ilofar_file)

        lims_i = []

        for ilofar_spec in ilofar_specs:
            lims = [np.percentile(s.data, (2, 98)) for s in ilofar_spec]
            lims_i.append(lims)
        
        for ilofar_spec, lims in zip(ilofar_specs, lims_i):
            [spec.plot(axes=axes[1], norm=LogNorm(vmin=l[0], vmax=l[1])) for spec, l in zip(ilofar_spec, lims)]
    
        
        
        
        axes[1].set_yscale('log')
        axes[1].invert_yaxis()
        axes[1].set_xlim(start_time_obj, end_time_obj)
        axes[1].set_ylabel(' ', fontsize='12')
        axes[1].tick_params(which='major', bottom=False, top=False, left=True, right=False, direction='inout', labelsize='12', length=3)



    
    # If WAVES returns no results
    
    elif len(query[3][0]) == 0:
        
        fig, axes = plt.subplots(4,1, sharex=True, figsize = (10,14), gridspec_kw={'height_ratios': [3, 3, 2, 2]})
        
        axes[0].set_yticks([])
        axes[0].tick_params(which='both', bottom=False, top=False, left=False, right=False, direction='out', labelsize='8', length=5)
    
        
        goes_data = Fido.fetch(query[0][0][0])
        goes = TimeSeries(goes_data)
        goes.plot(axes=axes[2], color=['black','red'], markersize=0.1)
        
        axes[2].set_xlim(start_time_obj, end_time_obj)
        
        axes[2].set_yscale('log')
        axes[2].set_ylabel("Flux ($W\;m^{-2}$)", fontsize='12')
        axes[2].legend([r'$1-8\:\AA$', r'$0.4-5\:\AA$'],fontsize=10, loc='upper right', frameon=True)
        axes[2].grid(axis='y')
        axes[2].set_ylim(10**(-8),10**(-3))
        axes[2].yaxis.set_major_locator(ticker.FixedLocator([10**-7, 10**-6, 10**-5, 10**-4]))
        axes[2].tick_params(which='major', bottom=False, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)
        axes[2].tick_params(which='minor', bottom=False, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)

        flare_class = axes[2].twinx()    
        flare_class.set_yticks(range(0,6))
        flare_class.yaxis.set_major_formatter(ticker.NullFormatter())
        flare_class.yaxis.set_minor_locator(ticker.FixedLocator([0.5,1.5,2.5,3.5,4.5]))

        flare_class.yaxis.set_minor_formatter(ticker.FixedFormatter(['A', 'B', 'C', 'M', 'X']))
        flare_class.tick_params(which='major', bottom=True, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)
        #flare_class.tick_params(which='minor', bottom=False, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)


        
        stix_data = Fido.fetch(query[1][0][0])
        stix = TimeSeries(stix_data)
        stix.plot(axes=axes[3], markersize=0.1)
        
        axes[3].set_xlim(start_time_obj, end_time_obj)
        axes[3].legend(fontsize=10, loc='upper right', frameon=True)
        axes[3].set_ylim(2*10**(0),3.162*10**(4))
        axes[3].set_ylabel('Counts ($s^{-1}\;keV^{-1}$)', fontsize='12')
        
        
        # S/WAVES is plotted if available instead
        
        swaves_data1 = Fido.fetch(query[2][1][0])
        swaves_data2 = Fido.fetch(query[2][1][1])

        swaves_spec1 = Spectrogram(swaves_data1)
        swaves_spec2 = Spectrogram(swaves_data2)
        
        sw_data = np.concatenate((swaves_spec1.data, swaves_spec2.data), axis=0)
        sw_freq = np.concatenate((swaves_spec1.frequencies, swaves_spec2.frequencies), axis=0)
        lims_sw = np.percentile(sw_data, (2, 98))
        
        axes[0].pcolormesh(swaves_spec1.times.datetime, sw_freq/1000, sw_data, vmin=lims_sw[0], vmax=lims_sw[1])
        
        axes[0].invert_yaxis()
        axes[0].set_xlim(start_time_obj, end_time_obj)
        axes[0].set_ylabel(' ')
        axes[0].set_yscale('log')
        axes[0].tick_params(which='major', bottom=False, top=False, left=True, right=False, direction='inout', labelsize='12', length=3)
      
        
        
#         ilofar_data = Fido.fetch(query[4][0])
#         ilofar_spec = Spectrogram(ilofar_data[0])
#         lims_i = [np.percentile(s.data, (2, 98)) for s in ilofar_spec]
        
#         [spec.plot(axes=axes[1], norm=LogNorm(vmin=l[0], vmax=l[1])) for spec, l in zip(ilofar_spec, lims_i)]
        
        
        ilofar_data = [Fido.fetch(query[4][0][i]) for i in np.arange(len(query[4][0])) if i % 2 == 0]
        
        ilofar_specs = []
    
        for i, file in enumerate(ilofar_data):
            ilofar_file = Spectrogram(ilofar_data[i])
            ilofar_specs.append(ilofar_file)

        lims_i = []

        for ilofar_spec in ilofar_specs:
            lims = [np.percentile(s.data, (2, 98)) for s in ilofar_spec]
            lims_i.append(lims)
        
        for ilofar_spec, lims in zip(ilofar_specs, lims_i):
            [spec.plot(axes=axes[1], norm=LogNorm(vmin=l[0], vmax=l[1])) for spec, l in zip(ilofar_spec, lims)]
    
        
        
        axes[1].set_yscale('log')
        axes[1].invert_yaxis()
        axes[1].set_xlim(start_time_obj, end_time_obj)
        axes[1].set_ylabel(' ', fontsize='12')
        axes[1].tick_params(which='major', bottom=False, top=False, left=True, right=False, direction='inout', labelsize='12', length=3)



        
    elif len(query[4][0]) == 0:
        
        fig, axes = plt.subplots(4,1, sharex=True, figsize = (10,14), gridspec_kw={'height_ratios': [3, 3, 2, 2]})
        
        axes[1].set_yticks([])
        axes[1].tick_params(which='both', bottom=False, top=False, left=False, right=False, direction='out', labelsize='8', length=5)
    
        
        goes_data = Fido.fetch(query[0][0][0])
        goes = TimeSeries(goes_data)
        goes.plot(axes=axes[2], color=['black','red'], markersize=0.1)
        
        axes[2].set_xlim(start_time_obj, end_time_obj)
        
        axes[2].set_yscale('log')
        axes[2].set_ylabel("Flux ($W\;m^{-2}$)", fontsize='12')
        axes[2].legend([r'$1-8\:\AA$', r'$0.4-5\:\AA$'],fontsize=10, loc='upper right', frameon=True)
        axes[2].grid(axis='y')
        axes[2].set_ylim(10**(-8),10**(-3))
        axes[2].yaxis.set_major_locator(ticker.FixedLocator([10**-7, 10**-6, 10**-5, 10**-4]))
        axes[2].tick_params(which='major', bottom=False, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)
        axes[2].tick_params(which='minor', bottom=False, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)

        flare_class = axes[2].twinx()    
        flare_class.set_yticks(range(0,6))
        flare_class.yaxis.set_major_formatter(ticker.NullFormatter())
        flare_class.yaxis.set_minor_locator(ticker.FixedLocator([0.5,1.5,2.5,3.5,4.5]))

        flare_class.yaxis.set_minor_formatter(ticker.FixedFormatter(['A', 'B', 'C', 'M', 'X']))
        flare_class.tick_params(which='major', bottom=True, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)
        #flare_class.tick_params(which='minor', bottom=False, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)


        
        stix_data = Fido.fetch(query[1][0][0])
        stix = TimeSeries(stix_data)
        stix.plot(axes=axes[3], markersize=0.1)
        
        axes[3].set_xlim(start_time_obj, end_time_obj)
        axes[3].legend(fontsize=10, loc='upper right', frameon=True)
        axes[3].set_ylim(2*10**(0),3.162*10**(4))
        axes[3].set_ylabel('Counts ($s^{-1}\;keV^{-1}$)', fontsize='12')
        
        
        
        swaves_data1 = Fido.fetch(query[2][1][0])
        swaves_data2 = Fido.fetch(query[2][1][1])

        swaves_spec1 = Spectrogram(swaves_data1)
        swaves_spec2 = Spectrogram(swaves_data2)
        
        sw_data = np.concatenate((swaves_spec1.data, swaves_spec2.data), axis=0)
        sw_freq = np.concatenate((swaves_spec1.frequencies, swaves_spec2.frequencies), axis=0)
        lims_sw = np.percentile(sw_data, (2, 98))
        
        #axes[0].pcolormesh(swaves_spec1.times.datetime, sw_freq/1000, sw_data, vmin=lims_sw[0], vmax=lims_sw[1])

        
        wind_data1 = Fido.fetch(query[3][0][0])
        wind_data2 = Fido.fetch(query[3][0][1])
        
        wind_spec1 = Spectrogram(wind_data1)
        wind_spec2 = Spectrogram(wind_data2)
        
        wind_data = np.concatenate((wind_spec1.data, wind_spec2.data), axis=0)
        wind_freq = np.concatenate((wind_spec1.frequencies, wind_spec2.frequencies), axis=0)
        lims_w = np.percentile(wind_data.data, (2, 98))
        
        axes[0].pcolormesh(wind_spec1.times.datetime, wind_freq/1000, wind_data, vmin=lims_w[0], vmax=lims_w[1])
        
        axes[0].invert_yaxis()
        axes[0].set_xlim(start_time_obj, end_time_obj)
        axes[0].set_ylabel(' ')
        axes[0].set_yscale('log')
        axes[0].tick_params(which='major', bottom=False, top=False, left=True, right=False, direction='inout', labelsize='12', length=3)
      


    # If GOES returns no results
    
    elif len(query[0][0]) == 0:
        
        
        
        fig, axes = plt.subplots(4,1, sharex=True, figsize = (10,14), gridspec_kw={'height_ratios': [3, 3, 2, 2]}) #gridspec_kw={'height_ratios': [1, 1, 1, 1, 2, 1, 1]}
        
        
        #axes[2].text(0.4, 0.4, 'MISSING DATA',color='black', size='16')
        axes[2].tick_params(which='both', bottom=False, top=False, left=False, right=False, direction='out', labelsize='8', length=5)
        axes[2].set_yticks([])
         
        
        stix_data = Fido.fetch(query[1][0][0])
        stix = TimeSeries(stix_data)
        stix.plot(axes=axes[3], markersize=0.1)
        
        axes[3].set_xlim(start_time_obj, end_time_obj)
        axes[3].legend(fontsize=10, loc='upper right', frameon=True)
        axes[3].set_ylim(2*10**(0),3.162*10**(4))
        axes[3].set_ylabel('Counts ($s^{-1}\;keV^{-1}$)', fontsize='12')
        
        
        
        swaves_data1 = Fido.fetch(query[2][1][0])
        swaves_data2 = Fido.fetch(query[2][1][1])

        swaves_spec1 = Spectrogram(swaves_data1)
        swaves_spec2 = Spectrogram(swaves_data2)
        
        sw_data = np.concatenate((swaves_spec1.data, swaves_spec2.data), axis=0)
        sw_freq = np.concatenate((swaves_spec1.frequencies, swaves_spec2.frequencies), axis=0)
        lims_sw = np.percentile(sw_data, (2, 98))
        
        #axes[0].pcolormesh(swaves_spec1.times.datetime, sw_freq/1000, sw_data, vmin=lims_sw[0], vmax=lims_sw[1])
        
        
        wind_data1 = Fido.fetch(query[3][0][0])
        wind_data2 = Fido.fetch(query[3][0][1])

        wind_spec1 = Spectrogram(wind_data1)
        wind_spec2 = Spectrogram(wind_data2)
        
        wind_data = np.concatenate((wind_spec1.data, wind_spec2.data), axis=0)
        wind_freq = np.concatenate((wind_spec1.frequencies, wind_spec2.frequencies), axis=0)
        lims_w = np.percentile(wind_data.data, (2, 98))
        
        axes[0].pcolormesh(wind_spec1.times.datetime, wind_freq/1000, wind_data, vmin=lims_w[0], vmax=lims_w[1])
        
        axes[0].invert_yaxis()
        axes[0].set_xlim(start_time_obj, end_time_obj)
        axes[0].set_ylabel(' ')
        axes[0].set_yscale('log')
        axes[0].tick_params(which='major', bottom=False, top=False, left=True, right=False, direction='inout', labelsize='12', length=3)
      
        
#         ilofar_data = Fido.fetch(query[4][0][0])
#         ilofar_spec = Spectrogram(ilofar_data)
#         lims_i = [np.percentile(s.data, (2, 98)) for s in ilofar_spec]
        
#         [spec.plot(axes=axes[1], norm=LogNorm(vmin=l[0], vmax=l[1])) for spec, l in zip(ilofar_spec, lims_i)]
        
    
        ilofar_data = [Fido.fetch(query[4][0][i]) for i in np.arange(len(query[4][0])) if i % 2 == 0]
        
        ilofar_specs = []
    
        for i, file in enumerate(ilofar_data):
            ilofar_file = Spectrogram(ilofar_data[i])
            ilofar_specs.append(ilofar_file)

        lims_i = []

        for ilofar_spec in ilofar_specs:
            lims = [np.percentile(s.data, (2, 98)) for s in ilofar_spec]
            lims_i.append(lims)
        
        for ilofar_spec, lims in zip(ilofar_specs, lims_i):
            [spec.plot(axes=axes[1], norm=LogNorm(vmin=l[0], vmax=l[1])) for spec, l in zip(ilofar_spec, lims)]
    
        
        axes[1].set_yscale('log')
        axes[1].invert_yaxis()
        axes[1].set_xlim(start_time_obj, end_time_obj)
        axes[1].set_ylabel(' ', fontsize='12')
        axes[1].tick_params(which='major', bottom=False, top=False, left=True, right=False, direction='inout', labelsize='12', length=3)


        
    
    
    else:
        
        fig, axes = plt.subplots(4,1, sharex=True, figsize = (10,14), gridspec_kw={'height_ratios': [3, 3, 2, 2]})
        
        goes_data = Fido.fetch(query[0][0][0])
        goes = TimeSeries(goes_data)
        goes.plot(axes=axes[2], color=['black','red'], markersize=0.1)
        
        axes[2].set_xlim(start_time_obj, end_time_obj)
        
        axes[2].set_yscale('log')
        axes[2].set_ylabel("Flux ($W\;m^{-2}$)", fontsize='12')
        axes[2].legend([r'$1-8\:\AA$', r'$0.4-5\:\AA$'],fontsize=10, loc='upper right', frameon=True)
        axes[2].grid(axis='y')
        axes[2].set_ylim(10**(-8),10**(-3))
        axes[2].yaxis.set_major_locator(ticker.FixedLocator([10**-7, 10**-6, 10**-5, 10**-4]))
        axes[2].tick_params(which='major', bottom=False, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)
        axes[2].tick_params(which='minor', bottom=False, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)

        flare_class = axes[2].twinx()    
        flare_class.set_yticks(range(0,6))
        flare_class.yaxis.set_major_formatter(ticker.NullFormatter())
        flare_class.yaxis.set_minor_locator(ticker.FixedLocator([0.5,1.5,2.5,3.5,4.5]))

        flare_class.yaxis.set_minor_formatter(ticker.FixedFormatter(['A', 'B', 'C', 'M', 'X']))
        flare_class.tick_params(which='major', bottom=True, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)
        #flare_class.tick_params(which='minor', bottom=False, top=False, left=True, right=True, direction='inout', labelsize='12', length=3)


        
        stix_data = Fido.fetch(query[1][0][0])
        stix = TimeSeries(stix_data)
        stix.plot(axes=axes[3], markersize=0.1)
        
        axes[3].set_xlim(start_time_obj, end_time_obj)
        axes[3].legend(fontsize=10, loc='upper right', frameon=True)
        axes[3].set_ylim(2*10**(0),3.162*10**(4))
        axes[3].set_ylabel('Counts ($s^{-1}\;keV^{-1}$)', fontsize='12')

        
        
        swaves_data1 = Fido.fetch(query[2][1][0])
        swaves_data2 = Fido.fetch(query[2][1][1])

        swaves_spec1 = Spectrogram(swaves_data1)
        swaves_spec2 = Spectrogram(swaves_data2)
        
        sw_data = np.concatenate((swaves_spec1.data, swaves_spec2.data), axis=0)
        sw_freq = np.concatenate((swaves_spec1.frequencies, swaves_spec2.frequencies), axis=0)
        lims_sw = np.percentile(sw_data, (2, 98))
        
        #axes[0].pcolormesh(swaves_spec1.times.datetime, sw_freq/1000, sw_data, vmin=lims_sw[0], vmax=lims_sw[1])  
        
        wind_data1 = Fido.fetch(query[3][0][0])
        wind_data2 = Fido.fetch(query[3][0][1])

        wind_spec1 = Spectrogram(wind_data1)
        wind_spec2 = Spectrogram(wind_data2)
        
        wind_data = np.concatenate((wind_spec1.data, wind_spec2.data), axis=0)
        wind_freq = np.concatenate((wind_spec1.frequencies, wind_spec2.frequencies), axis=0)
        lims_w = np.percentile(wind_data.data, (2, 98))
        
        axes[0].pcolormesh(wind_spec1.times.datetime, wind_freq/1000, wind_data, vmin=lims_w[0], vmax=lims_w[1])
        
        axes[0].invert_yaxis()
        axes[0].set_ylabel(' ', fontsize='12')
        axes[0].set_yscale('log')
        axes[0].set_xlim(start_time_obj, end_time_obj)
        axes[0].tick_params(which='major', bottom=False, top=False, left=True, right=False, direction='inout', labelsize='12', length=3)
        
        
#         ilofar_data = Fido.fetch(query[4][0][0])
#         ilofar_spec = Spectrogram(ilofar_data)
#         lims_i = [np.percentile(s.data, (2, 98)) for s in ilofar_spec]
        
#         [spec.plot(axes=axes[1], norm=LogNorm(vmin=l[0], vmax=l[1])) for spec, l in zip(ilofar_spec, lims_i)]
        
        ilofar_data = [Fido.fetch(query[4][0][i]) for i in np.arange(len(query[4][0])) if i % 2 == 0]
        
        ilofar_specs = []
    
        for i, file in enumerate(sorted(ilofar_data)):
            ilofar_file = Spectrogram(ilofar_data[i])
            ilofar_specs.append(ilofar_file)

        lims_i = []

        for ilofar_spec in ilofar_specs:
            lims = [np.percentile(s.data, (2, 90)) for s in ilofar_spec]
            lims_i.append(lims)
        
        for ilofar_spec, lims in zip(ilofar_specs, lims_i):
            [spec.plot(axes=axes[1], norm=LogNorm(vmin=l[0], vmax=l[1])) for spec, l in zip(ilofar_spec, lims)]
    
        
        
        axes[1].set_yscale('log')
        axes[1].invert_yaxis()
        axes[1].set_xlim(start_time_obj, end_time_obj)
        axes[1].set_ylabel(' ', fontsize='12')
        axes[1].tick_params(which='major', bottom=False, top=False, left=True, right=False, direction='inout', labelsize='12', length=3)

    
    for ax in axes:
#         ax.set_xlim(start_time_obj, end_time_obj)
        ax.set_xlim(ilofar_specs[0][0].times.datetime[0], ilofar_specs[-1][-1].times.datetime[-1])
        ax.set_title('')
#         date_format = mdates.DateFormatter('%H:%M')
#         ax.xaxis.set_major_formatter(date_format)
        
    
    
    date_format = mdates.DateFormatter('%H:%M')
    axes[3].xaxis.set_major_formatter(date_format)
    axes[3].xaxis.set_tick_params(rotation=0, labelbottom=True)
    axes[3].set_xlabel('Time (UT) '+ str(start_time_obj)[:10], fontsize='12')
    
    axes[3].tick_params(which='major', bottom=True, top=False, left=True, right=False, direction='inout', labelsize='12', length=3)
    axes[3].tick_params(which='minor', bottom=False, top=False, left=True, right=False, direction='inout', labelsize='12', length=3)
    plt.xticks(ha='center')
    
    axes[0].set_title('  (a)', y=1.0, pad=-9, loc='left', bbox=box, fontsize='12', color='black')
    axes[1].set_title('  (b)', y=1.0, pad=-9, loc='left', bbox=box, fontsize='12', color='black')
    axes[2].set_title('  (c)', y=1.0, pad=-9, loc='left', fontsize='12', color='black')
    axes[3].set_title('  (d)', y=1.0, pad=-9, loc='left', fontsize='12', color='black')
    
    #fig.suptitle('Summary Plot for ' + str(start_time_obj)[:10], fontsize='16')
    plt.subplots_adjust(wspace=0, hspace=0)
    fig.text(0.059,0.75, "Frequency ($MHz$)", ha="center", va="center", rotation=90, fontsize='12')
    fig.align_ylabels(axes)
    fig.subplots_adjust(top=0.95)
    fig.savefig("/Users/thomas/SS Research Project/week 4/ilofar_obs/"+str(start_time_obj)[:10]+".png")
    #plt.show()
    
    
    


# +
def get_dates(start, end):
    '''
    Takes a start and end date and returns a list of days where there
    were I-LOFAR observations
    '''
    
    # calling a Fido query from which the list of dates can be generated

    query_dates = Fido.search(a.Time(start, end), a.Instrument.ilofar)

    # appending the dates to a list which will store the start times

    start_list = []

    for i in np.arange(len(query_dates[0])):
        if i % 2 == 0:  # there are two files for each (X and Y polarisation)
            start_list.append(query_dates[0][i][0].datetime) 

    start_list = list(set(sorted(start_list)))
    
    return start_list


def multiple_dates(dates):
    '''
    Takes a list of dates and returns a list of Fido queries
    '''
    
    queries_list = []
    
    for date in dates:
        
        # creating start and end strings
        
        year = str(date.year)
        month = str(date.month)
        day = str(date.day)

        start = year + '/' + month + '/' + day + ' 00:00:00'
        end = year + '/' + month + '/' + day + ' 23:59:59'


        # calling a Fido search for each instrument
        
        result = query_data(start, end)

        queries_list.append(result)
    
    return queries_list


def multiple_summary_plots(start_time, end_time, queries):
    '''
    Creates summary plots using a list of queries 
    '''
    
    for query in queries:
    
        data = fetch_data(start_time, end_time, query)
#         processed_data = process_data(data[0],data[1],data[2],data[3],data[4],data[5],data[6])
        
#         summary_plot(data[7], start_time, end_time, processed_data[0], processed_data[1], processed_data[2], processed_data[3], processed_data[4], processed_data[5], processed_data[6], processed_data[7], processed_data[7], processed_data[9], processed_data[10], processed_data[11])
        plt.show()
# +
start_date = date(2021, 6, 22)
number_of_days = 2

date_list = []
for day in range(number_of_days):
    a_date = (start_date + timedelta(days = day)) #.isoformat()
    date_list.append(a_date)


#all = multiple_dates()
# -


queries = multiple_dates(date_list[:])

queries[109]

multiple_summary_plots('00:00:00', '23:59:59', queries)
plt.show()

queries[0][2][1][0]







# +
start = '2021/07/01 00:00:00'
end = '2021/07/01 23:59:00'



query = query_data(start, end)
# -









# +
# Example:

start_date = '2021/04/30'
end_date = '2021/05/31'

start_time = '00:00:00'
end_time = '23:59:00'

date_list = get_dates(start_date, end_date)
# data = multiple_dates(date_list)

#summary_plots = multiple_summary_plots(start_time, end_time, data)
# -
date_list

data = multiple_dates(sorted(date_list))
data

summary_plots = multiple_summary_plots(start_time, end_time, data[10:])













query = query_data('2021/09/28 00:00:00', '2021/09/28 23:59:00')


data = fetch_data(query)


