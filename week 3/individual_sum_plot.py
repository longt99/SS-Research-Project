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


# +
def get_dates(start, end):
    '''
    Takes a start and end date and returns a list of days where there
    was I-LOFAR observations
    '''
    
    # calling a Fido query from which the list of dates can be generated

    query_dates = Fido.search(a.Time(start, end), a.Instrument.ilofar)

    # appending the dates to a list which will store the start times

    start_list = []

    for i in np.arange(len(query_dates[0])):
        if i % 2 == 0:  # there are two files for each (X and Y polarisation)
            start_list.append(get_dates[0][i][0].datetime) 

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
        
        
        
        qs = []
        
        instruments = [a.Instrument.goes, a.Instrument.stix, a.Instrument.swaves, a.Instrument.waves, a.Instrument.ilofar]
        
        for inst in instruments:
            
            query = Fido.search(a.Time(start, end), inst)

            # if statement which removes any Fido search with no results

            if len(query[0]) == 0:
                continue

            qs.append(query)


        # if statement which ensures query for specific day has results for all instruments specified

        if len(qs) != len(instruments):
            Print('Missing Data on '+ start[:10])
            continue
    
        else:
            queries_list.append(qs)
    
    return queries_list


def multiple_summary_plots(start, end, queries):
    '''
    Creates summary plots using a list of queries 
    '''
    
    for query in queries:
    
        data = fetch_data(query)
        processed_data = process_data(data[0],data[1],data[2],data[3],data[4],data[5],data[6])
        
        summary_plot(start, end, processed_data[0], processed_data[1], processed_data[2], processed_data[3], processed_data[4], processed_data[5], processed_data[6], processed_data[7], processed_data[7], processed_data[9], processed_data[10], processed_data[11])
        plt.show()


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
        
            if len(query[0]) == 0:
                continue
            
            queries.append(query)
            
        except Exception:
            print('There was a prboblem with %s', inst)

    
    return queries


def fetch_data(queries):
    
    '''
    Fetching and assigning data...  need to add in for query in queries if queries is list
    '''
    
    goes_data = Fido.fetch(queries[0][0][0])
    stix_data = Fido.fetch(queries[1][0][0])
    swaves_data1 = Fido.fetch(queries[2][1][0])
    swaves_data2 = Fido.fetch(queries[2][1][1])
    wind_data1 = Fido.fetch(queries[3][0][0])
    wind_data2 = Fido.fetch(queries[3][0][1])
    ilofar_data = Fido.fetch(queries[4][0])
    
    day = str(queries[0][0][0][0].datetime.day)
    month = str(queries[0][0][0][0].datetime.month)
    year = str(queries[0][0][0][0].datetime.year)
    
    date = year+'-'+month+'-'+day
    
    return goes_data, stix_data, swaves_data1, swaves_data2, wind_data1, wind_data2, ilofar_data, date
    
     
def process_data(goes_data, stix_data, swaves_data1, swaves_data2, wind_data1, wind_data2, ilofar_data):
    
    '''
    Converting data from query to Timeseries and Spectrogram objects for plotting
    '''
    
    
    # creating timeseries objects for plotting
    
    goes = TimeSeries(goes_data[0])
    stix = TimeSeries(stix_data[0])
    
    # fetching data for both frequency ranges

    swaves_spec1 = Spectrogram(sorted(swaves_data1))
    swaves_spec2 = Spectrogram(sorted(swaves_data2))
    
    wind_spec1 = Spectrogram(sorted(wind_data1))
    wind_spec2 = Spectrogram(sorted(wind_data2))
    
    # concatenating data for waves rad1 and rad2 and swaves instruments

    sw_data = np.concatenate((swaves_spec1.data, swaves_spec2.data), axis=0)
    sw_freq = np.concatenate((swaves_spec1.frequencies, swaves_spec2.frequencies), axis=0)
    
    wind_data = np.concatenate((wind_spec1.data, wind_spec2.data), axis=0)
    wind_freq = np.concatenate((wind_spec1.frequencies, wind_spec2.frequencies), axis=0)

    # finding the upper and lower limits for normalisation for WAVES/SWAVES
    
    lims_sw = np.percentile( sw_data, (2, 98))
    lims_w = np.percentile(wind_data.data, (2, 98))
    
    
    # I-LOFAR

    ilofar_spec = Spectrogram(ilofar_data[0])
    
    #print(ilofar_data)
    
     # finding the upper and lower limits for normalisation for I-LOFAR
        
    lims_i = [np.percentile(s.data, (2, 98)) for s in ilofar_spec]
    
    return goes, stix, swaves_spec1, wind_spec1, sw_data, sw_freq, wind_data, wind_freq, lims_sw, lims_w, ilofar_spec, lims_i
    
    
def summary_plot(date, start_time, end_time, goes, stix, swaves_spec1, wind_spec1, sw_data, sw_freq, wind_data, wind_freq, lims_sw, lims_w, ilofar_spec, lims_i):  
    
    '''
    Creates a summary plot of the data
    '''
    
    # Creating a time object from start time arg.
    
    start = date + start_time
    end = date + end_time
    
    start_time_obj = datetime.strptime(start, '%Y-%m-%d %H:%M:%S')
    end_time_obj = datetime.strptime(end, '%Y-%m-%d %H:%M:%S')
    
    # Making the summary plot, composed of two spectrograms and two light curves
    
    fig, axes = plt.subplots(4,1, sharex=True, figsize = (10,15), gridspec_kw={'height_ratios': [2, 3, 2, 2]}) #gridspec_kw={'height_ratios': [1, 1, 1, 1, 2, 1, 1]}
    
    # defining the text box features
    
    box = dict(facecolor='white', alpha=1.0, pad=0.1)
    
    # plotting I-LOFAR spectrogram, using the lim_i list to normalise
    
    [spec.plot(axes=axes[1], norm=LogNorm(vmin=l[0], vmax=l[1])) for spec, l in zip(ilofar_spec, lims_i)]
    #axes[1].set_title('I-LOFAR', y=1.0, pad=-14, bbox=box, fontsize='12', loc='left')
    #axes[1].text(wind_spec1.times.datetime[10], 12.3, 'I-LOFAR',color='black', size='14', bbox=box)
    axes[1].set_yscale('log')
    axes[1].set_ylabel('MHz', fontsize='14')
    axes[1].tick_params(which='major', bottom=False, top=False, left=True, right=False, direction='out', labelsize='10', length=5)
    


    # plotting WAVES spectrogram
    
    axes[0].pcolormesh(wind_spec1.times.datetime, wind_freq/1000, wind_data, vmin=lims_w[0], vmax=lims_w[1]) #wind_spec1.times.datetime, wind_freq, vmin=lims_w[0], vmax=lims_w[1]
    axes[0].set_ylabel('MHz', fontsize='14')
    axes[0].set_yscale('log')
    axes[0].tick_params(which='major', bottom=False, top=False, left=True, right=False, direction='out', labelsize='10', length=5)
    

    
    # plotting SWAVES spectrogram

#     axes[0].pcolormesh(swaves_spec1.times.datetime, sw_freq/1000, sw_data, vmin=lims_sw[0], vmax=lims_sw[1])
#     axes[0].set_ylabel('MHz', fontsize='14')
#     axes[0].set_yscale('log')
#     axes[0].tick_params(which='major', bottom=False, top=False, left=True, right=False, direction='out', labelsize='10', length=5)
    
        
    
    # Inverting the axes of the spectrograms so that they are ordered from lower to higher frequency
    
    for ax in axes[:2]:
        ax.invert_yaxis()
    
    # plotting GOES and STIX lightcurves
    
    goes.plot(axes=axes[2], color=['black','red'], markersize=0.1)
    stix.plot(axes=axes[3], markersize=0.1)
    
    # Setting limits of x-axis to the start and end times
    # and removing the title from each subplot

    for ax in axes:
        ax.set_xlim(start_time_obj, end_time_obj)
        ax.set_title('')
    
    axes[2].set_yscale('log')
    axes[2].set_ylabel("Flux (Watts m$^{-2}$)", fontsize='14')
    axes[2].legend([r'$1-8\:\AA$', r'$0.4-5\:\AA$'],fontsize=10, loc='upper right', frameon=True)
    axes[2].grid(axis='y')
    #axes[2].text(wind_spec1.times.datetime[10], 5*10**(-6), 'GOES XRS',color='black', size='14', bbox=box)
    #axes[2].set_title('GOES', y=1.0, pad=-14, bbox=box, loc='left', fontsize='12')
    axes[2].set_ylim(10**(-9),1*10**(-5))
    axes[2].tick_params(which='major', bottom=False, top=False, left=True, right=False, direction='out', labelsize='10', length=5)
    
    
    date_format = mdates.DateFormatter('%H:%M')
    axes[3].xaxis.set_major_formatter(date_format)
    axes[3].xaxis.set_tick_params(rotation=0, labelbottom=True)
    axes[3].legend(fontsize=10, loc='upper right', frameon=True)
    axes[3].set_ylim(5,3*10**(3))
    #axes[3].set_title('STIX', y=1.0, pad=-14, loc='left', fontsize='12')
    #axes[3].text(wind_spec1.times.datetime[10], 1.35*10**(3), 'STIX',color='black', size='14', bbox=box)
    axes[3].set_xlabel('Time (UT) '+ str(start_time_obj)[:10], fontsize='12')
    axes[3].set_ylabel('Counts (s$^{-1}$ keV$^{-1}$)', fontsize='12')
    
    axes[3].tick_params(which='major', bottom=True, top=False, left=True, right=False, direction='out', labelsize='10', length=5)
    axes[3].tick_params(which='minor', bottom=False, top=False, left=True, right=False, direction='out', labelsize='10', length=5)
    #axes[3].xaxis.set_major_locator(plt.MaxNLocator(24))
    
    #axes[0].set_title(' WAVES ', y=1.01, pad=-14, bbox=box, loc='left', fontsize='14')
    axes[0].set_title('WAVES', y=1.0, pad=-14, bbox=box, loc='left', fontsize='14')
    axes[1].set_title(' I-LOFAR ', y=1.01, pad=-14, bbox=box, fontsize='14', loc='left')
    axes[2].set_title(' GOES ', y=1.01, pad=-14, bbox=box, loc='left', fontsize='14')
    axes[3].set_title(' STIX ', y=1.01, pad=-14, bbox=box, loc='left', fontsize='14')
    
    fig.suptitle('Summary Plot for ' + str(start_time_obj)[:10], fontsize='16')
    #plt.xticks(ha='left')
    plt.subplots_adjust(wspace=0, hspace=0)
    fig.subplots_adjust(top=0.95)
    fig.savefig("/Users/thomas/SS Research Project/week 3/plots/waves/sp_"+str(start_time_obj)[:10]+".png")
    #plt.show()

# +
start_day_sample = "2020/11/17 00:00:00"
end_day_sample = "2020/11/17 23:59:59"

test = query_data(start_day_sample, end_day_sample)
# -

test

data = fetch_data(test)
processed = process_data(data[0],data[1],data[2],data[3],data[4],data[5],data[6])

# +
summary_plot(data[7], start_day_sample[10:], end_day_sample[10:], processed[0], processed[1], processed[2], processed[3], processed[4], processed[5], processed[6], processed[7], processed[8], processed[9], processed[10], processed[11])

plt.show()



# -


