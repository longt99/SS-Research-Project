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
from astropy import constants as const

import radiospectra
import radiospectra.net
from radiospectra.net import sources

# from radiospectra.spectrum import Spectrum
# from radiospectra.spectrogram2 import Spectrogram
# from radiospectra.spectrogram2.spectrogram import GenericSpectrogram
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


from shapely.geometry import LineString

from scipy.optimize import curve_fit
# -

M_sun = const.M_sun
G = const.G
k = const.k_B
R_sun = const.R_sun
g_sun = 274 * u.meter / (u.second)**2
T = 1e6 * u.Kelvin
mu = 0.6
m_p = 1.6726219e-27 * u.kg


# +
def freq(n):
    
    '''
    Relationship between plasma frequency and electron density
    '''
    
    f = (8.93e-3)*(n)**(1/2)
    
    return f

def density(f):
    
    '''
    Relationship between plasma frequency and electron density
    '''
    
    n = (f/(8.93e-3))**2
    
    return n


# NEWKIRK MODEL

def newkirk(r):
    
    '''
    Newkirk density model
    '''
    
    n = ((4.2e4 * u.cm**-3)*(10**(4.32/r)))
    
    return n


def rad(f):
    
    '''
    Distance using Newkirk model
    '''
    
    r = 2.16/(np.log10(f/1.83))
    
    return r


def err_r(f1, f2):
    
    '''
    Getting error in Newkirk distance
    '''
    
    err = np.sqrt(((2.16*0.125)/(f1*(np.log10(f1/1.83))**2))**2 + ((2.16*0.125)/(f2*(np.log10(f2/1.83))**2))**2)
    #err = (2.16*0.125)/(f1*(np.log10(f1/1.83))**2)
    return err


def err_v(v, dr, r, t):
    
    '''
    Getting error in the Newkirk velocity
    '''
    
    err = v*np.sqrt((dr/r)**2 + (0.35355/t)**2)
    return err


# MANN MODEL

def mann(r):
    
    '''
    Mann density model
    '''
    
    n = (5.14e9 * u.cm**-3)*np.exp(((mu*m_p*G*M_sun)/(k*T*R_sun))*(1/r - 1))  # n_ss_0 = 5.14e9 * u.cm**-3
    
    return n



def r_m(f):
    
    '''
    Distance using Mann Model
    '''
    
    r = ((1)/(1+((2*k*T*R_sun)/(mu*m_p*G*M_sun))*np.log(f/644)))
    
    return r


def err_rm(f1, f2):
    
    '''
    Getting error in Mann distance
    '''
    e1 = ((2*k*T*R_sun*0.125)/(mu*m_p*G*M_sun))/(f1*((1+((2*k*T*R_sun)/(mu*m_p*G*M_sun))*np.log(f1/644)))**2)
    e2 = ((2*k*T*R_sun*0.125)/(mu*m_p*G*M_sun))/(f2*((1+((2*k*T*R_sun)/(mu*m_p*G*M_sun))*np.log(f2/644)))**2)
    
    err = np.sqrt((e1)**2 + (e2)**2)
    #err = (2.16*0.125)/(f1*(np.log10(f1/1.83))**2)
    return err


def err_vm(v, dr, r, t):
    
    '''
    Getting error in the Mann velocity
    '''
    
    err = v*np.sqrt((dr/r)**2 + (0.35355/t)**2)
    return err

# SAITO MODEL

def saito(r):
    
    '''
    Saito density model
    '''
    
    n = (1.36e6)*r**-2.14 + (1.68e8)*r**-6.13
    
    return n


def interpolate(x, f):
    
    '''
    Method for interpolating the heliocentric distance
    '''
    n = density(f)
    
    line = SG.LineString(list(zip(x, saito(x))))
    yline = SG.LineString([(min(x), n), (max(x), n)])
    coords = np.array(line.intersection(yline))
    
    return coords[0]


# def r_saito(r1, r2):
    
#     '''
#     Distance using Saito Model
#     '''
    
#     dr = r1 - r2
    
#     return dr

# def err_rst(dr):
    
#     '''
#     Getting error in Saito distance
#     '''
#     e1 = ((2*k*T*R_sun*0.125)/(mu*m_p*G*M_sun))/(f1*((1+((2*k*T*R_sun)/(mu*m_p*G*M_sun))*np.log(f1/644)))**2)
#     e2 = ((2*k*T*R_sun*0.125)/(mu*m_p*G*M_sun))/(f2*((1+((2*k*T*R_sun)/(mu*m_p*G*M_sun))*np.log(f2/644)))**2)
    
#     err = np.sqrt((e1)**2 + (e2)**2)
#     #err = (2.16*0.125)/(f1*(np.log10(f1/1.83))**2)
#     return err


def err_vs(v, v_s, dv_s):
    
    '''
    Getting error in the Saito velocity
    '''
    
    err = v*(dv/v_s)
    
    return err


# -

interpolate(x, 20)

# +
# PLOTTING THE DENSITY MODELS

# x = np.logspace(0,2.2,1000)
x = np.linspace(1,2.5,1000)

fig, axes = plt.subplots(figsize=(8,5))

axes.plot(x, mann(x), color='b', label='Mann')

axes.plot(x, newkirk(x), color='r', label='Newkirk')

axes.plot(x, saito(x), color='magenta', label='Saito')

# frequency = axes.twinx()
# frequency.plot(x, freq(mann(x)), color='white')
# frequency.set_yscale('log')

# axes.axhline(4.938e6, color='g')
# axes.axvline(1.56, color='g')

axes.set_xlabel('Distance ($R_{\odot}$)', fontsize='14')
axes.set_ylabel('Electron Density ($cm^{-3}$)', fontsize='14')

axes.set_xlim(1, 2.5)
axes.set_ylim(7*10**4, 10**10)

# axes.set_xscale('log')
axes.set_yscale('log')

axes.tick_params(which='both', bottom=True, top=True, left=True, right=True, direction='in', labelsize='14', length=5)

plt.legend(fontsize=14, loc='upper right', frameon=False)
#plt.grid()

#fig.savefig("/Users/thomas/SS Research Project/week 7/models.png")

plt.show()
# -

import shapely.geometry as SG

freq(4.93800000e+06)

# +
# # PLOTTING DISTANCE VS FREQUENCY FOR SAITO MODEL

# x = np.linspace(1,2.5,1000)

# fig, axes = plt.subplots(figsize=(8,5))

# axes.plot(x, freq(saito(x)), color='b', label='Saito')



# line = SG.LineString(list(zip(x, freq(saito(x)))))
# freq_0 = 20
# yline = SG.LineString([(min(x), freq_0), (max(x), freq_0)])
# coords = np.array(line.intersection(yline))
# print(coords[:])


# # frequency = axes.twinx()
# # frequency.plot(x, freq(mann(x)), color='white')
# # frequency.set_yscale('log')

# #axes.axhline(4.938e6, color='g')
# # axes.axvline(1.56, color='g')

# axes.set_xlabel('Distance ($R_{\odot}$)', fontsize='14')
# axes.set_ylabel('Frequency ($MHz$)', fontsize='14')

# axes.set_xlim(1, 2.5)
# #axes.set_ylim(7*10**4, 10**10)

# # axes.set_xscale('log')
# axes.set_yscale('log')

# axes.tick_params(which='both', bottom=True, top=True, left=True, right=True, direction='in', labelsize='14', length=5)

# plt.legend(fontsize=14, loc='upper right', frameon=False)
# #plt.grid()

# #fig.savefig("/Users/thomas/SS Research Project/week 7/models.png")

# plt.show()


# +
# Data

start = [80, 80, 80, 60, 60, 60, 60, 63, 58, 80, 75, 80, 80, 80, 80]
end = [45, 45, 55, 45, 25, 20, 25, 38, 20, 45, 30, 45, 45, 45, 30]
time = [1.25, 1.25, 1, 1.75, 4.5, 6.5, 2.5, 2, 3.5, 1.25, 3, 1.75, 2, 1.75, 3.75]
spec_ind = [4.785, 5.097, 5.6115, 6.209, 5.626, 6.148, 5.119, 5.634, 5.100, 5.110, 5.601, 5.944, 5.57, 5.282, 5.576]
sp_ind_err = [0.172 ,0.0642, 0.1826, 0.274, 0.152, 0.242, 0.0712, 0.203, 0.228, 0.153, 0.201, 0.182, 0.203, 0.209, 0.198]


# +
# Newkirk calculations

r_s = []
v_s = []
dr_s = []
dv_s = []

for s,e,t in zip(start, end, time):
    
    r = rad(e) - rad(s)
    r_s.append(r)

    v = (r*const.R_sun)/(const.c*t*u.second)
    v_s.append(v)

    dr = err_r(e, s)
    dr_s.append(dr)

    dv = err_v(v, dr, r, t)
    dv_s.append(dv)


print('Newkirk','\n')

for v,dv in zip(v_s, dv_s):
    print(v,'+/-',dv)

# for v,dv in zip(v_s, dv_s):
#     print(dv/v)

# +
# Mann Calculations

r_ms = []
v_ms = []
dr_ms = []
dv_ms = []

for s,e,t in zip(start, end, time):
    
    r_mann = r_m(e).value - r_m(s).value
    r_ms.append(r_mann)
    
    v_mann = (r_mann*const.R_sun)/(const.c*t*u.second)
    v_ms.append(v_mann)
    
    dr_mann = err_rm(e, s)
    dr_ms.append(dr_mann)

    dv_mann = err_vm(v_mann, dr_mann, r_mann, t)
    dv_ms.append(dv_mann)
    

print('Mann','\n')
    
for v, dv in zip(v_ms,dv_ms):
    
    print(v, '+/-', dv)

    
# for v,dv in zip(v_ms, dv_ms):
#     print(dv/v)

# +
# Saito Calculations

x = np.linspace(1,2.5,1000)

r_ss = []
v_ss = []
dr_ss = []
dv_ss = []

for s,e,t,v,dv in zip(start, end, time, v_s, dv_s):
    
    r_saito = interpolate(x, e) - interpolate(x, s)
    r_ss.append(r_saito)
    
    v_saito = (r_saito*const.R_sun)/(const.c*t*u.second)
    v_ss.append(v_saito)
    
#     dr_saito = err_rm(e, s)
#     dr_ms.append(dr_mann)

    dv_saito = err_vs(v_saito, v, dv)
    dv_ss.append(dv_saito)
    

print('Saito','\n')
    
for v, dv in zip(v_ss,dv_ss):
    
    print(v, '+/-', dv)
    
print('\n', 'Newkirk', '\n')

for v,dv in zip(v_s, dv_s):
    
    print(v,'+/-',dv)


# +
# Fitting the data to following equation 

def fit(x, m, c):
    
    return m*x + c
    
# def fit(x, m, c):    
    
#     return m * x **(-c)

# Newkirk model Fit

popt_n, pcov_n = curve_fit(fit, v_s, spec_ind, 
                       sigma=sp_ind_err)

fit_err_n = np.sqrt(np.diag(pcov_n))


# Mann model Fit

popt_m, pcov_m = curve_fit(fit, v_ms, spec_ind, 
                       sigma=sp_ind_err)

fit_err_m = np.sqrt(np.diag(pcov_m))


# Saito model Fit

popt_s, pcov_s = curve_fit(fit, v_ss, spec_ind, 
                       sigma=sp_ind_err)

fit_err_s = np.sqrt(np.diag(pcov_s))

# -

print(popt_n, '\n', popt_m, '\n', popt_s)

# +
# Plotting spectral index vs. drift velocity fits to data for each model 

x_val = np.linspace(0.1,0.6,200)

fig, axes = plt.subplots(figsize=(8,5))

# with error bars...

# plt.errorbar(v_s, spec_ind,  yerr=sp_ind_err, xerr=dv_s, linestyle='', color='black', capsize=2, zorder=0)
# plt.errorbar(v_ms, spec_ind,  yerr=sp_ind_err, xerr=dv_ms, linestyle='', color='black', capsize=2, zorder=0)

plt.scatter(v_s, spec_ind, color='r', label = 'Newkirk Model', zorder=5)
plt.scatter(v_ms, spec_ind, color='b', label = 'Mann Model', zorder=5)
plt.scatter(v_ss, spec_ind, color='g', label = 'Saito Model', zorder=5)


plt.plot(x_val, fit(x_val, *popt_n), color='red', linestyle='--', linewidth=2, label='Fit to Newkirk', zorder=10)
plt.plot(x_val, fit(x_val, *popt_m), color='blue', linestyle='--', linewidth=2, label='Fit to Mann', zorder=10)
plt.plot(x_val, fit(x_val, *popt_s), color='green', linestyle='--', linewidth=2, label='Fit to Saito', zorder=10)


axes.set_ylabel('Spectral Index (arb. u.)', fontsize='14')
axes.set_xlabel('Drift Velocity ($c$)', fontsize='14')

# axes.set_xscale('log')
# axes.set_yscale('log')

axes.tick_params(which='both', bottom=True, top=True, left=True, right=True, direction='in', labelsize='14', length=5)

axes.set_xlim(0.1,0.6)
axes.set_ylim(4.5,6.6)

plt.legend(fontsize=14, loc='upper right', frameon=False)

#fig.savefig("/Users/thomas/SS Research Project/week 7/sp_in_vs_v.png")

plt.show()

# +
# Plotting spectral index vs. drift velocity fits to data for each model 

x_val = np.linspace(0.1,0.6,200)

fig, axes = plt.subplots(figsize=(8,5))

# with error bars...

plt.errorbar(v_s, spec_ind,  yerr=sp_ind_err, xerr=dv_s, linestyle='', color='black', capsize=2, zorder=0)
plt.errorbar(v_ms, spec_ind,  yerr=sp_ind_err, xerr=dv_ms, linestyle='', color='black', capsize=2, zorder=0)
plt.errorbar(v_ss, spec_ind,  yerr=sp_ind_err, xerr=dv_ss, linestyle='', color='black', capsize=2, zorder=0)

plt.scatter(v_s, spec_ind, color='r', label = 'Newkirk Model', zorder=5)
plt.scatter(v_ms, spec_ind, color='b', label = 'Mann Model', zorder=5)
plt.scatter(v_ss, spec_ind, color='g', label = 'Saito Model', zorder=5)


plt.plot(x_val, fit(x_val, *popt_n), color='red', linestyle='--', linewidth=2, label='Fit to Newkirk', zorder=10)
plt.plot(x_val, fit(x_val, *popt_m), color='blue', linestyle='--', linewidth=2, label='Fit to Mann', zorder=10)
plt.plot(x_val, fit(x_val, *popt_s), color='green', linestyle='--', linewidth=2, label='Fit to Saito', zorder=10)


axes.set_ylabel('Spectral Index (arb. u.)', fontsize='14')
axes.set_xlabel('Drift Velocity ($c$)', fontsize='14')

# axes.set_xscale('log')
# axes.set_yscale('log')

axes.tick_params(which='both', bottom=True, top=True, left=True, right=True, direction='in', labelsize='14', length=5)

axes.set_xlim(0.1,0.6)
axes.set_ylim(4.5,6.6)

plt.legend(fontsize=14, loc='upper right', frameon=False)

#fig.savefig("/Users/thomas/SS Research Project/week 7/sp_in_vs_v.png")

plt.show()
# -

print('\n', popt_n[0],'+/-', fit_err_n[1], 
      '\n', popt_m[0],'+/-', fit_err_m[1], 
      '\n', popt_s[0],'+/-', fit_err_s[1])

# +

print('\n', popt_n[1],'+/-', fit_err_n[1], '\n', popt_m[1],'+/-', fit_err_m[1])
# -


