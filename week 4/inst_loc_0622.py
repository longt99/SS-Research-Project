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
from sunpy.coordinates import frames, get_body_heliographic_stonyhurst, get_horizons_coord
import sunpy.timeseries
import sunpy.map
from sunpy.time import parse_time

from astropy import units as u
from astropy.coordinates import SkyCoord

import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np 

from astropy import constants as const

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
# -

earth_pos = get_body_heliographic_stonyhurst("earth", "2021-06-22")
mars_pos = get_body_heliographic_stonyhurst("mars", "2021-06-22")
solo_coord = get_horizons_coord('Solar Orbiter', "2021-06-30")

obstime_seq = parse_time('2021-06-22 04:30')


earth_seq = SkyCoord(get_body_heliographic_stonyhurst("earth", obstime_seq))
sun_seq = SkyCoord(get_body_heliographic_stonyhurst("sun", obstime_seq))
solo_coord_seq = get_horizons_coord('solo', obstime_seq)

solo_coord_seq

# +
travel_distance = (earth_seq.radius - solo_coord_seq.radius).to(u.m)

time_offset = travel_distance/const.c

print(time_offset)
# -

solo_coord_seq.heliocentricinertial.lon.to(u.rad)

# +
fig = plt.figure(dpi=120, figsize=(8,5))

ax = plt.subplot(projection='polar')

ax.set_facecolor('black')

ax.set_xticks(np.arange(0,2.0*np.pi,np.pi/6.0))
ax.tick_params(axis='x', colors='black')


ax.set_rlim(0,1.2)
ax.set_rticks(np.arange(0,1.2,0.2))
ax.tick_params(axis='y', colors='lightgrey')
ax.set_rlabel_position(0)

ax.plot(earth_seq.heliocentricinertial.lon.to(u.rad), earth_seq.heliocentricinertial.distance, marker='o', markersize = 7, linewidth = 0, label="Earth", color = 'blue')
ax.plot(sun_seq.heliocentricinertial.lon.to(u.rad), sun_seq.heliocentricinertial.distance, marker='o', markersize = 14, linewidth = 0, label="Sun", color="gold")
ax.plot(solo_coord_seq.heliocentricinertial.lon.to(u.rad), solo_coord_seq.heliocentricinertial.distance, 'o', markersize = 7, label='Solar Orbiter', color='red')

plt.annotate('', xy=(-2.8829,0.99), xytext=(-2.8829,0.04), arrowprops=dict(arrowstyle="<|-|>", color='white'))
plt.annotate('', xy=(1.6535,0.92), xytext=(1.6535,0.04), arrowprops=dict(arrowstyle="<|-|>", color='white'))

# plt.annotate('', xy=(earth_seq.heliocentricinertial.lon.to(u.rad),0.99), xytext=(earth_seq.heliocentricinertial.lon.to(u.rad),0.05), arrowprops=dict(arrowstyle="<->"))
# plt.annotate('', xy=(sun_seq.heliocentricinertial.lon.to(u.rad),0.85), xytext=(sun_seq.heliocentricinertial.lon.to(u.rad),0.05), arrowprops=dict(arrowstyle="<->"))




ax.legend(bbox_to_anchor=(1.4,1), loc='upper right', fontsize=10, facecolor = 'lightgrey')
#ax.set_title("Positions in HCI frame", fontsize=16)

fig.savefig("/Users/thomas/SS Research Project/week 4/event_1/polar_plot_0622.png")


# -

ax.get_rlabel_position()


