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
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy import timeseries as ts

import astropy.units as u
# -

result1 = Fido.search(a.Time('2021/08/28 21:57:00', '2021/08/28 21:59:00'),
                     a.Instrument.hmi, a.Physobs.los_magnetic_field)

downloaded_file1 = Fido.fetch(result1)
print(downloaded_file1)

# +
hmi_map = sunpy.map.Map(downloaded_file1[1])

hmi_map.plot_settings['cmap'] = "hmimag"
hmi_map.plot_settings['norm'] = plt.Normalize(-1500, 1500)
# -

result2 = Fido.search(a.Time('2021/08/28 21:57:00', '2021/08/28 21:59:00'),
                     a.Instrument.aia,
                     a.Physobs.intensity,
                     a.Wavelength(171*u.angstrom))

result2

downloaded_file2 = Fido.fetch(result2[0][0])
#print(downloaded_file2)

aia_map = sunpy.map.Map(downloaded_file2[0])

# +
fig = plt.figure(figsize=(12, 5))
ax1 = fig.add_subplot(1, 2, 1, projection=aia_map)
aia_map.plot(axes=ax1, clip_interval=(1, 99.9)*u.percent)
ax2 = fig.add_subplot(1, 2, 2, projection=aia_map)
hmi_map.plot(axes=ax2, autoalign=True, title='HMI image in AIA reference frame')
ax2.axis(ax1.axis())

#fig.savefig("/Users/thomas/SS Research Project/week 2/solar_monitor.png")
plt.show()

# -


