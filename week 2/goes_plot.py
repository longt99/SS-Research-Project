# %matplotlib inline

# +
import matplotlib.pyplot as plt
import numpy as np

from sunpy import timeseries as ts
from sunpy.net import Fido
from sunpy.net import attrs as a

# +
tstart = "2021-08-22"
tend = "2021-09-04"

# tstart = "2021-08-22 00:00"
# tend = "2021-08-22 23:59"

result_goes16 = Fido.search(a.Time(tstart, tend), a.Instrument("XRS"), a.goes.SatelliteNumber(16))
print(result_goes16)
# -

file_goes16 = Fido.fetch(result_goes16)

# +
import matplotlib.dates as mdates

sample_path = '/Users/thomas/SS Research Project/week 2'

goes_conc = []

for file in file_goes16:
    t = ts.TimeSeries(file)
    goes_conc.append(t)

# +
fig = plt.figure(figsize = (18,10)) # figsize = (12,18)

for timeseries in goes_conc:
    timeseries.plot(color=['black','red'])

#goes_conc.plot()

plt.yscale('log')
plt.ylim(10**(-8),10**(-4))
plt.xlim('2021-08-22','2021-09-05')
plt.title('GOES XRS PLOT FROM ' + tstart + ' to ' + tend, fontsize = 15)
plt.grid(axis='y')

plt.legend(['1-8 A', '0.4-5 A'], loc='upper right')

fig.savefig(sample_path + "/solar_monitor_goes.png")
# -


