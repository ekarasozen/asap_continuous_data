from obspy import read
from obspy import Stream
from obspy import Trace
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates

sto = read("../keskin_data/2018/*BR10**.064.00.00.00")
stp = read("results/keskin_processed_2018_064")
print(stp)


fig = plt.figure(figsize=(8,8))
gs = gridspec.GridSpec(6, 1)
gs.update(wspace=0.01, hspace=0.01) # set the spacing between axes. 
fig.suptitle('Original waveforms')


for i in range(len(sto)):
    ax0 = fig.add_subplot(gs[i])
    t = sto[i].stats.starttime # does this stay here or moved to second loop? does it matter?
    sto[i].trim((t), (t + (12*5*60))) # 5 minutes of signal to perform spectral subtraction on
    ax0.plot(sto[i].times("matplotlib"), sto[i].data, 'k', linewidth=0.2)
    text = (sto[i].stats.station)
    ax0.text(0.85, 0.90, text, transform=ax0.transAxes, fontsize=10, fontweight='bold', color='blue', verticalalignment='top')
    ax0.xaxis_date()
    ax0.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
fig.autofmt_xdate()
plt.tight_layout()
plt.draw()
fig.savefig('plots/original_enitre_data_wf.png', bbox_inches='tight')



fig1 = plt.figure(figsize=(8,8))
gs1 = gridspec.GridSpec(6, 1)
gs1.update(wspace=0.01, hspace=0.01) # set the spacing between axes. 
fig1.suptitle('Processed waveforms')
for i in range(len(stp)):
    ax1 = fig1.add_subplot(gs1[i])
    ax1.plot(stp[i].times("matplotlib"), stp[i].data, 'k', linewidth=0.2)
    text = (stp[i].stats.station)
    ax1.text(0.85, 0.90, text, transform=ax1.transAxes, fontsize=10, fontweight='bold', color='blue', verticalalignment='top')
    ax1.xaxis_date()
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
fig1.autofmt_xdate()
plt.tight_layout()
plt.draw()
fig1.savefig('plots/processed_entire_data_wf.png', bbox_inches='tight')
