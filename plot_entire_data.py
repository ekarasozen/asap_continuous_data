from obspy import read
from obspy import Stream
from obspy import Trace
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import matplotlib.cm as cm
import pycwt
import ss_functions as func
import os

outpath = '2018/65/'   #save processed waveforms
if not os.path.exists(outpath):
   os.makedirs(outpath)

std = read("../keskin_data/2018/*BR10**.065.00.00.00")
stp = read(outpath + "/keskin_2018_065_processed")
strs = read(outpath + "/keskin_2018_065_residual")
print(stp)
print(strs)

fig = plt.figure(figsize=(8,8))
gs = gridspec.GridSpec(6, 1)
gs.update(wspace=0.01, hspace=0.01) 
fig.suptitle('Degraded waveforms')
fig1 = plt.figure(figsize=(8,8))
gs1 = gridspec.GridSpec(6, 1)
gs1.update(wspace=0.01, hspace=0.01) 
fig1.suptitle('Processed waveforms')
fig2 = plt.figure(figsize=(8,8))
gs2 = gridspec.GridSpec(6, 1)
gs2.update(wspace=0.01, hspace=0.01) 
fig2.suptitle('Residual waveforms')
t = std[0].stats.starttime 
#std.trim((t), (t + (12*5*60))-std[0].stats.delta) # when you want to work on 1 hour data. 
std.detrend("linear")
std.detrend("demean")
y_min = np.nanmin(std[:]) #for y axis to share same min and max values within all wfs
y_max = np.nanmax(std[:])

#Plot waveforms
for i in range(len(std)):
    ax0 = fig.add_subplot(gs[i])
    ax1 = fig1.add_subplot(gs1[i])
    ax2 = fig2.add_subplot(gs2[i])
    ax0.plot(std[i].times("matplotlib"), std[i].data, 'k', linewidth=0.2)
    ax1.plot(stp[i].times("matplotlib"), stp[i].data, 'k', linewidth=0.2)
    ax2.plot(strs[i].times("matplotlib"), strs[i].data, 'k', linewidth=0.2)
    text = (std[i].stats.station)
    ax0.text(0.85, 0.90, text, transform=ax0.transAxes, fontsize=10, fontweight='bold', color='blue', verticalalignment='top')
    ax1.text(0.85, 0.90, text, transform=ax1.transAxes, fontsize=10, fontweight='bold', color='blue', verticalalignment='top')
    ax2.text(0.85, 0.90, text, transform=ax2.transAxes, fontsize=10, fontweight='bold', color='blue', verticalalignment='top')
    ax0.set_ylim(ymin=y_min, ymax=y_max)
    ax1.set_ylim(ymin=y_min, ymax=y_max)
    ax2.set_ylim(ymin=y_min, ymax=y_max)
    ax0.xaxis_date()
    ax1.xaxis_date()
    ax2.xaxis_date()
    ax0.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
fig.autofmt_xdate()
fig1.autofmt_xdate()
fig2.autofmt_xdate()
plt.tight_layout()
plt.draw()
fig.savefig(outpath + 'degraded_wfs.png', bbox_inches='tight')
fig1.savefig(outpath + 'processed_wfs.png', bbox_inches='tight')
fig2.savefig(outpath + 'residual_wfs.png', bbox_inches='tight')

#Plot scaleograms

fig3 = plt.figure(figsize=(8,8))
gs3 = gridspec.GridSpec(6, 1)
gs3.update(wspace=0.01, hspace=0.01) # set the spacing between axes. 
fig3.suptitle('Time window:5 min., degraded')
fig4 = plt.figure(figsize=(8,8))
gs4 = gridspec.GridSpec(6, 1)
gs4.update(wspace=0.01, hspace=0.01) # set the spacing between axes. 
fig4.suptitle('Time window:5 min., processed')
fig5 = plt.figure(figsize=(8,8))
gs5 = gridspec.GridSpec(6, 1)
gs5.update(wspace=0.01, hspace=0.01) # set the spacing between axes. 
fig5.suptitle('Time window:5 min., residual')
fig6 = plt.figure(figsize=(8,8))
gs6 = gridspec.GridSpec(6, 1)
gs6.update(wspace=0.01, hspace=0.01) # set the spacing between axes. 
fig6.suptitle('Time window:5 min., noise')
freqs_ne = np.load(outpath + 'freqs_array_noise_estimate.npy')
for i in range(len(std)):
    amp_Xna = np.load(outpath + 'noise_estimate_BR10' + str(i+1) + '.npy')
    td = np.arange(std[i].stats.npts) / std[i].stats.sampling_rate
    ax11 = fig3.add_subplot(gs3[i])
    ax12 = fig4.add_subplot(gs4[i])
    ax13 = fig5.add_subplot(gs5[i])
    ax14 = fig6.add_subplot(gs6[i])
    amp_Xd, Xd, freqs, scales, dt, dj = func.do_wavelet(std[i]) #cwt
    amp_Xp, Xp, freqs, scales, dt, dj = func.do_wavelet(stp[i]) #cwt
    amp_Xrs, Xrs, freqs, scales, dt, dj = func.do_wavelet(strs[i]) #cwt
    maxamp = abs(Xp).max()/2             # these two lines just adjust the color scale
    minamp = 0
    tX, fX = np.meshgrid(std[i].times(), freqs)
    tXne, fXne = np.meshgrid(strs[i].times(), freqs_ne)
    img_d = ax11.pcolormesh(tX, fX, np.abs(Xd), cmap=cm.hot, vmin=minamp, vmax=maxamp)
    img_p = ax12.pcolormesh(tX, fX, np.abs(Xp), cmap=cm.hot, vmin=minamp, vmax=maxamp)
    img_rs = ax13.pcolormesh(tX, fX, np.abs(Xrs), cmap=cm.hot, vmin=minamp, vmax=maxamp)
    img_ne = ax14.pcolormesh(tXne, fXne, amp_Xna, cmap=cm.hot, vmin=minamp, vmax=maxamp)
    text = (std[i].stats.station)
    ax11.text(0.85, 0.90, text, transform=ax11.transAxes, fontsize=10, fontweight='bold', color='white', verticalalignment='top')
    ax12.text(0.85, 0.90, text, transform=ax12.transAxes, fontsize=10, fontweight='bold', color='white', verticalalignment='top')
    ax13.text(0.85, 0.90, text, transform=ax13.transAxes, fontsize=10, fontweight='bold', color='white', verticalalignment='top')
    ax14.text(0.85, 0.90, text, transform=ax14.transAxes, fontsize=10, fontweight='bold', color='white', verticalalignment='top')
    ax11.set_xlim(td[0], td[-1])
    ax11.set_ylim(freqs[-1], freqs[0])
    ax12.set_xlim(td[0], td[-1])
    ax12.set_ylim(freqs[-1], freqs[0])
    ax13.set_xlim(td[0], td[-1])
    ax13.set_ylim(freqs[-1], freqs[0])
    ax14.set_xlim(td[0], td[-1])
    ax14.set_ylim(freqs_ne[-1], freqs_ne[0])
    ax11.xaxis_date()
    ax12.xaxis_date()
    ax13.xaxis_date()
    ax14.xaxis_date()
    ax11.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax12.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax13.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax14.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
fig3.autofmt_xdate()
fig3.savefig(outpath + 'degraded_scaleograms.png', bbox_inches='tight', pad_inches = 0)
fig4.autofmt_xdate()
fig4.savefig(outpath + 'processed_scaleograms.png', bbox_inches='tight', pad_inches = 0)
fig5.autofmt_xdate()
fig5.savefig(outpath + 'residual_scaleograms.png', bbox_inches='tight', pad_inches = 0)
fig6.autofmt_xdate()
fig6.savefig(outpath + 'noise_estimate_scaleograms.png', bbox_inches='tight', pad_inches = 0)






