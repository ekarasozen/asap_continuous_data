from obspy import read
from obspy import Stream
from obspy import Trace
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import matplotlib.cm as cm
import pycwt

std = read("../keskin_data/2018/*BR10**.064.00.00.00")
stp = read("results/*SHZ.2018.3.5.0.0.0.processed*")
stn = read("results/*SHZ.2018.3.5.0.0.0.noise*")

fig = plt.figure(figsize=(8,8))
gs = gridspec.GridSpec(6, 1)
gs.update(wspace=0.01, hspace=0.01) # set the spacing between axes. 
fig.suptitle('5 minute window original waveforms')
for i in range(len(std)):
    ax0 = fig.add_subplot(gs[i])
    t = std[i].stats.starttime # does this stay here or moved to second loop? does it matter?
    std[i].trim((t), (t + (1*5*60))) # 5 minutes of signal to perform spectral subtraction on
    ax0.plot(std[i].times("matplotlib"), std[i].data, 'k', linewidth=0.2)
    text = (std[i].stats.station)
    ax0.text(0.85, 0.90, text, transform=ax0.transAxes, fontsize=10, fontweight='bold', color='blue', verticalalignment='top')
    ax0.xaxis_date()
    ax0.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
fig.autofmt_xdate()
plt.tight_layout()
plt.draw()
fig.savefig('plots/windowed_original_waveforms.png', bbox_inches='tight')



fig1 = plt.figure(figsize=(8,8))
gs1 = gridspec.GridSpec(6, 1)
gs1.update(wspace=0.01, hspace=0.01) # set the spacing between axes. 
fig1.suptitle('5 minute processed waveforms')
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
fig1.savefig('plots/windowed_processed_waveforms.png', bbox_inches='tight')


fig2 = plt.figure(figsize=(8,8))
gs2 = gridspec.GridSpec(6, 1)
gs2.update(wspace=0.01, hspace=0.01) # set the spacing between axes. 
fig2.suptitle('5 minute window noise sample')
for i in range(len(stn)):
    ax2 = fig2.add_subplot(gs2[i])
    ax2.plot(stn[i].times("matplotlib"), stn[i].data, 'k', linewidth=0.2)
    text = (stn[i].stats.station)
    ax2.text(0.85, 0.90, text, transform=ax2.transAxes, fontsize=10, fontweight='bold', color='blue', verticalalignment='top')
    ax2.xaxis_date()
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
fig2.autofmt_xdate()
plt.tight_layout()
plt.draw()
fig2.savefig('plots/windowed_noise_waveforms.png', bbox_inches='tight')


fig3 = plt.figure(figsize=(8,8))
gs3 = gridspec.GridSpec(6, 1)
gs3.update(wspace=0.01, hspace=0.01) # set the spacing between axes. 
fig3.suptitle('5 minute window original scaleograms')
fig4 = plt.figure(figsize=(8,8))
gs4 = gridspec.GridSpec(6, 1)
gs4.update(wspace=0.01, hspace=0.01) # set the spacing between axes. 
fig4.suptitle('5 minute window processed scaleograms')
fig5 = plt.figure(figsize=(8,8))
gs5 = gridspec.GridSpec(6, 1)
gs5.update(wspace=0.01, hspace=0.01) # set the spacing between axes. 
fig5.suptitle('5 minute window noise sample scaleograms')
for i in range(len(std)):
    td = np.arange(std[i].stats.npts) / std[i].stats.sampling_rate
    tp = np.arange(stp[i].stats.npts) / stp[i].stats.sampling_rate
    tn = np.arange(stn[i].stats.npts) / stn[i].stats.sampling_rate
    ax11 = fig3.add_subplot(gs3[i])
    ax12 = fig4.add_subplot(gs4[i])
    ax13 = fig5.add_subplot(gs5[i])
    dtd = std[i].stats.delta
    dtp = stp[i].stats.delta
    dtn = stp[i].stats.delta
    dj = 0.05 #scale spacing
    s0 =  0.096801331 # smallest scale
    omega0 = 6
    Jd = (np.log2(len(std[i]) * dtd / s0)) / dj  # number of scales-1 
    Jp = (np.log2(len(stp[i]) * dtp / s0)) / dj  # number of scales-1 
    Jn = (np.log2(len(stn[i]) * dtn / s0)) / dj  # number of scales-1 
    mother = pycwt.Morlet(omega0)              # See https://github.com/regeirk/pycwt/blob/master/pycwt/wavelet.py 
    Xd, scales_d, freqs_d, coi_d, fft_d, fftfreqs_d = pycwt.cwt(std[i].data, dtd, dj, s0, Jd, mother)    #degraded signal wavelet transform
    Xp, scales_p, freqs_p, coi_p, fft_p, fftfreqs_p = pycwt.cwt(stp[i].data, dtp, dj, s0, Jp, mother)    #degraded signal wavelet transform
    Xn, scales_n, freqs_n, coi_n, fft_n, fftfreqs_n = pycwt.cwt(stn[i].data, dtn, dj, s0, Jp, mother)    #degraded signal wavelet transform
    amp_Xd = abs(Xd) #amplitude spectrum of the degraded signal for spectral subtraction
    amp_Xp = abs(Xp) #amplitude spectrum of the degraded signal for spectral subtraction
    amp_Xn = abs(Xn) #amplitude spectrum of the degraded signal for spectral subtraction
    maxamp_d = abs(Xd).max()/2           # these two lines just adjust the color scale
    maxamp_p = abs(Xp).max()/2           # these two lines just adjust the color scale
    maxamp_n = abs(Xn).max()/2           # these two lines just adjust the color scale
    minamp = 0
    tXd, fd = np.meshgrid(std[i].times(), freqs_d)
    tXp, fp = np.meshgrid(stp[i].times(), freqs_p)
    tXn, fn = np.meshgrid(stn[i].times(), freqs_n)
    img_d = ax11.pcolormesh(tXd, fd, np.abs(Xd), cmap=cm.hot, vmin=minamp, vmax=maxamp_p)
    img_p = ax12.pcolormesh(tXp, fp, np.abs(Xp), cmap=cm.hot, vmin=minamp, vmax=maxamp_p)
    img_n = ax13.pcolormesh(tXn, fn, np.abs(Xn), cmap=cm.hot, vmin=minamp, vmax=maxamp_p)
    ax11.set_xlim(td[0], td[-1])
    ax11.set_ylim(freqs_d[-1], freqs_d[0])
    ax11.set_xlabel("time after %s [s]" % std[i].stats.starttime)
    ax12.set_xlim(tp[0], tp[-1])
    ax12.set_ylim(freqs_p[-1], freqs_p[0])
    ax12.set_xlabel("time after %s [s]" % stp[i].stats.starttime)
    ax13.set_xlim(tn[0], tn[-1])
    ax13.set_ylim(freqs_n[-1], freqs_n[0])
    ax13.set_xlabel("time after %s [s]" % stn[i].stats.starttime)
fig3.autofmt_xdate()
plt.tight_layout()
plt.draw()
fig3.savefig('plots/windowed_original_scaleograms.png', bbox_inches='tight')
fig4.autofmt_xdate()
plt.tight_layout()
plt.draw()
fig4.savefig('plots/windowed_processed_scaleograms.png', bbox_inches='tight')
fig5.autofmt_xdate()
plt.tight_layout()
plt.draw()
fig5.savefig('plots/windowed_noise_scaleograms.png', bbox_inches='tight')






