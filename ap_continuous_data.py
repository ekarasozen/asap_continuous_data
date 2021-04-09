import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import numpy as np
from obspy import read, Stream, Trace
from obspy.core.util import AttribDict
from obspy.imaging.cm import obspy_sequential
from obspy.signal.array_analysis import array_processing
import os

outpath = '2018/90/'   #save processed waveforms
if not os.path.exists(outpath):
   os.makedirs(outpath)

# Load data
#st = read("../keskin_data/2018/*BR10**.090.00.00.00")
st = read(outpath + "ss_prcs_2018_090.mseed")
day = "_2018_090"
#dtype =  "dgrd"
dtype =  "prcs"
# Set coordinates for all 5 channels
st[0].stats.coordinates = AttribDict({
    'latitude': 39.725346,
    'elevation': 1463.7,
    'longitude': 33.639096})    
st[1].stats.coordinates = AttribDict({
    'latitude': 39.735643,
    'elevation': 1552.4,
    'longitude': 33.648541})    
st[2].stats.coordinates = AttribDict({
    'latitude': 39.719633,
    'elevation': 1515.0,
    'longitude': 33.657081})    
st[3].stats.coordinates = AttribDict({
    'latitude': 39.707446,
    'elevation': 1374.9,
    'longitude': 33.641485})    
st[4].stats.coordinates = AttribDict({
    'latitude': 39.718039,
    'elevation': 1458.5,
    'longitude': 33.617315})    
st[5].stats.coordinates = AttribDict({
    'latitude': 39.733703,
    'elevation': 1400.2,
    'longitude': 33.617991})    

# Execute array_processing
stime = st[0].stats.starttime
etime = st[0].stats.endtime
#etime = stime + 61 #comment out to test minor edits quickly
#st.trim((stime), (stime + (1*1*62)) - st[0].stats.delta) #comment out to test minor edits quickly
kwargs = dict(
   sll_x=-0.5, slm_x=0.5, sll_y=-0.5, slm_y=0.5, sl_s=0.005,    # slowness grid:   min,   max, Y min, Y max, Slow Step
   win_len=6.0, win_frac=0.083333334,    # every 0.5 sec
   frqlow=0.5, frqhigh=6, prewhiten=0,    # frequency properties
   semb_thres=-1e9, vel_thres=-1e9, timestamp='mlabday',    # restrict output
   stime=stime, etime=etime, 
   method=0) #the method to use 0 == bf, 1 == capon
out = array_processing(st, **kwargs)

#Save output as a numpy array
ap_params = out[:, :]
np.save(outpath + 'ap_' + dtype + day, ap_params) #Column orders: Time, Rel. Power, Abs. Power, BAZ, Slowness
#Save output as a MSEED 
rp = ap_params[:,1]
ap = ap_params[:,2]
baz = ap_params[:,3]
slw = ap_params[:,4]
#PRCS
#DGRD
rp_stats = {'network': '', 'station': 'PRCS', 'location': '', 'channel': 'rpw', 'npts': len(rp), 'sampling_rate': 2, 'mseed': {'dataquality': 'D'}}
rp_stats['starttime'] = stime
ap_stats = {'network': '', 'station': 'PRCS', 'location': '', 'channel': 'apw', 'npts': len(ap), 'sampling_rate': 2, 'mseed': {'dataquality': 'D'}}
ap_stats['starttime'] = stime
baz_stats = {'network': '', 'station': 'PRCS', 'location': '', 'channel': 'baz', 'npts': len(baz), 'sampling_rate': 2, 'mseed': {'dataquality': 'D'}}
baz_stats['starttime'] = stime
slw_stats = {'network': '', 'station': 'PRCS', 'location': '', 'channel': 'slw', 'npts': len(slw), 'sampling_rate': 2, 'mseed': {'dataquality': 'D'}}
slw_stats['starttime'] = stime
stap = Stream([Trace(data=rp, header=rp_stats), Trace(data=ap, header=ap_stats), Trace(data=baz, header=baz_stats), Trace(data=slw, header=slw_stats)])
stap.write(outpath + "ap_" + dtype + day + '.mseed', format='MSEED', encoding="FLOAT64")
stap = read(outpath + "ap_" + dtype + day + '.mseed')
print(stap)

# Plot
fig1 = plt.figure()
labels = ['rel_power', 'abs_power', 'baz', 'slowness']
xlocator = mdates.AutoDateLocator()
for i, lab in enumerate(labels):
     ax1 = fig1.add_subplot(4, 1, i + 1)
     ax1.scatter(out[:, 0], out[:, i + 1], c=out[:, 1], alpha=0.6,
                edgecolors='none', cmap=obspy_sequential)
     ax1.set_ylabel(lab)
     ax1.set_xlim(out[0, 0], out[-1, 0]) 
     ax1.set_ylim(out[:, i + 1].min(), out[:, i + 1].max()) 
     ax1.xaxis.set_major_locator(xlocator)
     ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
fig1.suptitle('Keskin Array %s' % (stime.strftime('%Y-%m-%d'), ))
fig1.subplots_adjust(left=0.15, top=0.95, right=0.95, bottom=0.2, hspace=0)
fig1.autofmt_xdate()
plt.draw()
fig1.savefig(outpath + 'ap_' + dtype + day + '_scatter.png', bbox_inches='tight')

#Plot array parameter mseed data
fig2 = plt.figure(figsize=(8,8))
gs = gridspec.GridSpec(4, 1)
gs.update(wspace=0.01, hspace=0.01) # set the spacing between axes. 
fig2.suptitle('Keskin Array %s' % (stime.strftime('%Y-%m-%d'), ))

for i in range(len(stap)):
    ax2 = fig2.add_subplot(gs[i])
    ax2.plot(stap[i].times("matplotlib"), stap[i].data, 'k', linewidth=0.2)
    text = (stap[i].stats.station)
    ax2.text(0.85, 0.90, text, transform=ax2.transAxes, fontsize=10, fontweight='bold', color='blue', verticalalignment='top')
    ax2.xaxis_date()
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
fig2.autofmt_xdate()
plt.tight_layout()
plt.draw()
fig2.savefig(outpath + 'ap_' + dtype + day + '_mseed.png', bbox_inches='tight')

