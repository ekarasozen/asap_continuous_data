from obspy import read
from obspy import Stream
import numpy as np
import os
import sys
import ss_functions as func

filename = input("Parameters file: ")
exec(open(filename).read())
st = read(event_file)
stp = Stream() #processed stream
nos = len(st)
for s in range(nos): #station loop
    t = st[s].stats.starttime 
    inc = (st[s].stats.npts/st[s].stats.sampling_rate)/(tw) #calculate how many 5 minute time increments are needed for this data 
    inc = inc/24 #for now, 1 hour, comment out for the entire data
    for i in range(0, int(inc)): # 5 minutes increments
        trd = st[s].copy() #[d]egraded data (input)
        trn = st[s].copy() #[n]oise signal (estimated)
        trd.trim((t + (i*tw)), (t + ((i+1)*tw))) # 5 minutes of signal to perform spectral subtraction on
        trn.trim((t + (i*tw)), (t + ((i+1)*tw))) # 5 minutes of signal to estimate noise from
        trn.detrend("linear")
        trn.detrend("demean")
        trp = trd.copy() #[p]rocessed signal (output)
        amp_Xd, amp_Xn, Xd, scales_d, dt, dj = func.do_wavelet(trd,trn) #cwt
        amp_Xp, alpha, rho, a, b = func.do_subtraction(amp_Xd,amp_Xn, 0.75, 0.005) # spectral subtraction
        trp.data = func.do_inv_wavelet(Xd,amp_Xp,trp,scales_d,dt,dj)      
        #save trimmed files for plotting, these following lines will be deleted in the final version of the code. 
        outpath = 'results/'   #save processed waveforms
        if not os.path.exists(outpath):
           os.makedirs(outpath)
        trn.write(outpath + trn.stats.station + "." + trn.stats.channel + "." + str(trn.stats.starttime.year) + "." + str(trn.stats.starttime.month) + "." + str(trn.stats.starttime.day) +
        "." + str(trn.stats.starttime.hour) + "." + str(trn.stats.starttime.minute) + "." + str(i) +".noise", format="MSEED")
        trp.write(outpath + trp.stats.station + "." + trp.stats.channel + "." + str(trp.stats.starttime.year) + "." + str(trp.stats.starttime.month) + "." + str(trp.stats.starttime.day) +
        "." + str(trp.stats.starttime.hour) + "." + str(trp.stats.starttime.minute) + "." + str(i) +".processed", format="MSEED")
        stp += trp
    stp.sort(['station']) #otherwise BR106 becomes the first trace in the stream
    stp.merge(method=1)
stp.write(output_file, format="MSEED")