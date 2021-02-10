from obspy import read
from obspy import Stream
import numpy as np
import os
import sys
sys.path.append("/Users/ezgikarasozen/Documents/Research/Array_processing/spectral_subtraction/")
import ss_functions as func


filename = input("Parameters file: ")
exec(open(filename).read())
st = read(event_file)
stp = st.copy() #processed stream
stp.clear()
nos = len(st)
for s in range(nos): #station loop
    t = st[s].stats.starttime # does this stay here or moved to second loop? does it matter?
    for i in range(0, 12): # 5 minutes increments; later define this in the parameters file
        trd = st[s].copy() #[d]egraded data (input)
        trn = st[s].copy() #[n]oise signal (estimated)
        trd.trim((t + (i*tws*60)), (t + ((i+1)*tws*60))) # 5 minutes of signal to perform spectral subtraction on
        trn.trim((t + (i*twn*60)), (t + ((i+1)*twn*60))) # 1 minutes of signal to estimate noise from
#        print(trd.stats.starttime)
        trn.detrend("linear")
        trn.detrend("demean")
        trp = trd.copy() #[p]rocessed signal (output)
        amp_Xd, amp_Xn, Xd, scales_d, dt, dj = func.do_wavelet(trd,trn) #cwt
        amp_Xp, alpha, rho, a, b = func.do_subtraction(amp_Xd,amp_Xn, 0.75, 0.005) # spectral subtraction
        trp.data = func.do_iwavelet(Xd,amp_Xp,trp,scales_d,dt,dj)      
        outpath = 'results/'   #save processed waveforms
        if not os.path.exists(outpath):
           os.makedirs(outpath)
        print(trp.stats.starttime)
        print(trp.stats.endtime)
        trn.write(outpath + trn.stats.station + "." + trn.stats.channel + "." + str(trn.stats.starttime.year) + "." + str(trn.stats.starttime.month) + "." + str(trn.stats.starttime.day) +
        "." + str(trn.stats.starttime.hour) + "." + str(trn.stats.starttime.minute) + "." + str(i) +".noise", format="MSEED")
        trp.write(outpath + trp.stats.station + "." + trp.stats.channel + "." + str(trp.stats.starttime.year) + "." + str(trp.stats.starttime.month) + "." + str(trp.stats.starttime.day) +
        "." + str(trp.stats.starttime.hour) + "." + str(trp.stats.starttime.minute) + "." + str(i) +".processed", format="MSEED")
    stp.append(trp)  
    #st[s].merge(method=1)
print(stp)
stp.write("trial.mseed", format="MSEED")