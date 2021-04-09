from obspy import read, Stream
import numpy as np
import os
import sys
import ss_functions as func

filename = input("Parameters file: ")
exec(open(filename).read())
st = read(event_file)
stp = Stream() #processed stream
strs = Stream() #residual stream
nos = len(st)
st.detrend("linear")
st.detrend("demean")
if not os.path.exists(outpath):
   os.makedirs(outpath)
freqs_array = np.empty((233,)) #233 is freq number, not sure how to hardwire this from do_wavelet
for s in range(nos): #station loop
    t = st[s].stats.starttime 
    inc = int((st[s].stats.npts/st[s].stats.sampling_rate)/(tw)) #calculate how many 5 minute time increments are needed for this data 
    #inc = inc/24 #for now, 1 hour, comment out for the entire data
    windowlength = int(tw*st[s].stats.sampling_rate)
    noise_estimate = np.zeros([233,inc*windowlength])
  #  noise_estimate = np.empty((233,0)) #233 is freq number, not sure how to hardwire this from do_wavelet
    for i in range(0, int(inc)): # 5 minutes increments
        trd = st[s].copy() #[d]egraded data (input)
        trd.trim((t + (i*tw)), (t + ((i+1)*tw) - trd.stats.delta)) # 5 minutes of signal to perform spectral subtraction on, stat.delta is necessary to have equal 6000 pts trimmed data. 
        trp = trd.copy() #[p]rocessed signal (output)
        trrs = trd.copy() #residual data
        amp_Xd, Xd, freqs, scales, dt, dj = func.do_wavelet(trd) #cwt
        amp_Xp, amp_Xna, alpha, rho, a, b = func.do_subtraction(amp_Xd) # spectral subtraction
        trp.data = func.do_inv_wavelet(Xd,amp_Xp,trp,scales,dt,dj)      
        trrs.data = trd.data-trp.data
        noise_estimate[:,i*windowlength:i*windowlength+windowlength] = amp_Xna
        #noise_estimate = np.concatenate((noise_estimate, amp_Xna), axis=1) #you can delete this once you do a successful run
        #save trimmed files for plotting, these following lines will be deleted in the final version of the code. 
        ##trrs.write(outpath + output_name + '_' + trrs.stats.station + '_' + str(i) + "_residual", format="MSEED")
        ##trp.write(outpath + output_name + '_' +trp.stats.station + '_' + str(i) + "_processed", format="MSEED")
        ##np.savetxt(outpath + output_name + '_' + trp.stats.station + '_' + str(i)  + '_noise_estimate', amp_Xna, delimiter=',', newline="\n")  
        stp += trp
        strs += trrs
    stp.merge(method=1)
    stp.sort(['station']) #otherwise BR106 becomes the first trace in the stream
    strs.merge(method=1)
    strs.sort(['station']) #otherwise BR106 becomes the first trace in the stream
    np.save(outpath + 'noise_estimate_' + trp.stats.station, noise_estimate)  
freqs_array = freqs
np.save(outpath + 'freqs_array_noise_estimate', freqs_array)
stp.write(outpath + output_name + '_processed', format="MSEED")
strs.write(outpath + output_name + '_residual', format="MSEED")