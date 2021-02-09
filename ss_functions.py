import numpy as np
from obspy import Stream
import pycwt

def do_wavelet(trd, trn): 
    dt = trd.stats.delta
    dj = 0.05 #scale spacing
    s0 =  0.096801331 # smallest scale
    omega0 = 6
    J = (np.log2(len(trd) * dt / s0)) / dj  # number of scales-1 
    mother = pycwt.Morlet(omega0)              # See https://github.com/regeirk/pycwt/blob/master/pycwt/wavelet.py 
    #########start cwt#########
    Xd, scales_d, freqs_d, coi_d, fft_d, fftfreqs_d = pycwt.cwt(trd.data, dt, dj, s0, J, mother)    #degraded signal wavelet transform
    amp_Xd = abs(Xd) #amplitude spectrum of the degraded signal for spectral subtraction
    Xn, scales_n, freq_n, coi_n, fft_n, fftfreqs_n = pycwt.cwt(trn.data, dt, dj, s0, J, mother)     #noise signal wavelet transform
    amp_Xn = abs(Xn) #amplitude spectrum of the noise signal for spectral subtraction
    return amp_Xd, amp_Xn, Xd, scales_d, dt, dj
    
def do_subtraction(amp_Xd, amp_Xn, a, b): 
    m, n = amp_Xd.shape 
    rho = np.zeros((m,n))
    amp_Xp = np.zeros((m,n))
    ones_array=np.ones((1,n))  
    amp_Xna = np.array([np.median(amp_Xn,axis=1)]).transpose() 
    amp_Xna = np.dot(amp_Xna,ones_array)
    rho = amp_Xd / amp_Xna
    alpha =  a * (rho * (1-np.tanh(b*rho**(4)))) 
    amp_Xp = amp_Xd - (alpha*amp_Xna)
    return amp_Xp, alpha, rho, a, b

def do_iwavelet(Xd,amp_Xp,trp,scales_d,dt,dj): #wavelet parameters (scales, dt, dj are taken from do_wavelet outputs)
    phase_Xd = np.angle(Xd) #phase spectrum of the processed signal, from degraded
    Xp = amp_Xp*(np.exp(1.j*phase_Xd)) #combine phase and amplitude spectrum of the processed signal
    trp.data = pycwt.icwt(Xp, scales_d, dt, dj, wavelet='morlet') #inverse cwt to compute processed signal 
    print('----> maximum imaginary value in "processed signal" is: ', np.max(np.imag(trp.data)))
    if np.max(np.imag(trp.data)) == 0:
        trp.data = np.real(trp.data)
    else: #this needs explanation & testing.
        print('maximum imaginary value in "processed signal" is not zero, therefore not outputted')
    return trp.data
