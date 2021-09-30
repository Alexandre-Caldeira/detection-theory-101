from numpy import abs, angle, unwrap
from scipy.fft import fft, fftfreq

def getFreq(sinal, fs):
    S = fft(sinal).transpose()
    SINALf = abs(S)
    SINALa = unwrap(angle(S))
    FREQs = fftfreq(len(SINALf), d=1/fs)
    return [FREQs,SINALf, SINALa]