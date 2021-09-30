from numpy import sum, reshape
from scipy.fft import fft

def ord_MSC(sinal, tamanhoJanela, M):
    sinal = reshape(sinal[0:tamanhoJanela*M], (M,tamanhoJanela))

    SINAL = fft(sinal).transpose()
    
    MSC = (abs(sum(SINAL,axis=1))**2) / (M*sum(abs(SINAL)**2,axis=1))

    return [SINAL,MSC]