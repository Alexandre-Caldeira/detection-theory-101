from numpy import sum, reshape
from scipy.fft import fft

def ordMSC(sinal, tamanhoJanela, M):
    sinal = reshape(sinal[0:tamanhoJanela*M], (M,tamanhoJanela))

    SINAL = fft(sinal).transpose()
    
    MSC = (abs(sum(SINAL,axis=1))**2) / (M*sum(abs(SINAL)**2,axis=1))

    return [SINAL,MSC]
    
def vcMSC(M, alpha = 0.05, VERBOSE = 0):
    vc = 1-alpha**(1/(M-1))
    if VERBOSE==1:
        print('Significância desejada: ',alpha*100,'\b%')
    if VERBOSE==2:
        print('Significância desejada: ',alpha*100,'\b%')
        print('Valor crítico MSC: ',vc)
    return vc # valor crítico