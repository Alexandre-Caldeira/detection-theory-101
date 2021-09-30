from numpy import angle, unwrap, sin, cos, sum, reshape
from scipy.fft import fft

def ordCSM(sinal, tamanhoJanela, M):
    if len(sinal)-tamanhoJanela*M>=0:
        sinal = reshape(sinal[0:tamanhoJanela*M], (M,tamanhoJanela))

        FFT_SINAL = fft(sinal).transpose()
        angulo = unwrap(angle(FFT_SINAL))

        s = ((1/M)*(sum(sin(angulo),1)))**2
        c = ((1/M)*(sum(cos(angulo),1)))**2

        CSM =  c+s

        return [FFT_SINAL,CSM]
    else:
        print('Erro no número de janelas', round(tamanhoJanela),'(ou amostras, M =', M,'\b) escolhido.')
        print('Tamanho do Sinal menos M*tamanhoJanela:', len(sinal)-M*tamanhoJanela)
        print('(Retorna 0)')
        return 0 

from scipy.stats.distributions import chi2

def vcCSM(M, alpha = 0.05, VERBOSE=0):
    vc = chi2.ppf(1-alpha,2)/(2*M)
    
    if VERBOSE==1:
        print('Significância desejada: ',alpha*100,'\b%')
    if VERBOSE==2:
        print('Significância desejada: ',alpha*100,'\b%')
        print('Valor crítico CSM: ',vc)

    return vc