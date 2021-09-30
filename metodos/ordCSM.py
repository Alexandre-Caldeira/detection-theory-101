from numpy import angle, unwrap, sin, cos, sum, reshape
from scipy.fft import fft

def ord_CSM(sinal, tamanhoJanela, M):
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