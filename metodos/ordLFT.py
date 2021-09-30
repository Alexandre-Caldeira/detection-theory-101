from numpy import abs, delete, floor, sum
from scipy.fft import fft

def ord_SFT(sinal, L, BIN):
    if (round(BIN) +L/2 ) > len(sinal): print('warning!')

    # L: tamanho das laterais
    SINAL = abs(fft(sinal))

    central = round(BIN) # certifica inteiro (não float)
    DEN =(SINAL[central])**2

    lateralMenor = round(central-L/2)
    lateralMaior = round(central+L/2)       
    SINAL_lateral = SINAL[lateralMenor:lateralMaior+1]

    if len(SINAL_lateral)>0:    
        SINAL_lateral = delete(SINAL_lateral,round(floor(len(SINAL_lateral)/2)))
    else: print('len(SINAL_lateral) = 0, revisar:','L',L,'bin',BIN,'!')

    NUM = (1/L)*sum((SINAL_lateral)**2)

    SFT = DEN/NUM
    
    return [SINAL,SFT]