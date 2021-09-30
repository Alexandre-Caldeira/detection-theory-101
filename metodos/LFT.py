from numpy import abs, delete, floor, sum
from scipy.fft import fft

def ordLFT(sinal, L, BIN):
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


from scipy.stats import f as FDIST
def vcLFT(L, alpha = 0.05, VERBOSE = 0):
    # valor crítico teórico, dado tomando inversa da distribuição acumulada F
    # (Fisher-Snedecor)
    # com significancia 1-alpha (95% se não for alterada)
    vc = FDIST.ppf(1-alpha, 2,2*L) 
    
    if VERBOSE==1:
        print('Significância desejada: ',alpha*100,'\b%')
    if VERBOSE==2:
        print('Significância desejada: ',alpha*100,'\b%')
        print('Valor crítico SFT: ',vc)

    return vc