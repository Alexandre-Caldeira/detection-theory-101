from numpy import pi, sin, cos, linspace

def geraSenoide(
    amplitude = 1,      # Valor de pico
    theta0 = 0,         # Fase em t=0
    f0 = 100,           # Frequência do sinal
    fs = 1000,          # Frequência de amostragem
    N = -1,             # Número de amostras desejado
    M = -1,             # Número de janelas desejado
    tj = -1,            # Tamanho de janelas (em amostras)
    cosseno = -1,       # Trocar para cosseno
    ):

    fl = 2*pi*f0/fs     # Frequência observada

    # Se informado, calcular número de pontos através do 
    # número e tamanho de janelas desejado:
    if M>0 and tj>0: N = M*tj

    # Gerar tempo discreto:
    t = linspace(0,N-1,num=N)
    arg = fl*t

    # Calcular pontos da senóide e retornar:
    if cosseno<0: sinal = amplitude*sin(arg +theta0)
    else: sinal = amplitude*cos(arg +theta0); print('Aviso: cosseno não validado!');

    return sinal