from numpy import quantile

def vc_pratico(ORD, alpha = 0.05):
    # Calcula quantil a partir de 
    # limiar de significância desejado:
    return quantile(a= ORD, q = 1-alpha)