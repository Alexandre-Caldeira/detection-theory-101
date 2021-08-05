from numpy import mean, quantile

def limiar_MEDIA(espontaneo,alpha):
    # media do sinal espontâneo (ruído)
    return quantile(mean(espontaneo), 1-alpha)