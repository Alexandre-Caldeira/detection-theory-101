function [ORD] = LFT_freq(y,L,fs,f)
%
%Local spectral F-test.   
tamanho_sinal = size(y,1); 
bin = round(f*tamanho_sinal/fs+1);

Y = fft(y);
Yfo = Y(bin,1);
Yfn = Y([(bin-L/2):(bin-1);(bin+1):(bin+L/2)],1);
ORD = (abs(Yfo).^2)./(1/L*sum(abs(Yfn).^2,1));



