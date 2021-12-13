function [ORD] = LFT(y,L,bin)
%
%Local spectral F-test.   
Y = fft(y);
Yfo = Y(bin,1);
Yfn = Y([(bin-L/2):(bin-1);(bin+1):(bin+L/2)],1);
ORD = (abs(Yfo).^2)./(1/L*sum(abs(Yfn).^2,1));



