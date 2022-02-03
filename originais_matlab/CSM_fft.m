function [ORD] = CSM_fft(Y,M)
teta = angle(Y); %pegar o angulo; 
ORD =  (1/M^2)*sum(cos(teta),2).^2+(1/M^2)*sum(sin(teta),2).^2;


