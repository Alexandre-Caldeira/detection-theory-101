function [MLphi, MLphicrit] = MLFT(y,L,fs,alpha,fo)
%
% Multivariate local F-test.
%
%Extension of the Local F-test to multi channel analysis. This is a MORD 
%technique that makes use only magnitude information and can be used as a 
%detector of hidden periodicities in noise, since N real-valued signals are
%available. The method is described in "Rocha, P.F.F., Felix, L.B.,
%Miranda de Sá, A.M.F.L. and Mendes, E.M.A.M. . Multivariate Evoked 
%Response Detection based on the Spectral F-test. Journal of Neuroscience 
%Methods, p. 113-118, 2016."
%
%Sintaxes:
%
%[MLphi, MLphicrit] = MLFT(y,fo,L,fs,alpha) => detetor's value at fo and
%respective theoretical critical value.
%
%Input Parameters:
%
%y => matrix whose columns are the signals to be considered. 
%fo => frequency to be tested.
%L => number of closest neighbouring frequencies to fo. L must be even.
%fs => sample rate.
%alpha => significance level, e.g. alpha = 0.05. 
% 
%Example
%
% fs=100;M=10;fo=21;L=38;
% y1 = awgn(sin(2*pi*fo*(linspace(0,M,M*fs))),-15,'measured','db')';
% y2 = awgn(sin(2*pi*fo*(linspace(0,M,M*fs))),-16,'measured','db')';
% Rfo = fo-(L/2):1/M:fo+(L/2);
% MLphi=zeros(size(Rfo));
% for i=1:length(Rfo)
% [MLphi(i),MLphicrit] = MLFT([y1 y2],L,fs,0.05,Rfo(i));
% end
% figure;plot(Rfo,MLphi,'b',fo,1.05*MLphi(fix(length(Rfo)/2)+1),'-rv',[Rfo(1) Rfo(end)],[MLphicrit MLphicrit],'k--')
% xlabel('Frequency (Hz)');ylabel('MLFT')
if bitget(L,1) %odd
    clc
    L = L+1;
    disp('L must be even. Now L=L+1');
end

[t,n] = size(y);   % length and number of channel
nfft = fix(t/2);   % Number of points for FFT
pfo = round(fo*nfft/(fs/2))+1; % Position of fo in FFT spectrum

% Amplitude spectrum in fo and closest neighbouring frequencies
Y = fft(y);
Y = abs(Y(1:nfft+1,:)); 
Yfo = Y(pfo,:); 
Yfn = Y(pfo-L/2:pfo+L/2,:);
Yfn(L/2+1,:) = []; 

% Compute F value
MLphi = sum(Yfo.^2)/(sum(1/L*sum(Yfn.^2,1)));

% Critical value
    MLphicrit = finv((1-alpha),2*n,2*n*L);
end
