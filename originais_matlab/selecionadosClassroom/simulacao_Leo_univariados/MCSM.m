function [csmN,F,csmNcrit] = MCSM(y,tj,fs,alpha)
%
%Multiple componenent synchrony measure (MCSM)  
%
%MORD technique that makes use only phase information and can be used as a
%detector of hidden periodicities
%in noise, since N real-valued signals are available. The method is 
%described in "Miranda de Sa, AMFL and Felix, LB (2003). Multi-channel 
%evoked response detection using only phase information. Journal of 
%Neuroscience Methods. 129(1):1-10."
%
%Sintaxes:
%
%[csmN] = MCSM(y,tj) => MMSC spectrum.
%[csmN,F] = MCSM(y,tj,fs) => MMSC spectrum plus frequency vector. 
%[csmN,F,csmNcrit] = MCSM(y,tj,fs,alpha) =>  returns also the theoretical 
%critical value.
%
%Input Parameters:
%
%y => matrix whose columns are the signals to be considered. 
%tj => number of points of each epoch.
%fs => sample rate.
%alpha => significance level, e.g. alpha = 0.05. 
%
%Example: MCSM using two signals
%
%y1 = awgn(sin(2*pi*21*(linspace(0,10,1000))),-15,'measured','db')';
%y2 = awgn(sin(2*pi*21*(linspace(0,10,1000))),-16,'measured','db')';
%[csmN,F,csmNcrit] = MCSM([y1 y2],100,100,0.05);
%figure;plot(F,csmN,'b',21,1.05*csmN(22),'-rv', [0 50],[csmNcrit csmNcrit],'k--')
%axis([0 F(end) 0 1]);xlabel('Frequency (Hz)');ylabel('MCSM')

%Leonardo Bonato Felix - Feb/2019



[tamsinal,N] = size(y);
nfft = fix(tj/2); %number of points in FFT
M = fix(tamsinal/tj); %number of windows  
y = y(1:M*tj,:); %limit the signal for a entire number of windows. Prevents shorter windows

for i = 1:N
    Y(:,:,i)=fft(reshape(y(:,i),tj,M)); %3D matrix ith-slice is the FFT of yi[n], i=1,2,...N 
end
Y = Y(1:nfft+1,:,:); %only half of the FFT-values are returned, because of FFT's simmetry in real signals

%Algorithm's start 
teta = angle(Y);
C = cos(teta);
S = sin(teta);
Cmed = mean(C,3);
Smed = mean(S,3);
Cmed(1,:)=NaN; Cmed(end,:)=NaN ; Smed(1,:)=NaN; Smed(end,:)=NaN; %desconsider the first and last freq. values

%temp1 = zeros(size(Cmed));
%temp2 = temp1;
temp1 = atan( (Smed.*(Cmed<0))./Cmed)+pi*(Cmed<0); %calculates mean teta matrix only for Cmed<0, zero-padding other positions
temp2 = atan( (Smed.*(Cmed>=0))./Cmed); %calculates mean teta matrix only for Cmed>0, zero-padding other positions
teta_med = temp1 + temp2; 

csmN=(1/M^2)*sum(cos(teta_med')).^2+(1/M^2)*sum(sin(teta_med')).^2; %multiple CSM
%End 

if nargin > 2
    F = 0:fs/tj:fs/2; %frequency vector
    if nargin == 4
        csmNcrit = chi2inv(1-alpha,2)/(2*M); %associated critical value
    end
end
end

