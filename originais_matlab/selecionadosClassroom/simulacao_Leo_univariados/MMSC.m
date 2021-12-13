function [k2N,F,k2Ncrit] = MMSC(y,tj,fs,alpha)
%
%Multiple magnitude-squared coherence (MMSC)  
%
%MORD technique that makes use of both magnitude and phase and can be used
%as a detector of hidden periodicities in noise, since N real-valued 
%signals are available. If only one signal is used, than it reduces to the
%standart MSC. The method was first described in "Miranda de Sa, AMFL, Felix, LB 
%and Infantosi, AFC (2004). A Matrix-based Algorithm for Estimating 
%Multiple Coherence of a Periodic Signal and its Application to the 
%Multi-channel EEG during Sensory Stimulation. IEEE Transactions on 
%Biomedical Engineering. 51(7):1140-6." The implementation seen in this 
%code was done by refering to "Netto AD., Infantosi AFC, Miranda de Sá 
%AMFL(2015). A Sweep Operator-Based Algorithm for Multiple Coherence 
%Estimation in BCI Applications. In: Lackovi? I., Vasic D. (eds) 6th 
%European Conference of the International Federation for Medical and 
%Biological Engineering. IFMBE Proceedings, vol 45. Springer, Cham 
%
%Sintaxes:
%
%k2N = MMSC(y,tj) => MMSC spectrum.
%[k2N,F] = MMSC(y,tj,fs) => MMSC spectrum plus frequency vector. 
%[k2N,F,k2Ncrit] = MMSC(y,tj,fs,alpha) =>  returns also the theoretical 
%critical value.
%
%Input Parameters:
%
%y => matrix whose columns are the signals.  
%tj => number of points of each epoch in which the signals will be divided.
%fs => sample rate of signals.
%alpha => significance level, e.g. alpha = 0.05.
%
%Example: MMSC using two signals 
%
%y1 = awgn(sin(2*pi*21*(linspace(0,10,1000))),-15,'measured','db')'; 
%y2 = awgn(sin(2*pi*21*(linspace(0,10,1000))),-16,'measured','db')';
%[K2N,F,K2Ncrit] = MMSC([y1 y2],100,100,0.05);
%figure;plot(F,K2N,'b',21,1.05*K2N(22),'-rv', [0 50],[K2Ncrit K2Ncrit],'k--')
%axis([0 F(end) 0 1]);xlabel('Frequency (Hz)');ylabel('MMSC')

%Bonato mar/2019

nf = fix(tj/2)+1;
N = size(y,2); %numero de sinais
M = fix(size(y,1)/(tj));   %determina numero de segmentos
Sfft = zeros(tj,N,M);   %acumulador para fft em anel
k2N = zeros(nf,1);  %calculo da coerencia
idl = 1;

%monta matriz de espectros
for k=1:M
    
    %retira tendência linear
    %Sfft(:,:,k) = fft(detrend(canais(idl:(idl+L-1),:),'linear'));
    Sfft(:,:,k) = fft(y(idl:(idl+tj-1),:));
    
    idl = idl+tj;
end

for kf =1:nf
    %reinicia matrizes para calculo
    %V = zeros(N,1);
    %Specm = zeros(N,N);
    
    %matriz de espectros aumentados, parte com V e VH
    Specm_a = zeros(N+1,N+1);  
    for p = 1:N
        for ks = 1: M
            %V(p) = V(p) + Sfft(kf,p,ks);
            Specm_a(N+1,p) = Specm_a(N+1,p) + Sfft(kf,p,ks);
        end
        %VH - hermitiano de V
        Specm_a(p,N+1) = conj(Specm_a(N+1,p));
    end
    
    
    %monta mtriz de espectro aumentado
    for p = 1:N
        for q = 1:p
            for ks=1:M
                Specm_a(p,q) = Specm_a(p,q) + conj(Sfft(kf,p,ks)).*Sfft(kf,q,ks);
            end
            Specm_a(q,p) = conj(Specm_a(p,q));
        end
    end
    Specm_a(N+1,N+1) = 1;
    Specm_as = Msweep(Specm_a,N);
    k2N(kf) = (1-real(Specm_as(N+1,N+1)))/M;
end

if nargin > 2
    F = (0:(nf-1))*fs/tj; %frequency vector
    if nargin == 4
        Fcrit = finv(1-alpha,2*N,2*M-2*N);% F CDF
        k2Ncrit = Fcrit/(((M-N)/N)+Fcrit); %critical value
    end
end   


%%%%%%%%%%%%%%%%%%%%%%%%%%%


function M = Msweep(M,r)
%
%operador sweep para matrizes quadradas
%
ordem = size(M,1);
N = ordem;

for k=1:r

    %para m,n ~= r
    for m=1:N %linha
        for n=1:N %coluna
            if  (m ~= k && n~=k)
                M(m,n)=M(m,n)-M(m,k)*M(k,n)/(M(k,k)+eps);
            end
        end
    end

    %para (m ~= r && n == r)|| (n ~= r && m == r)
    for m=1:N
        for n=1:N
            if  (n ~= k && m==k)
                M(k,n)=M(k,n)/(M(k,k)+eps);
            end

            if  (m ~= k && n==k)
                M(m,k)=-M(m,k)/(M(k,k)+eps);
            end
        end
    end

    %para rr - > pivo
    M(k,k) = 1/(M(k,k)+eps);
end
