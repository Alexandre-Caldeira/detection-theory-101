% Código adaptado do material de aula de Thiago Zanotelli
clearvars; close all; clc ;

% numero de simulações:
NRuns = 1e4; 

%parametros do detector ------------------- 
energia_ruido = 1; 
alpha = 0.05;
f = 100; 
fs = 1000; 
tj = 100;
M = 300;
N=M*tj+1;
bin = round(tj*f/fs +1);
vetor_SNR =-60:2:-20; %vetor da SNR 

%% 
disp('-')
%obter o limiar 
[~,limiar] = VC_MSC(M,alpha,1);

for jj = 1:size(vetor_SNR,2) %fazer para cada SNR 
    nd = 0;
    SNRi = vetor_SNR(jj);
    A = 1/(10^(SNRi/-20));
    
    for ii = 1:NRuns %fazer por experimento 
        s(:,1) = sin(2*pi*f/fs*[0:(N-1)]);
        x = A*s + energia_ruido*randn(N,1); 
     
        ORD = msc(x,tj,M);

        if ORD(bin) > limiar
            nd = nd+1;
        end
    end %fim de cada experimento
    
    PD(jj) = nd/NRuns;
    disp(nd/NRuns)
    disp(vetor_SNR(jj))
    
end

%%
figure 
lw  =2; 
plot(vetor_SNR,PD,'LineWidth',lw) 
xlabel('SNR','fontsize',14)
ylabel('PD','fontsize',14)
title('Curva PD da MSC')
% legend({num2str(vN')},'edgecolor','none')
hold on
grid on 
pd_msc_python = readtable('../notebooks/PD_MSC.csv');
snr_p = pd_msc_python.SNR;
pd_p = pd_msc_python.ProbabilidadeDeDetec____o___;
plot(snr_p,pd_p,'-o') 
hold off
grid on
legend('MATLAB','PYTHON')

