%gerar a curva FA para os detectores 
clear all
% close all
clc

% c = parcluster('local'); % build the 'local' cluster object
% nw = c.NumWorkers;        % get the number of workers
% parpool('local',nw,'IdleTimeout', Inf);
% 

M =30;alpha=0.05;Nruns=1e6;
L = 12;
N =1;
vetor_SNR = [-30:0.5:5]; 
vetor_SNR = [-5:0.5:15]; 




% %% Simulação no tempo --------------------------
%Nmax =2; 
%for N =1:Nmax
N =1; 
NFFT = 32; %lowest base two with good results. 
freq = round(NFFT*.1);
s = cos(2*pi*freq*[1:(M*NFFT)]/(NFFT));
s = repmat(s,N,1);
sigma_n = 2/NFFT;

   
Ptempo_MSC = zeros(Nruns,length(vetor_SNR));
Ptempo_CSM = zeros(Nruns,length(vetor_SNR));
Ptempo_LFT = zeros(Nruns,length(vetor_SNR));

%valor críticos 
MMSCcrit = betainv(1-alpha,N,M-N);
MCSMcrit = chi2inv(1-alpha,2)/(2*M);
MLFTcrit = finv((1-alpha),2*N,2*N*L);
% 
% 
% for jj =1:size(vetor_SNR,2) %quicker
%     
%         SNR = 10^(vetor_SNR(jj)/10);
%         %parfor
%         for ii = 1:Nruns %número de repetições
%               
%             %Mesmo Sinal para os três detectores 
%             y = sqrt(4*sigma_n*SNR/NFFT)*s + sqrt(sigma_n).*randn(size(s));       
%             
%             ORD_MSC = MMSC(y',NFFT);
%             ORD_CSM = MCSM(y',NFFT);
%             ORD_LFT = MLFT(y',L,NFFT,alpha,freq);
%             
%             Ptempo_MSC(ii,jj) = ORD_MSC(freq+1)>MMSCcrit;
%             Ptempo_CSM(ii,jj) = ORD_CSM(freq+1)>MCSMcrit;
%             Ptempo_LFT(ii,jj) = ORD_LFT>MLFTcrit;
%                                 
%         end         
% end
% 
% PD_MSCt = mean(Ptempo_MSC);
% PD_CSMt = mean(Ptempo_CSM);
% PD_LFTt = mean(Ptempo_LFT);
%%
% figure
% lw = 1.5;
% mz = 1;
% hold on 
% plot(vetor_SNR, PD_MSCt,'k','LineWidth',lw,'MarkerSize',mz);
% plot(vetor_SNR, PD_LFTt,'r','LineWidth',lw,'MarkerSize',mz);
% plot(vetor_SNR, PD_CSMt,'b','LineWidth',lw,'MarkerSize',mz);
% legend({'MSC','LFT','CSM'},'fontsize',12,'EdgeColor','none','Location','best')
% ylabel('PD','fontsize',12); 
% xlabel('SNR','fontsize',12); 
% set(gca,'XMinorTick','on')
% set(gca,'YMinorTick','on')
% hold off 
% box off
% % xlim([-30 15])
% title(['Tempo N=' num2str(N) 'M=' num2str(M) '- Mesmo sinal para todos os detectores'])



%% simulação na Frequência 


Pfreq_MSC = zeros(Nruns,length(vetor_SNR));
Pfreq_CSM = zeros(Nruns,length(vetor_SNR));
Pfreq_LFT = zeros(Nruns,length(vetor_SNR));
for jj =1:size(vetor_SNR,2) %quicker
    
            SNR = 10^(vetor_SNR(jj)/10);
                    
            Y = sqrt(2*SNR) +(randn(M,Nruns)+j*randn(M,Nruns));
            %MSC 
            ORD_MSC = abs(sum(Y).^2)./(M*sum(abs(Y).^2));
            
            
            %CSM 
            teta = angle(Y); 
            ORD_CSM =   (1/M^2)*sum(cos(teta)).^2+(1/M^2)*sum(sin(teta)).^2;
            
            
            %LFT 
            %1 simulação: Igual ao artigo 2016 
            Y = sqrt(2*SNR) +(randn(1,Nruns)+j*randn(1,Nruns));
            ruido = (randn(L,Nruns)+j*randn(L,Nruns));
            
            ORD_LFT1 = abs(Y).^2./(1/L*sum(abs(ruido).^2));                      
                        
            %LFT 
            %2 simulação: Minha alteração para que na simulação no domínio da freq. 
            %a SNR seja equivalente ao domínio do tempo 
            %multiplica por M 
            Y = sqrt(2*SNR*M) +(randn(1,Nruns)+j*randn(1,Nruns));
            ruido = (randn(L,Nruns)+j*randn(L,Nruns));
            
            ORD_LFT2 = abs(Y).^2./(1/L*sum(abs(ruido).^2));

            Pfreq_MSC(:,jj) = ORD_MSC>MMSCcrit;            
            Pfreq_CSM(:,jj) = ORD_CSM>MCSMcrit;
            Pfreq_LFT1(:,jj) = ORD_LFT1>MLFTcrit;
            Pfreq_LFT2(:,jj) = ORD_LFT2>MLFTcrit;                                      
end


PD_MSCf = mean(Pfreq_MSC);
PD_CSMf = mean(Pfreq_CSM);
PD_LFT1f = mean(Pfreq_LFT1);
PD_LFT2f = mean(Pfreq_LFT2);


%% plotar na frequência 
figure
lw = 1.5;
mz = 1;
hold on 
plot(vetor_SNR, PD_MSCf,'k','LineWidth',lw,'MarkerSize',mz);
% plot(vetor_SNR, PD_LFT1f,'g','LineWidth',lw,'MarkerSize',mz);
plot(vetor_SNR, PD_LFT2f,'r','LineWidth',lw,'MarkerSize',mz);
plot(vetor_SNR, PD_CSMf,'b','LineWidth',lw,'MarkerSize',mz);

% legend({'MSC','LFT (2106)','LFT (Mod)','CSM'},'fontsize',12,'EdgeColor','none','Location','best')
legend({'MSC','LFT ','CSM'},'fontsize',12,'EdgeColor','none','Location','best')

ylabel('PD','fontsize',12); 
xlabel('SNR','fontsize',12); 
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
hold off 
box off
% xlim([-30 15])
title(['Frequência - N=' num2str(N) ' M=' num2str(M)])


%% Diferença - Tempo - Frew 
% figure
% lw = 1.5;
% mz = 1;
% hold on 
% plot(vetor_SNR, PD_MSCt - PD_MSCf,'k','LineWidth',lw,'MarkerSize',mz);
% plot(vetor_SNR, PD_LFTt-PD_LFT1f,'g','LineWidth',lw,'MarkerSize',mz);
% plot(vetor_SNR, PD_LFTt-PD_LFT2f,'r','LineWidth',lw,'MarkerSize',mz);
% plot(vetor_SNR, PD_CSMt -PD_CSMf,'b','LineWidth',lw,'MarkerSize',mz);
% legend({'MSC','LFT (2106)','LFT (Mod)','CSM'},'fontsize',12,'EdgeColor','none','Location','best')
% ylabel('PD','fontsize',12); 
% xlabel('SNR','fontsize',12); 
% set(gca,'XMinorTick','on')
% set(gca,'YMinorTick','on')
% hold off 
% box off
% % xlim([-30 15])
% title(['Tempo - Frequência - N=' num2str(N) ' M=' num2str(M)])
% 



