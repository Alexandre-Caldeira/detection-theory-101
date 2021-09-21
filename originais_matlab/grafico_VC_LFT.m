%Gráfico LFT 

clear all, close all, clc



%plotar valores críticos
alfa = 0.05;
nRuns = 1e3; 

vL = 2:2:50; 
Nl = size(vL,2);

VC_MC = zeros(Nl,1); 
VC_teorico = zeros(Nl,1);


for ii = 1:Nl 
    [VC_MC(ii),VC_teorico(ii)] = VC_LFT(vL(ii),alfa,nRuns);
end


%%
figure 
subplot(211)
lw = 1.5;
plot(vL,VC_MC,'-k','LineWidth',lw)
hold on 
plot(vL,VC_teorico,':k','LineWidth',lw)
hold off 
grid on
xlabel('L','fontsize',14)
ylabel('Valor Crítico','fontsize',14)
legend({'MC';'Teorico'},'edgecolor','none','fontsize',14)
title('Análise de convergência: VC LFT','fontsize',14)

subplot(212)
plot(vL,abs(VC_MC-VC_teorico),'-r','LineWidth',lw)
ylabel('Erro absoluto','fontsize',14)
grid on
xlabel('Tamanho das bandas laterais [L]','fontsize',14)

