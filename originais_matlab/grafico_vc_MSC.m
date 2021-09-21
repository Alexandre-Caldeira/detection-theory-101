%Gráfico LFT 
clearvars; close all; clc;

%plotar valores críticos
alfa = 0.05;
nRuns = 1e3; 

vM = 10:5:400; 
nM = size(vM,2);

VC_MC = zeros(nM,1); 
VC_teorico = zeros(nM,1);


for ii = 1:nM 
    [VC_MC(ii),VC_teorico(ii)] = VC_MSC(vM(ii),alfa,nRuns);
end

%%
figure 
subplot(211)
lw = 1.5;
plot(vM,VC_MC,'-k','LineWidth',lw)
hold on 
plot(vM,VC_teorico,':k','LineWidth',lw)
hold off 
grid on
xlabel('M','fontsize',14)
ylabel('Valor Crítico','fontsize',14)
legend({'MC';'Teorico'},'edgecolor','none','fontsize',14)
title('Análise de convergência: VC MSC','fontsize',14)

subplot(212)
plot(vM,abs(VC_MC-VC_teorico),'-r','LineWidth',lw)
ylabel('Erro absoluto','fontsize',14)
grid on
xlabel('Número de janelas [M]')
