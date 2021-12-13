%Gráfico LFT 

clear all, close all, clc

%plotar valores críticos
alfa = 0.05;
nRuns = 30000; 

vM = 5:5:50; 
nM = size(vM,2);

VC_MC = zeros(nM,1); 
VC_teorico = zeros(nM,1);


for ii = 1:nM 
    [VC_MC(ii),VC_teorico(ii)] = VC_CSM(vM(ii),alfa,nRuns);
end


figure 
lw = 1.5;
plot(vM,VC_MC,'-k','LineWidth',lw)
hold on 
plot(vM,VC_teorico,':k','LineWidth',lw)
hold off 
xlabel('M','fontsize',12)
ylabel('Valor Crítico','fontsize',12)
legend({'MC';'Teorico'},'edgecolor','none','fontsize',12)

