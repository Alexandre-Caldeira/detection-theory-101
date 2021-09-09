%Gráfico LFT 

clear all, close all, clc



%plotar valores críticos
alfa = 0.05;
nRuns = 30000; 

vL = 2:2:20; 
Nl = size(vL,2);

VC_MC = zeros(Nl,1); 
VC_teorico = zeros(Nl,1);


for ii = 1:Nl 
    [VC_MC(ii),VC_teorico(ii)] = VC_LFT(vL(ii),alfa,nRuns);
end


figure 
lw = 1.5;
plot(vL,VC_MC,'-k','LineWidth',lw)
hold on 
plot(vL,VC_teorico,':k','LineWidth',lw)
hold off 
xlabel('L','fontsize',12)
ylabel('Valor Crítico','fontsize',12)
legend({'MC';'Teorico'},'edgecolor','none','fontsize',12)

