%% Tra�ando graficos
close all;clearvars;clc;
% Carrega variavies (t�m que estar na pasta local do caminho)
load('Workspace_M1','FP');
fp1 = FP;
load('Workspace_M2','FP');
fp2 = FP;
load('Workspace_M3','FP');
fp3 = FP;
load('Workspace_M4','FP');
fp4 = FP;
load('Workspace_M1','linf');
load('Workspace_M1','lsup');
%%
% cria legendas
leg0 = 'FP Desejado = 5%';
legc1 = ['Confian�a Inf. = ',num2str(linf),'%'];
legc2 = ['Confian�a Sup. = ',num2str(lsup),'%'];
legC = ['Confian�a Inf. = [',num2str(linf),',',num2str(lsup),'] %'];
leg1 = ['M�todo 1, m�dia = ',num2str(round(mean(fp1)*100,2)),'%'];
leg2 = ['M�todo 2, m�dia = ',num2str(round(mean(fp2)*100,2)),'%'];
leg3 = ['M�todo 3, m�dia = ',num2str(round(mean(fp3)*100,2)),'%'];
leg4 = ['M�todo 4, m�dia = ',num2str(round(mean(fp4)*100,2)),'%'];
tamanho = 20;


%% plota 
fig = figure(1);
x = [0 1400];
y = 0.05*100*ones(1,length(x));
ysup = lsup*ones(1,length(x));
yinf = linf*ones(1,length(x));

subplot(221)
title('M�todo 1','FontName','Times New Roman','FontSize',16);
scatter(1:size(fp1,1),fp1*100,tamanho,'filled','CData',[0, 0.4470, 0.7410]);
grid on
hold on
plot(x,yinf, 'Color','k','LineStyle','--','LineWidth',2);
plot(x,ysup, 'Color','k','LineStyle','--','LineWidth',2);


subplot(222)
title('M�todo 2','FontName','Times New Roman','FontSize',16);
scatter(1:size(fp2,1),fp2*100,tamanho,'filled','CData',[0, 0.4470, 0.7410]);
grid on
hold on
plot(x,yinf, 'Color','k','LineStyle','--','LineWidth',2);
plot(x,ysup, 'Color','k','LineStyle','--','LineWidth',2);
legend(leg0,legc1,legc2,'Times New Roman','FontSize',10)

subplot(223)
title('M�todo 3','FontName','Times New Roman','FontSize',16);
scatter(1:size(fp3,1),fp3*100,tamanho,'filled','CData',[0, 0.4470, 0.7410]);
grid on
hold on
plot(x,yinf, 'Color','k','LineStyle','--','LineWidth',2);
plot(x,ysup, 'Color','k','LineStyle','--','LineWidth',2);

subplot(224)
title('M�todo 4','FontName','Times New Roman','FontSize',16);
scatter(1:size(fp4,1),fp4*100,tamanho,'filled','CData',[0, 0.4470, 0.7410]);
grid on
hold on
plot(x,yinf, 'Color','k','LineStyle','--','LineWidth',2);
plot(x,ysup, 'Color','k','LineStyle','--','LineWidth',2);


han=axes(fig,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
a = ylabel(han,'Taxa de Falso Positivo [%]','FontName','Times New Roman','FontSize',20);
xlabel(han,'�ndice do Conjunto de Par�metros','FontName','Times New Roman','FontSize',20);
