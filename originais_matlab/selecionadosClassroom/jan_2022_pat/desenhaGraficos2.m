%% Traçando graficos
close all;clearvars;clc;
% Carrega variavies (têm que estar na pasta local do caminho)
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
legc1 = ['Confiança = ',num2str(linf),'%'];
legc2 = ['Confiança = ',num2str(lsup),'%'];
legC = ['Confiança = [',num2str(linf),',',num2str(lsup),'] %'];
leg1 = ['Método 1, média = ',num2str(round(mean(fp1)*100,2)),'%'];
leg2 = ['Método 2, média = ',num2str(round(mean(fp2)*100,2)),'%'];
leg3 = ['Método 3, média = ',num2str(round(mean(fp3)*100,2)),'%'];
leg4 = ['Método 4, média = ',num2str(round(mean(fp4)*100,2)),'%'];
tamanho = 20;

%% plota 1
figure(1)
lsup = 6.73;
linf = 3.37;

x = [0 length(FP)];
y = 0.05*100*ones(1,length(x));
ym = mean(fp1)*100*ones(1,length(x));
ysup = lsup*ones(1,length(x));
yinf = linf*ones(1,length(x));

plot(x,y,'Color','b','LineStyle','--','LineWidth',2);
hold on

plot(x,ym,'Color','k','LineStyle',':','LineWidth',2);
plot(x,yinf, 'Color','k','LineStyle','--','LineWidth',2);
plot(x,ysup, 'Color','k','LineStyle','--','LineWidth',2);

scatter(1:size(fp1,1),fp1*100,tamanho,'filled','CData',[0, 0.4470, 0.7410]);

legend(leg0,leg1,legC)

ylabel('Falso Positivo','fontsize',12);
xlim([0,length(FP)])

%% plota metodo 1 E 2
figure(2)

plot(x,y,'Color','b','LineStyle','--','LineWidth',2);
hold on

plot(x,yinf, 'Color','k','LineStyle','--','LineWidth',2);
plot(x,ysup, 'Color','k','LineStyle','--','LineWidth',2);

scatter(1:size(fp1,1),fp1*100,tamanho,'filled','CData',[0, 0.4470, 0.7410]);
scatter(1:size(fp2,1),fp2*100,tamanho,'filled','CData',[0.8500, 0.3250, 0.0980]);

legend(leg0,legc1,legc2,leg1,leg2)
ylabel('Falso Positivo','fontsize',12);
xlim([0,length(FP)])

%% plota metodo 3 E 4
figure(3)

plot(x,y,'Color','b','LineStyle','--','LineWidth',2);
hold on

plot(x,yinf, 'Color','k','LineStyle','--','LineWidth',2);
plot(x,ysup, 'Color','k','LineStyle','--','LineWidth',2);

scatter(1:size(fp3,1),fp3*100,tamanho,'filled','CData',[0, 0.4470, 0.7410]);
scatter(1:size(fp4,1),fp4*100,tamanho,'filled','CData',[0.8500, 0.3250, 0.0980]);

legend(leg0,legc1,legc2,leg3,leg4)
ylabel('Falso Positivo','fontsize',12);
xlim([0,length(FP)])

%% plota histogramas separados
figure(4)
subplot(2,2,1)
histogram(fp1,20);
legend(leg1,'Location','northwest')
ylabel('Frequência Observada')


subplot(2,2,2)
histogram(fp2,15);
legend(leg2)

subplot(2,2,3)
histogram(fp3,20);
legend(leg3)
ylabel('Frequência Observada')
xlabel('Falso Positivo')

subplot(2,2,4)
histogram(fp4,15);
legend(leg4)
xlabel('Falso Positivo')

%% plota histogramas juntos
figure(5)
hold on
h1 = histogram(fp1,20);
h2 = histogram(fp2,15);
h3 = histogram(fp3,20);
h4 = histogram(fp4,15);
legend(leg1,leg2,leg3,leg4)

title('Distribuição dos FP')
ylabel('Frequência Observada')
xlabel('Falso Positivo')



