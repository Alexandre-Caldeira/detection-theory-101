clear all
close all
clc

nRuns  = 10000;
Mmax = 50; %n�mero m�ximo de janela 
alfa = 0.05;
FP_desejado =0.05;


[alfa_corrigido,NDC_minimo,cost_alfa, P] = funcao_NDC_alfaCorrigido_Mmax(nRuns,Mmax ,alfa,FP_desejado);
%salvar as vari�veis 
save(['NDC_AlfaCorrigido_Mmax' num2str(Mmax) '_alfa_'  num2str(alfa) '_FPdesejado' num2str(FP_desejado) '.mat'],'alfa_corrigido', ...
    'NDC_minimo','P', 'nRuns','alfa','FP_desejado')


t1 = table(P, NDC_minimo,alfa_corrigido,cost_alfa);

%% Avaliar a varian�a na estima��o do NDC e alfa 

%repetir isso v�rias vezes 
% [alfa_corrigido,NDC_minimo,cost_alfa, P] = funcao_NDC_alfaCorrigido_Mmax(nRuns,Mmax ,alfa_teste,FP_desejado);
% comparar os valor estimados e c�lcular a vari�ncia

%% 



