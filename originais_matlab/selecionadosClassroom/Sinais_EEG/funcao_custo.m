function [J,grad] = funcao_custo(alfa, MM, det, FP_desejado)

nRuns = size(det,1);
%1-Critério de parada ------------------------
valor_critico = VC_MSC(MM,alfa);
det = det>repmat(valor_critico,nRuns,1);
FP = mean((sum(det,2)>0)); %TAXA DE FALSO POSITIVO DA SEÇÃO (EXAME)

%2 - Cálculo da função de custo ---------------------------
erro = (FP-FP_desejado);
%J = erro.^2;
J = erro.^2+100*(alfa<0);


%3 - Cálculo do gradiente -----------------------------
%grad = erro*FP;
grad = erro*FP; %~porporcional ao gradiente (não é o verdadeiro - mas porporcional)
 
