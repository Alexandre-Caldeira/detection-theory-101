function [J,grad] = funcao_custo(alfa, MM, det, FP_desejado)

nRuns = size(det,1);
%1-Crit�rio de parada ------------------------
valor_critico = VC_MSC(MM,alfa);
det = det>repmat(valor_critico,nRuns,1);
FP = mean((sum(det,2)>0)); %TAXA DE FALSO POSITIVO DA SE��O (EXAME)

%2 - C�lculo da fun��o de custo ---------------------------
erro = (FP-FP_desejado);
%J = erro.^2;
J = erro.^2+100*(alfa<0);


%3 - C�lculo do gradiente -----------------------------
%grad = erro*FP;
grad = erro*FP; %~porporcional ao gradiente (n�o � o verdadeiro - mas porporcional)
 
