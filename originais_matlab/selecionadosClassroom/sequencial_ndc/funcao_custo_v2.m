function [J,grad] = funcao_custo_v2(alfa, NDC, MM, ord, FP_desejado)

nRuns = size(ord,1);
alfa = max(alfa,eps);
alfa = min(alfa,1-eps); 
%1-Crit�rio de parada ------------------------
dr = zeros(nRuns,1); 
for ii = 1:nRuns   
   dr(ii,1) = ETS(ord(ii,:),MM,alfa,NDC);   
end
FP = mean(dr); 

%2 - C�lculo da fun��o de custo ---------------------------
erro = (FP-FP_desejado);
%J = erro.^2;
J = 1/2*erro.^2;


%3 - C�lculo do gradiente -----------------------------
grad = erro;
    
%grad = erro*FP; %~porporcional ao gradiente (n�o � o verdadeiro - mas porporcional)
 