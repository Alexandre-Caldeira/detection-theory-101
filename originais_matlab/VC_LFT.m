function [VC_MC,VC_teorico] = VC_LFT(L,alfa,nRuns)

%valor crítico 
%parâmetros 

% L = 12; %tamanho L
% alfa = 0.05 %nível de significância 
% nRuns = 10000; 

%parametros para a simulação ----------------------
N = 256; %tamanho da janela
bin = 50; 

ORD = zeros(nRuns,1); 

for ii = 1:nRuns
    x = randn(N,1);
    ORD(ii) = LFT(x,L,bin);
end

VC_MC  = quantile(ORD,1-alfa);

%valor teorico 
VC_teorico = finv(1-alfa,2,2*L);



