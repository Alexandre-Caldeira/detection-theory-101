function [VC_MC,VC_teorico] = VC_CSM(M,alfa,nRuns)


%parametros para a simulação ----------------------
N = 32; %tamanho da janela
bin = 7;

ORD = zeros(nRuns,1); 

for ii = 1:nRuns
    x = randn(N*M,1);
    aux = CSM(x,N,M);
    ORD(ii) = aux(bin); 
end

VC_MC  = quantile(ORD,1-alfa);

%valor teorico 
VC_teorico = chi2inv(1-alfa,2)/(2*M);



