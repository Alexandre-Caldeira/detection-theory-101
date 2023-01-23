function [VC_teorico] = VC_MSC(M,alfa)

%function [VC_MC,VC_teorico] = VC_MSC(M,alfa,nRuns)


%parametros para a simulação ----------------------
% N = 32; %tamanho da janela
% bin = 7;
% 
% ORD = zeros(nRuns,1); 
% 
% for ii = 1:nRuns
%     x = randn(N*M,1);
%     aux = MSC(x,N,M);
%     ORD(ii) = aux(bin); 
% end
% 
% VC_MC  = quantile(ORD,1-alfa);

%valor teorico 
VC_teorico = 1 - alfa.^(1./(M-1));



