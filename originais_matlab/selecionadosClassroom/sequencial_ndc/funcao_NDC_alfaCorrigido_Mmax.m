%função obter o alfa para todos os Mmax ------
%------------------------------------------------------------
%%function funcao_alfaCorrigido_Mmax
% Aplicação do protocolo de detecção na condição de Ho -------------
function [alfa_corrigido,NDC_minimo,cost_alfa, P] =funcao_NDC_alfaCorrigido_Mmax(nRuns,Mmax,alfa_teste,FP_desejado)
%parametros de entrada ----------------
% clear all, close all, clc
% nRuns  = 10000;
% Mmax = 20; %número máximo de janela 
% alfa_teste = 0.05;
% FP_desejado =0.05;
%----------------------------------------

%parâmetros defaul 
% fs = 64; 
tj = 32; %cada janela um segundo 
bin = 8; 


Ntotal = Mmax*tj; %número de pontos totais 

%Na simulação iremos estimar a aplicação do detector a cada janela
ord = zeros(nRuns,Mmax); %armazena os valores dos detectores a cada experimento.

for ii = 1: nRuns    
    x = randn(Ntotal,1); 
    x = reshape(x,tj,Mmax); %dividir em janelas 
    %aplicar o detector a cada janela ------------------
    xfft = fft(x); %aplico uma ´única vez a FFT.  
    for M = 2:Mmax %fazer para cada acrescimo de uma janela       
        ord(ii,M) = msc_fft(xfft(bin,1:M),M);        
    end   
end


%% obter todos os parâmetros 
P = parametros_protocolo(Mmax);
alfa_corrigido = nan*ones(size(P,1),1);
cost_alfa = nan*ones(size(P,1),1);

NDC_minimo = nan*ones(size(P,1),1); %NDC_mínimo

%% 
[~,aux] = sort(ceil((Mmax - P(:,1))./P(:,2))+1);
P = P(aux,:);


for ii = 1:size(P,1)
    Mmin = P(ii,1);
    Mstep = P(ii,2); 
    Mmax = P(ii,3);
    MM = Mmin:Mstep:Mmax;
    
    
    %NDC:sugerir um bom 
    if ii == 1
        Ninicial = 1;
    else
        Ninicial =  max(round(NDC_minimo(ii-1))-3,1);
    end
    
    
    [NDC,FP]  = estimarNDC(Ninicial,alfa_teste,FP_desejado, ord, Mmin, Mstep, Mmax);
    NDC_minimo(ii) = NDC;
    
    %ajustar os valores crítico
    MM = Mmin:Mstep:Mmax;
    options = optimset('MaxIter', 50);
    cc = @(alfa) funcao_custo_v2(alfa, NDC_minimo(ii), MM, ord, FP_desejado);                               
    [alfa, cost] = fmincg(cc,alfa_teste, options);
    alfa_corrigido(ii) = alfa;
    
%     if isempty(cost)       
%         y= 1;        
%     end
        
    if ~isempty(cost)
        cost_alfa(ii) = cost(end);
    end
    

end

% figure 
% plot(alfa_corrigido,'ok')
% figure
% plot(cost_alfa,'ok')
%ylim([min(cost_alfa) max(cost_alfa)])





