%função obter o alfa para todos os Mmax ------
%------------------------------------------------------------
%%function funcao_alfaCorrigido_Mmax
% Aplicação do protocolo de detecção na condição de Ho -------------
function [alfa_corrigido,cost_alfa, P] =funcao_alfaCorrigido_Mmax(nRuns,Mmax,FP_desejado)
%parametros de entrada ----------------

% nRuns  = 10000;
% Mmax = 30; %número máximo de janela 
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

for ii = 1:size(P,1)
    Mmin = P(ii,1);
    Mstep = P(ii,2); 
    Mmax = P(ii,3);
    MM = Mmin:Mstep:Mmax;
    
    det = ord(:,MM);
    alfa = 0.05;  %TAXA DE FALSO POSITIVO DE CADA TESTES
    options = optimset('MaxIter', 50);
    cc = @(alfa) funcao_custo(alfa ,MM, det, FP_desejado);                               
    [alfa, cost] = fmincg(cc,alfa, options);
    alfa_corrigido(ii) = alfa; 
    
    if ~isempty(cost)
        cost_alfa(ii) = cost(end);
    end
    
end

% figure 
% plot(alfa_corrigido,'ok')
% figure
% plot(cost_alfa,'ok')
%ylim([min(cost_alfa) max(cost_alfa)])





