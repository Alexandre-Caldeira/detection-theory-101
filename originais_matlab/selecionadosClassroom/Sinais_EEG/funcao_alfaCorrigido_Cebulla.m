%% Correção de Cebulla, 2005
% Alpha (taxa de falsos positivos) é corrigido em função do número de "test 
% steps", ou seja, número de testes a serem aplicados em sequencia.

function [alfa_corrigido, P] =funcao_alfaCorrigido_Cebulla(nRuns,Mmax,FP_desejado)
%parâmetros default
% fs = 64; 
tj = 32; %cada janela um segundo 
Ntotal = Mmax*tj; %número de pontos totais 
bin = 8;

%Na simulação iremos estimar a aplicação do detector a cada janela
ord = zeros(nRuns,Mmax); %armazena os valores dos detectores a cada experimento.
maximo = zeros(nRuns,Mmax);

for ii = 1: nRuns    
    x = randn(Ntotal,1); 
    x = reshape(x,tj,Mmax); %dividir em janelas 
    
    %aplicar o detector a cada janela ------------------
    xfft = fft(x); %aplico uma ´única vez a FFT. 

    for M = 2:Mmax %fazer para cada acrescimo de uma janela       
%         ord(ii,M) = msc_fft(xfft(bin,1:M),M);  

%         ord(ii,M) = CSM_fft(xfft,M);
        ord(ii,M) = CSM_fft(xfft(bin,1:M),M);
        
    end   
end

%% obter todos os parâmetros 
P = parametros_protocolo(Mmax);
alfa_corrigido = nan*ones(size(P,1),1); % pegar o valor maximo

for ii = 1:size(P,1)
    Mmin = P(ii,1);
    Mstep = P(ii,2); 
    Mmax = P(ii,3);
    MM = Mmin:Mstep:Mmax;
       
    det = ord(:,MM);
    nTests = numel(MM);   
    
    maior = max(det);
    if maior > maximo(M)
        maximo(ii,M) = maior;
    end
    
    
    alfa_corrigido(ii) = maximo_med(Mmax);
    
end



end
