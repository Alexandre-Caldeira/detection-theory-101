%Fun��o que avalia o NDC e a taxa de FP da sess�o em fun��o do NDC 


%Paramentos da estrat�gia 
Mmax = 30;
Mmin = 10;
Mstep = 1; 
alfa = 0.05;

%%
nRuns = 10000;
tj = 32; %cada janela um segundo 
bin = 8; 

Ntotal = Mmax*tj; %n�mero de pontos totais 

%Na simula��o iremos estimar a aplica��o do detector a cada janela
ord = zeros(nRuns,Mmax); %armazena os valores dos detectores a cada experimento.

for ii = 1: nRuns    
    x = randn(Ntotal,1); 
    x = reshape(x,tj,Mmax); %dividir em janelas 
    %aplicar o detector a cada janela ------------------
    xfft = fft(x); %aplico uma ��nica vez a FFT.  
    for M = 2:Mmax %fazer para cada acrescimo de uma janela       
        ord(ii,M) = msc_fft(xfft(bin,1:M),M);        
    end   
end


%% 
%aplica��o do Estrat�gia para cada uma das simula��o 
MM = Mmin:Mstep:Mmax;
MM = MM';

dr = zeros(nRuns,1); 
vNDC = 1:20;
for NDC = vNDC 
    
    for ii = 1:nRuns
    
        dr(ii,1) = ETS(ord(ii,:),MM,alfa,NDC);
    
    end
    FP(NDC) = mean(dr); 
    dr = dr*0; 
end

figure 
plot(vNDC,FP(vNDC)*100,'-k','LineWidth',1.5)
hold on
plot([vNDC(1) vNDC(end)], [alfa alfa]*100,'-r','LineWidth',1.5)
xlabel('NDC','fontsize',12)
ylabel('FP (%) da sess�o','fontsize',12)



%% Obten��o do NDC "�timo"  ---PROBLEMA ---------
% FP_desejado = 0.05;
% 
% NDC = 2;  %TAXA DE FALSO POSITIVO DE CADA TESTES
% options = optimset('MaxIter', 50);
% cc = @(NDC) funcao_custo_v2(alfa, NDC, MM, ord, FP_desejado);                             
% [NDC, cost] = fmincg(cc,NDC, options);



%% Obten��o do NDC ---
Mmin = 10; 
Mstep = 1;
Mmax = 30;
alfa = 0.05;
FP_desejado = 0.05; 
Ninicial = 1;
[NDC,FP]  = estimarNDC(Ninicial,alfa,FP_desejado, ord, Mmin, Mstep, Mmax);


%ajustar os valores cr�tico
MM = Mmin:Mstep:Mmax;
options = optimset('MaxIter', 50);
cc = @(alfa) funcao_custo_v2(alfa, NDC, MM, ord, FP_desejado);                               
[alfa, cost] = fmincg(cc,alfa, options);
alfa_corrigido = alfa; 


%calcular o falso positivo
dr = zeros(nRuns,1);
FP =0; 
for ii = 1:nRuns    
        dr(ii,1) = ETS(ord(ii,:),MM,alfa_corrigido,NDC); 
end
FP = mean(dr); 












