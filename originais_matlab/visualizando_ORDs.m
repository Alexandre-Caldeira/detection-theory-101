% 09/Set/21 - Alexandre Caldeira
%% Boas pr�ticas
clearvars;close all; clc;

%% Constantes:
% Amplitude do Ru�do:
Ar = 5;

% M�dia do ru�do:
Mr = 10;

% Frequ�ncia do Sinal: [Hz]
f0 = 100;

% Fase inicial do sinal: [rad]
th0 = 0;

% Frequ�ncia de Amostragem: [Hz]
fs = 1000;

% Tamanho da janela:
tj = fs; % 1seg/janela

% N�mero de janelas: 
M = 10;

% N�mero de pontos amostrados:
N = M*tj; % 

% Tempo:
t = linspace(0,N-1,N);

%% Gera��o de sinais
% Distribui��o normal (gaussiana):
dist = randn(N,1);  

% "Ru�do" aleat�rio (m�dia 1, amplitude 5):
r = (Ar*randn(N,1) + Mr);

% Sen�ide:
s = 5*sin(2*pi*f0/fs*t+th0)';

% Sen�ide + Ru�do 1:
j = s+dist;

% Sen�ide + Ru�do 2:
k = s+r;

%% Gera��o das ORDs:
alfa = 0.05; L = 20;
[~,vc_lft] = VC_LFT(L,alfa, 1);
[~,vc_csm] = VC_CSM(M,alfa, 1);
[~,vc_msc] = VC_MSC(M,alfa, 1);

% Distribui��o normal (gaussiana):
csmD = CSM(dist,tj,M);  
mscD = msc(dist,tj,M);  


% "Ru�do" aleat�rio (m�dia 1, amplitude 5):
csmR = CSM(r,tj,M);  
mscR = msc(r,tj,M);  

% Sen�ide:
csmS = CSM(s,tj,M);  
mscS = msc(s,tj,M); 

% Sen�ide + Ru�do 1:
csmSR1 = CSM(j,tj,M);  
mscSR1 = msc(j,tj,M); 

% Sen�ide + Ru�do 2:
csmSR2 = CSM(k,tj,M);  
mscSR2 = msc(k,tj,M); 

%% Gr�fico das CSMs:
figure(2)
% subplot(211)
subplot(221)
hold on;
sMD = stem(csmD);
sMD.Color = [0.8500, 0.3250, 0.0980];
sMR = stem(csmR);
sMR.Color = [0, 0.4470, 0.7410];
yline(vc_csm,'--r','LineWidth', 2.5);
hold off;
grid on;
title('CSM Ruidos')
legend('Ru�do 1','Ru�do 2','Limiar')

% subplot(212)
subplot(222)
sMS = stem(csmS);
yline(vc_csm,'--r','LineWidth', 2.5);
grid on;
title('CSM Sen�ide')

% figure(2)
% subplot(211)
subplot(223)
sMSR1 = stem(csmSR1);
yline(vc_csm,'--r','LineWidth', 2.5);
grid on;
title('CSM Sen�ide + Ru�do 1')


% subplot(212)
subplot(224)
sMSR2 = stem(csmSR2);
yline(vc_csm,'--r','LineWidth', 2.5);
grid on;
title('CSM Sen�ide + Ru�do 2')

%% Gr�fico das MSCs:
figure(3)
% subplot(211)
subplot(221)
hold on;
sMD = stem(mscD);
sMD.Color = [0.8500, 0.3250, 0.0980];
sMR = stem(mscR);
sMR.Color = [0, 0.4470, 0.7410];
yline(vc_msc,'--r','LineWidth', 2.5);
hold off;
grid on;
title('MSC Ruidos')
legend('Ru�do 1','Ru�do 2','Limiar')

% subplot(212)
subplot(222)
sMS = stem(mscS);
yline(vc_msc,'--r','LineWidth', 2.5);
grid on;
title('MSC Sen�ide')

% figure(2)
% subplot(211)
subplot(223)
sMSR1 = stem(mscSR1);
yline(vc_msc,'--r','LineWidth', 2.5);
grid on;
title('MSC Sen�ide + Ru�do 1')

% subplot(212)
subplot(224)
sMSR2 = stem(mscSR2);
yline(vc_msc,'--r','LineWidth', 2.5);
grid on;
title('MSC Sen�ide + Ru�do 2')