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

% N�mero de pontos amostrados:
N = 1e3; % 

% Tamanho das bandas laterais = L/2:
L = 20;

% Tempo:
t = linspace(0,N-1,N)';

%% Gera��o de sinais
% Distribui��o normal (gaussiana):
dist = randn(N,1);  

% "Ru�do" aleat�rio (m�dia 1, amplitude 5):
r = (Ar*randn(N,1) + Mr);

% Sen�ide:
s = 5*sin(2*pi*f0/fs*t+th0);

% Sen�ide + Ru�do 1:
j = s+dist;

% Sen�ide + Ru�do 2:
k = s+r;

%% Gr�fico das LFTs:
bin_f0 = round(N*f0/fs+1);
freqs_LFT = L:1:round(N/2 - L+1);
lftD = zeros(round(N/2 - L+1),1);  
lftR = zeros(round(N/2 - L+1),1);   
lftS = zeros(round(N/2 - L+1),1); 
lftSR1 = zeros(round(N/2 - L+1),1); 
lftSR2 = zeros(round(N/2 - L+1),1); 

alfa = 0.05;
[~,vc_lft] = VC_LFT(L,alfa, 1);

for ff = freqs_LFT
    lftD(ff) = LFT(dist, L,ff);  
    lftR(ff) = LFT(r, L,ff);  
    lftS(ff) = LFT(s, L,ff); 
    lftSR1(ff) = LFT(j, L,ff);  
    lftSR2(ff) = LFT(k, L,ff); 
end

figure(1)
% subplot(211)
subplot(221)
hold on;
s_lftD = stem(lftD);
s_lftD.Color = [0.8500, 0.3250, 0.0980];
s_lftR = stem(lftR);
s_lftR.Color = [0, 0.4470, 0.7410];
yline(vc_lft,'--r','LineWidth', 2.5);
hold off;
grid on;
xlim([0,length(lftD)])
title('LFT Ruidos')
legend('Ru�do 1','Ru�do 2','Limiar')

% subplot(212)
subplot(222)
s_lftS = stem(lftS);
yline(vc_lft,'--r','LineWidth', 2.5);
grid on;
title('LFT Sen�ide')
xlim([0,length(lftD)])

% figure(2)
% subplot(211)
subplot(223)
s_lftSR1 = stem(lftSR1);
yline(vc_lft,'--r','LineWidth', 2.5);
grid on;
title('LFT Sen�ide + Ru�do 1')
xlim([0,length(lftD)])

% subplot(212)
subplot(224)
s_lftSR2 = stem(lftSR2);
yline(vc_lft,'--r','LineWidth', 2.5);
grid on;
title('LFT Sen�ide + Ru�do 2')
xlim([0,length(lftD)])