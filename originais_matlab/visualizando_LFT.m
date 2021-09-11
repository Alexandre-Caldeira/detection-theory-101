% 09/Set/21 - Alexandre Caldeira
%% Boas práticas
clearvars;close all; clc;

%% Constantes:
% Amplitude do Ruído:
Ar = 5;

% Média do ruído:
Mr = 10;

% Frequência do Sinal: [Hz]
f0 = 100;

% Fase inicial do sinal: [rad]
th0 = 0;

% Frequência de Amostragem: [Hz]
fs = 1000;

% Número de pontos amostrados:
N = 1e3; % 

% Tamanho das bandas laterais = L/2:
L = 20;

% Tempo:
t = linspace(0,N-1,N)';

%% Geração de sinais
% Distribuição normal (gaussiana):
dist = randn(N,1);  

% "Ruído" aleatório (média 1, amplitude 5):
r = (Ar*randn(N,1) + Mr);

% Senóide:
s = 5*sin(2*pi*f0/fs*t+th0);

% Senóide + Ruído 1:
j = s+dist;

% Senóide + Ruído 2:
k = s+r;

%% Gráfico das LFTs:
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
legend('Ruído 1','Ruído 2','Limiar')

% subplot(212)
subplot(222)
s_lftS = stem(lftS);
yline(vc_lft,'--r','LineWidth', 2.5);
grid on;
title('LFT Senóide')
xlim([0,length(lftD)])

% figure(2)
% subplot(211)
subplot(223)
s_lftSR1 = stem(lftSR1);
yline(vc_lft,'--r','LineWidth', 2.5);
grid on;
title('LFT Senóide + Ruído 1')
xlim([0,length(lftD)])

% subplot(212)
subplot(224)
s_lftSR2 = stem(lftSR2);
yline(vc_lft,'--r','LineWidth', 2.5);
grid on;
title('LFT Senóide + Ruído 2')
xlim([0,length(lftD)])