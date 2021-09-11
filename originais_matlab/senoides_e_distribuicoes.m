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

% Tamanho da janela:
tj = fs; % 1seg/janela

% Número de janelas: 
M = 5;

% Número de pontos amostrados:
N = M*tj; % 

% Tempo:
t = linspace(0,N-1,N);

%% Geração de sinais
% Distribuição normal (gaussiana):
dist = randn(N,1)';  

% "Ruído" aleatório (média 1, amplitude 5):
r = (Ar*randn(N,1) + Mr)';

% Senóide:
s = 5*sin(2*pi*f0/fs*t+th0);

% Senóide + Ruído 1:
j = s+dist;

% Senóide + Ruído 2:
k = s+r;

%% Gráficos:
figure(1)
hold on

% Distribuição normal:
hd = histfit(dist);
hd(1).FaceAlpha = 0.6;
hd(1).EdgeAlpha = 0.6;
hd(1).FaceColor = [0.8500, 0.3250, 0.0980];
hd(1).EdgeColor = [0.8500, 0.3250, 0.0980];
hd(2).Color = [0.8500, 0.3250, 0.0980];

% Ruído aletório:
hr = histfit(r);
hr(1).FaceAlpha = 0.6;
hr(1).EdgeAlpha = 0.6;
hr(1).FaceColor = [0, 0.4470, 0.7410];
hr(1).EdgeColor = [0, 0.4470, 0.7410];
hr(2).Color = [0, 0.4470, 0.7410];

limiar = quantile(dist,1-0.05);
px = xline(limiar);
px.Color = 'r';
px.LineStyle = '--';
px.LineWidth = 1.8;

hold off
grid on
Ldists = legend('Ruído 1:','$\mathcal{N}(\mu = 0,\sigma^2 = 1^2)$',...
           'Ruído 2:','$\mathcal{N}(\mu = 10,\sigma^2 = 5^2)$',...
           'Limiar (teste F da média)');
       
warning('off')
Ldists.Interpreter = 'latex';
Ldists.FontSize = 14;

xlabel('Amplitude [V]','FontSize',14)
ylabel('Quantidade de ocorrências','FontSize',14)
title('Visualização Distribuições aleatórias','FontSize',14)

% Senóide:
figure(2)
subplot(221)
ps = plot(t(1:101),s(1:101));
ylabel('Amplitude [V]','FontSize',14)
xlabel('Tempo [s]','FontSize',14)
title('Senóide amostrada sem ruído','FontSize',14)
Ls = legend('$f_0 = 100$ [Hz], $f_s = 1000$ [Hz]','FontSize',14);
grid on
Ls.Interpreter = 'latex';

% FFT:
subplot(222)
S = abs(fft(s));
% bin_f0 = round(N*f0/fs+1);
freqS = round((fs/N).*(linspace(0,length(S),N)),2);
pfs = stem(freqS(1:round(size(S,2)/2)),S(1:round(size(S,2)/2)));
xlim([0 max(freqS(1:round(size(S,2)/2)))])
ylabel('Magnitude','FontSize',14)
xlabel('Frequência [Hz]','FontSize',14)
title('FFT do sinal ao lado','FontSize',14)
grid on

% Senóide + Ruído 1:
subplot(223)
ps2 = plot(t(1:101),j(1:101),'LineWidth',2);
ps2.Color = [0.8500, 0.3250, 0.0980];
hold on
ps12 = plot(t(1:101),s(1:101));
ps12.Color = [0, 0.4470, 0.7410];

ylabel('Amplitude [V]','FontSize',14)
xlabel('Tempo [s]','FontSize',14)
title('Senóide amostrada com ruído 1','FontSize',14)
Ls = legend('$sen(\cdot) +\mathcal{N}(\mu = 0,\sigma^2 = 1^2)$','$f_0 = 100$ [Hz], $f_s = 1000$ [Hz]','FontSize',14);
grid on
Ls.Interpreter = 'latex';

subplot(224)
J = abs(fft(j));
% bin_f0 = round(N*f0/fs+1);
freqS = round((fs/N).*(linspace(0,length(J),N)),2);
pfs2 = stem(freqS(1:round(size(J,2)/2)),J(1:round(size(J,2)/2)));
xlim([0 max(freqS(1:round(size(J,2)/2)))])
ylabel('Magnitude','FontSize',14)
xlabel('Frequência [Hz]','FontSize',14)
title('FFT do sinal ao lado','FontSize',14)
grid on

figure(3)
subplot(121)
ps2 = plot(t(1:101),k(1:101),'LineWidth',2);
ps2.Color = [0.8500, 0.3250, 0.0980];
hold on
ps12 = plot(t(1:101),s(1:101));
ps12.Color = [0, 0.4470, 0.7410];

ylabel('Amplitude [V]','FontSize',14)
xlabel('Tempo [s]','FontSize',14)
title('Senóide amostrada com ruído 2','FontSize',14)
Ls = legend('$sen(\cdot) + \mathcal{N}(\mu = 10,\sigma^2 = 5^2)$','$f_0 = 100$ [Hz], $f_s = 1000$ [Hz]','FontSize',14);
grid on
Ls.Interpreter = 'latex';

subplot(122)
K = abs(fft(k));
% bin_f0 = round(N*f0/fs+1);
freqS = round((fs/N).*(linspace(0,length(K),N)),2);
pfs3 = stem(freqS(1:round(size(K,2)/2)),K(1:round(size(K,2)/2)));
xlim([0 max(freqS(1:round(size(K,2)/2)))])
ylabel('Magnitude','FontSize',14)
xlabel('Frequência [Hz]','FontSize',14)
title('FFT do sinal ao lado','FontSize',14)
grid on

