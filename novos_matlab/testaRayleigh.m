%% 03/09/22
% Replica Cebulla,2005:
% -------------------------------------------------------------------------
% Automated auditory response detection: 
% Statistical problems with repeated testing Evaluación repetida en la 
% detección de respuestas auditivas. 
% 
% International Journal of Audiology, 44(2), 110–117 
% doi: 10.1080/14992020400029228
% https://sci-hub.se/10.1080/14992020400029228
%
% Referencias:
% -------------------------------------------------------------------------
% A Modification of the Rayleigh Test for Vector Data. 
% Biometrika, 67(1), 175 
% doi: 10.2307/2335330
% https://sci-hub.se/10.2307/2335330
% -------------------------------------------------------------------------
% Sci-Hub | Objective Detection of the Amplitude Modulation Following 
% Response (AMFR):Detectión objetiva de la respuesta consecuente de 
% amplitud modulada (AMFR). 
% 
% International Journal of Audiology, 40(5), 245–252 
% doi: 10.3109/00206090109073118
% https://sci-hub.se/10.3109/00206090109073118 
% -------------------------------------------------------------------------
% Limpar workspace
clearvars; clc; close all;


% Parametros de simulacao
t = tic();
nRuns=1/2 *10000;
N = 100;
tj = N; M = N;

Rm0 = zeros(nRuns,N);

aux91= zeros(91,N);
Rm91 = zeros(nRuns,N);

aux16= zeros(16,N);
Rm16 = zeros(nRuns,N);

% approx 115.6759s => Rm0 e Rm16 @ nRuns = 5000
for ii = 1:nRuns
    y = randn(1,N^2);
    [Y,Rm0(ii,:)] = rayleigh(y,N,N);
    
    for M = 6:6:100
        [Y,aux16(M/6,:)] = rayleigh(y,N,M);
    end
    Rm16(ii,:) = max(aux16,[],1);
    
    for M = 10:1:100
        [Y,aux91(M-9,:)] = rayleigh(y,N,M);
    end
    Rm91(ii,:) = max(aux91,[],1);
end

disp(toc(t))
%% Apresenta resultados
figure

[f,xi] = ksdensity(Rm0(:,end));
plot(xi,f,'LineWidth', 2);
hold on

[f2,xi2] = ksdensity(max(Rm0(:,:),[],2));
plot(xi2,f2,'LineWidth', 2);

[f3,xi3] = ksdensity(max(Rm16,[],2));
plot(xi3,f3,'LineWidth', 2);

[f4,xi4] = ksdensity(max(Rm91,[],2));
plot(xi4,f4,'LineWidth', 2);

hold off
grid on
xlim([0,2])
legend('H_{01}','H_{02}','H_{0216}' ,'H_{0291}')