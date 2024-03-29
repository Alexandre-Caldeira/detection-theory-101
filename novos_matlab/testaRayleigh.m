%% 03/09/22
% Replica Cebulla,2005:
% -------------------------------------------------------------------------
% Automated auditory response detection: 
% Statistical problems with repeated testing Evaluaci�n repetida en la 
% detecci�n de respuestas auditivas. 
% 
% International Journal of Audiology, 44(2), 110�117 
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
% Response (AMFR):Detecti�n objetiva de la respuesta consecuente de 
% amplitud modulada (AMFR). 
% 
% International Journal of Audiology, 40(5), 245�252 
% doi: 10.3109/00206090109073118
% https://sci-hub.se/10.3109/00206090109073118 
% -------------------------------------------------------------------------
% Limpar workspace
clearvars; clc; close all;


% Parametros de simulacao
t = tic();
nRuns=1/2 *10000;nRuns=1000;
N = 100;
tj = N; M = N;
alpha = 0.01;

Rm0 = zeros(nRuns,N);

aux91= zeros(91,N);
Rm91 = zeros(nRuns,N);

aux16= zeros(16,N);
Rm16 = zeros(nRuns,N);

% approx 115.6759s => Rm0 e Rm16 @ nRuns = 5000
% approx 718.6909 => Rm0 e Rm16 @ nRuns = 5000
for ii = 1:nRuns
    if rem(ii,100)==0
        disp(ii)
    end
    y = randn(1,N^2);
    Rm0(ii,:) = rayleigh(y,N,N);
    
    for M = 6:6:100
        aux16(M/6,:) = rayleigh(y,N,M);
    end
    Rm16(ii,:) = max(aux16,[],1);
    
    %(Mmax-Min+1)/Mstep
    for M = 10:1:100
        aux91(M-9,:) = rayleigh(y,N,M);
    end
    Rm91(ii,:) = max(aux91,[],1);
end

disp(toc(t))
%% Apresenta resultados
figure

[f,xi] = ksdensity(Rm0(:,end));
limiar0 = quantile(f, 1-alpha);
plot(xi,f,'LineWidth', 2);
hold on
plot([limiar0 limiar0],[0,3],'r--')

[f2,xi2] = ksdensity(max(Rm0(:,:),[],2));
limiar02 = quantile(f2, 1-alpha);
plot(xi2,f2,'LineWidth', 2);
plot([limiar02 limiar02],[0,3],'r--')

% [f3,xi3] = ksdensity(max(Rm16,[],2));
[f3,xi3] = ksdensity(Rm16(:,end));
limiar16 = quantile(f3, 1-alpha);
plot(xi3,f3,'LineWidth', 2);
plot([limiar16 limiar16],[0,3],'r--')

% [f4,xi4] = ksdensity(max(Rm91,[],2));
[f4,xi4] = ksdensity(Rm91(:,end));
limiar91 = quantile(f4, 1-alpha);
plot(xi4,f4,'LineWidth', 2);
plot([limiar91 limiar91],[0,3],'r--')


hold off
grid on
xlim([0,2])

legenda = {'H_{01}',['Limiar_{01} =' num2str(limiar0)],...
           'H_{02}',['Limiar_{02} =' num2str(limiar02)],...
           'H_{0216}',['Limiar_{0216} =' num2str(limiar16)],...
           'H_{0291}',['Limiar_{0291} =' num2str(limiar91)]};

legend(legenda,'Location','NorthWest')