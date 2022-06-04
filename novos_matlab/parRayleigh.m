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
nRuns= 5000;
N = 100;
tj = N; M = N;
alpha = 0.01;

Rm0 = zeros(nRuns,N);

m91 = 10:1:100;
Rm91 = zeros(nRuns,N);

m16 = 6:6:100;
% m = Mmin:Mstep:Mmax;
Rm16 = zeros(nRuns,N);

% approx 7.1548 => Rm0 e Rm16 @ nRuns = 1 000
% approx 36.6971 => Rm0 e Rm91 @ nRuns = 1 000
% approx 41.7027 => Rm0, Rm91 e Rm16 @ nRuns = 1 000

% approx  => Rm0 e Rm16 @ nRuns = 5 000
% approx  => Rm0 e Rm91 @ nRuns = 5 000
% approx 244.8044 => Rm0, Rm91 e Rm16 @ nRuns = 5 000

% approx 72.3477 => Rm0 e Rm16 @ nRuns = 10 000
% approx  => Rm0 e Rm91 @ nRuns = 10 000
% approx 497.3268 => Rm0, Rm91 e Rm16 @ nRuns = 10 000

% estimativa: O(nRuns*16*91*N^2*)
parfor ii = 1:nRuns
    y = randn(1,N^2);
    
    Rm0(ii,:) = rayleigh(y,N,N);
    
    aux16= zeros(16,N);
    for m = 1:numel(m16)
        aux16(m,1:N) = rayleigh(y,N,m16(m));       
    end
    Rm16(ii,:) = max(aux16,[],1);
    
    aux91= zeros(91,N);
    for m = 1:numel(m91)
        aux91(m,1:N) = rayleigh(y,N,m91(m));        
    end
    Rm91(ii,:) = max(aux91,[],1);
    
end

disp(toc(t))
%%
save('workspaceRayleighDDMMAA.mat')

%% Apresenta resultados
figure

[f,xi] = ksdensity(Rm0(:,end));
limiar0 = quantile(f, 1-alpha);
p1 = plot(xi,f,'LineWidth', 2);
hold on
% h1  = histogram(Rm0(:,end));
plot([limiar0 limiar0],[0,3],'r--')

[f2,xi2] = ksdensity(max(Rm0(:,:),[],2));
limiar02 = quantile(f2, 1-alpha);
p2 = plot(xi2,f2,'LineWidth', 2);
% histogram(Rm0(:,end));
plot([limiar02 limiar02],[0,3],'--', 'Color', p2.Color)

% [f3,xi3] = ksdensity(max(Rm16,[],2));
[f3,xi3] = ksdensity(Rm16(:,end));
limiar16 = quantile(f3, 1-alpha);
p3 = plot(xi3,f3,'LineWidth', 2);
% histogram(Rm0(:,end));
plot([limiar16 limiar16],[0,3],'--', 'Color', p3.Color)

% [f4,xi4] = ksdensity(max(Rm91,[],2));
[f4,xi4] = ksdensity(Rm91(:,end));
limiar91 = quantile(f4, 1-alpha);
p4 = plot(xi4,f4,'LineWidth', 2);
% histogram(Rm0(:,end));
plot([limiar91 limiar91],[0,3],'--', 'Color', p4.Color)


hold off
grid on
xlim([0,2])

legenda = {'H_{01}',['Limiar_{H_{01}} =' num2str(limiar0)],...
           'H_{02}',['Limiar_{H_{02}} =' num2str(limiar02)],...
           'H_{0216}',['Limiar_{H_{0216}} =' num2str(limiar16)],...
           'H_{0291}',['Limiar_{H_{0291}} =' num2str(limiar91)]};

legend(legenda,'Location','NorthWest')
title('Rayleigh MCSim para N=5 000, alpha = 5%')
% /usr/local/Polyspace/R2021a/bin/matlab