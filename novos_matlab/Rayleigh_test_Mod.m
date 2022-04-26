%% Valor crítico para testes sequenciais (Stürzebecher, E., Cebulla, M., & Elberling, C. (2005))
clear
clc
M = 100;    % Número de janelas
Mmin = 10;  % Número mínimo de janelas para realizar o primeiro teste
Mstep = 1;  % Quantidade de janelas a ser acrescentada para realizar outro teste
MM = Mmin:Mstep:M;  % Quantidades de janelas que ocorrerão os testes
QT = length(MM);    % Quantidade de testes
a = 0.01   % Nível de significância
nRuns = 100000;
Rm = zeros(nRuns,QT);
t = tic();
parfor k = 1:nRuns
    % Modified Rayleigh Test
    YY = randn(1,M) + 1j*randn(1,M);       % Simulando as componenetes espectrais de cada janela na frequência de interesse.
    for j = 1:QT
        Y = sort(YY(1:MM(j)));%,'ComparisonMethod','abs'); % Colocar em ordem crescente em relação ao módulo
        Y = (Y./abs(Y)).*(1:MM(j)); % Trocando o módulo pelo rank
        Cm = mean(real(Y));
        Sm = mean(imag(Y));
        Rm(k,j) = (sqrt(Cm^2 + Sm^2))/sqrt(MM(j));
    end
end
disp(toc(t))

%% Apresenta resultados

alpha = 0.01;
figure
[f,xi] = ksdensity(Rm(:,end));
%[f,xi] = ksdensity(max(Rm0(:,:),[],2));
limiar0 = quantile(Rm(:,end), 1-alpha);
f = f/trapz(f);
plot(xi,f,'LineWidth', 2);
hold on
plot([limiar0 limiar0],[0,0.05],'b--')

%[f2,xi2] = ksdensity(max(Rm0(:,:),[],2));
%limiar02 = quantile(f2, 1-alpha);
%plot(xi2,f2,'LineWidth', 2);
%plot([limiar02 limiar02],[0,3],'g--')

%[f3,xi3] = ksdensity(max(Rm162,[],2));
[f3,xi3] = ksdensity(max(Rm(:,1:6:end),[],2));
limiar16 = quantile(max(Rm(:,1:6:end),[],2), 1-alpha);
f3 = f3/trapz(f3);
plot(xi3,f3,'LineWidth', 2);
plot([limiar16 limiar16],[0,0.05],'y--')

%[f4,xi4] = ksdensity(max(Rm91,[],2));
[f4,xi4] = ksdensity(max(Rm,[],2));
limiar91 = quantile(max(Rm,[],2), 1-alpha);
f4 = f4/trapz(f4);
plot(xi4,f4,'LineWidth', 2);
plot([limiar91 limiar91],[0,0.05],'r--')


%limiarBon = quantile(Rm0(:,end),1-alpha/91);
%plot([limiarBon limiarBon],[0,0.05],'b--')

hold off
grid on
xlim([0,2])

%legenda = {'H_{01}',['Limiar_{01} =' num2str(limiar0)],...
%           'H_{0216}',['Limiar_{0216} =' num2str(limiar16)],...
%           'H_{0291}',['Limiar_{0291} =' num2str(limiar91)]};
%           'H_{01Bon}',['Limiar_{01Bon} =' num2str(limiarBon)]};

%legend(legenda,'Location','NorthWest')