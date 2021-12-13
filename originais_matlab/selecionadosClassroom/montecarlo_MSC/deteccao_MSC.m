clear all, close all, clc

%plotar valores críticos
alfa = 0.05;
M= 50;
tamanho_janela = 100; 
f = 10; 
fs = 100; 
N = M *tamanho_janela; %duração do sinal

amplitude_ruido = 20; 
A =1;

%
[~,VC_teorico] = VC_MSC(M,alfa,2);
s(:,1) = sin(2*pi*f/fs*[0:(N-1)]);
x = A *s + amplitude_ruido*randn(N,1); 

ord_MSC = MSC(x,tamanho_janela,M);


f = [0:tamanho_janela/2]*fs/tamanho_janela;
figure 
lw = 1.5;
subplot(2,1,1)
plot([1:100],x(1:100),'-k','LineWidth',lw)
subplot(2,1,2)
plot(f,ord_MSC(1:(tamanho_janela/2+1)),'-k','LineWidth',lw)
hold on 
plot([0 N],[ VC_teorico,VC_teorico],'-r','LineWidth',lw)
xlim([f(2) f(end)])
xlabel('Frequência','fontsize',12)
ylabel('MSC','fontsize',12)
legend({'MSC';'Valor Crítico'},'edgecolor','none','fontsize',12)

% FP = mean(ord_LFT > VC_teorico)


% 

