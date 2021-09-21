clear all, close all, clc

%plotar valores críticos
alfa = 0.05;

alpha = 0.05;
f = 100; 
fs = 1000; 
tj = 100;
M = 300;
N=M*tj+1;
bin = round(tj*f/fs +1);

energia_ruido = 1; 
SNRi = -30;
A = 10.^(SNRi/20)*(energia_ruido^2);

%
[~,VC_teorico] = VC_CSM(M,alfa,1);
s(:,1) = sin(2*pi*f/fs*[0:(N-1)]);
x = A *s + energia_ruido*randn(N,1); 

% ord_LFT = zeros(N/2+1,1); 

ord_LFT = msc(x,tj,M);


%%
f = [0:N/2]*fs/N;
figure 
lw = 1.5;
subplot(1,2,1)
plot([1:100],x(1:100),'-r','LineWidth',lw)
hold on
grid on
plot([1:100],s(1:100),'-b','LineWidth',lw)
ylabel('Tensão [V]')
xlabel('Tempo discreto [amostras]')

subplot(1,2,2)
stem(ord_LFT,'-k','LineWidth',lw)
grid on
hold on 
plot([0 length(ord_LFT)],[ VC_teorico,VC_teorico],'-r','LineWidth',lw)
% xlim([f(L+1) f(N/2-L)])
xlabel('Frequência Observada','fontsize',12)
ylabel('MSC','fontsize',12)
legend({'MSC';'Valor Crítico'},'edgecolor','none','fontsize',12)

% FP = mean(ord_LFT > VC_teorico)


% 

