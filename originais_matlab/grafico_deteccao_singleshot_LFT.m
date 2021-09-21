clear all, close all, clc

%plotar valores críticos
alfa = 0.05;
L = 20; 

f = 100; 
fs = 1000; 
N = 1e3; 
energia_ruido = 1; 
SNRi = -15;
A = 10.^(SNRi/20)*(energia_ruido^2);

%
[~,VC_teorico] = VC_LFT(L,alfa,2);
s(:,1) = sin(2*pi*f/fs*[0:(N-1)]);
x = A *s + energia_ruido*randn(N,1); 

ord_LFT = zeros(N/2+1,1); 

for ii = L:1:(N/2-L)
    ord_LFT(ii) = LFT(x,L,ii);
end

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
stem(f,ord_LFT,'-k','LineWidth',lw)
grid on
hold on 
plot([0 N],[ VC_teorico,VC_teorico],'-r','LineWidth',lw)
xlim([f(L+1) f(N/2-L)])
xlabel('Frequência','fontsize',12)
ylabel('LFT','fontsize',12)
legend({'LFT';'Valor Crítico'},'edgecolor','none','fontsize',12)

% FP = mean(ord_LFT > VC_teorico)


% 

