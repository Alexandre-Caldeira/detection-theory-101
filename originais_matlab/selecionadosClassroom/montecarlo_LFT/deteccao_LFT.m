clear all, close all, clc

%plotar valores críticos
alfa = 0.05;
L = 10; 

f = 10; 
fs = 100; 
N = 1000; 
amplitude_ruido = 10; 
A =1;

%
[~,VC_teorico] = VC_LFT(L,alfa,2);
s(:,1) = sin(2*pi*f/fs*[0:(N-1)]);
x = A *s + amplitude_ruido*randn(N,1); 

ord_LFT = zeros(N/2+1,1); 

for ii = L:1:(N/2-L)
    ord_LFT(ii) = LFT(x,L,ii);
end

f = [0:N/2]*fs/N;
figure 
lw = 1.5;
subplot(2,1,1)
plot([1:100],x(1:100),'-k','LineWidth',lw)
subplot(2,1,2)
plot(f,ord_LFT,'-k','LineWidth',lw)
hold on 
plot([0 N],[ VC_teorico,VC_teorico],'-r','LineWidth',lw)
xlim([f(L+1) f(N/2-L)])
xlabel('Frequência','fontsize',12)
ylabel('LFT','fontsize',12)
legend({'LFT';'Valor Crítico'},'edgecolor','none','fontsize',12)

% FP = mean(ord_LFT > VC_teorico)


% 

