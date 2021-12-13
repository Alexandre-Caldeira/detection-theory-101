%Exemplo 2 - Frequencia de amostragem 

clear all
close all
clc

%% Código para gerar um final senoidal 
t_final = 1; %tempo final 
tc = 0:0.0001:t_final; 
f1 = 12; 
fs = 10;
Ts = 1/fs; %fs=1/Ts   TS->0 Tempo discreto ->tempo contínuo
N = t_final/Ts; 
td = 0:Ts:(N-1)*Ts;

%gerar senoide
figure

A = 1;
teta0 = pi/3;
xc = A*cos(2*pi*tc*f1+teta0);
xd = A*cos(2*pi*f1*Ts*[0:(N-1)]+teta0);
%%
% for ii = 1:N
%     xd(ii) = A*cos(2*pi*f1*Ts*(ii)+teta0);
% end


subplot(3,1,1)
plot(tc,xc,'k','LineWidth',2)
xlabel('Tempo(s)','fontsize',12)
ylabel('Amplitude','fontsize',12)
 
subplot(3,1,2)
stem(0:(N-1),xd,'k','LineWidth',2)
xlabel('Tempo(s)','fontsize',12)
ylabel('Amplitude','fontsize',12)

subplot(3,1,3)
plot(tc,xc,'b','LineWidth',1)
xlabel('Tempo(s)','fontsize',12)
ylabel('Amplitude','fontsize',12)
hold on
stem(td,xd,'k','LineWidth',.5)
plot(td,xd,'-r','LineWidth',2)
xlabel('Tempo(s)','fontsize',12)
ylabel('Amplitude','fontsize',12)
hold off
