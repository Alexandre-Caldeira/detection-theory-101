function [ORD] = CSM(y,tamanho_janela,M)


y = y(1:tamanho_janela*M,1); 
%dividir em janela; 
y =reshape(y,tamanho_janela,M); 

%aplicar a fft; 
Y =fft(y); %

teta = angle(Y); %pegar o angulo; 
ORD =  (1/M^2)*sum(cos(teta),2).^2+(1/M^2)*sum(sin(teta),2).^2;


