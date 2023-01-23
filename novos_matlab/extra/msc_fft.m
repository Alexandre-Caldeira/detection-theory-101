function [ORD] = msc_fft(Y,M)

if (size(Y,2) ~= M)
    error('Tamanho da janela diferente'); 
end
% y = y(1:tamanho_janela*M,1); 
% %dividir em janela; 
% y =reshape(y,tamanho_janela,M); 
% 
% %aplicar a fft; 
% Y =fft(reshape(y,tamanho_janela,M)); %

%MSC
ORD =  abs(sum(Y,2)).^2./(M*sum(abs(Y).^2,2));