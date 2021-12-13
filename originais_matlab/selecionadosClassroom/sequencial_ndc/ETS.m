%fun��o da estrat�gia para testes sequenciais 
function [dr,time] = ETS(ord,MM,alfa,NDC)

%parametros: 
% MM = Mmin:Mstep:Mmax; %parametros aonde vai aplicar o teste


%% 

NDC = ceil(NDC);
if size(MM,2)>1
    MM = MM';
end

if size(ord,2)>1
    ord = ord';
end


%VALOR CR�TICO - MUDA PARA CADA DETECTOR--------
valor_critico = VC_MSC(MM,alfa);  
det = ord(MM)>valor_critico; 
%-----------------------------------------------------

%avaliar se atende o NDC 
cont = 0; 
dr =0; 
for ii = 1:size(MM,1) 
   
    cont = det(ii)+cont*det(ii); %conta o n�mero de detec��es consecutivas 
    
    if cont == NDC
        dr = 1;
        time = MM(ii);
        break
    end
end

if dr == 0 
    time = MM(ii);
end

 