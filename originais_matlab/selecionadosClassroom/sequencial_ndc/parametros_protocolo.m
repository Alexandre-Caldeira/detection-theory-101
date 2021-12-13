function P = parametros_protocolo(Mmax)

%fun��o para obter todos os par�metros poss�veis de um protocolo de
%detec��o a partir de um Mmax de entrada 

%entrada: Mmax 

%sa�da: Matriz P onde: 
        %cada linha indica um par�metro. 
        %Coluna [ Min, Mstep, Mmax]; 

%%C�digo ------------------------------------------------

P = []; 
 
for M_step = 1:(Mmax-1)
    for Mmin = 2:(Mmax-1)%MMentrada                 
         k = (Mmax -Mmin)/M_step;
         if (ceil(k) == floor(k) & k>=0) %verificar se o valor de k � um n�mero real    
             P = [P;Mmin,M_step,Mmax];             
         end       
       end
end

%Incluir o teste �nico tamb�m 
P = [P;Mmax 1 Mmax]; 




