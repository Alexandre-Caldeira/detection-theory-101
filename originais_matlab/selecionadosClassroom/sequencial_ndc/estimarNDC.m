function [NDC,FP] = estimarNDC(NDCinicial, alfa, FP_desejado, ord, Mmin, Mstep, Mmax)



%parametros de entrada 
% NDCinicial  = 10; 
% Mmin = 10; 
% Mstep = 1; 
% Mmax = 30; 
%ord ->vetor com as ORDs - cada linha uma simulação 


MM = Mmin:Mstep:Mmax;
NNTmax = size(MM,2); %número de teste
nRuns =  size(ord,1);
fprintf('Mmin =%d - Mstep=%d - Mmax=%d - Nruns=%d  -NDCinicial=%d\n',Mmin,Mstep, Mmax,nRuns,NDCinicial); 

%cálcular o falso positivo
dr = zeros(nRuns,1);
FP =0; 
for ii = 1:nRuns    
        dr(ii,1) = ETS(ord(ii,:),MM,alfa,NDCinicial); 
end
FP = mean(dr);    

if (FP<FP_desejado) 
    if NDCinicial>1
        warning('NDC inicial grande --- Isso pode aumentar o tempo do código');
        NDCinicial = 1; 
    end 
end

FP = 0; 
for NDC = NDCinicial:NNTmax %número máximo de teste           
            
    for ii = 1:nRuns    
        dr(ii,1) = ETS(ord(ii,:),MM,alfa,NDC); 
    end
    FP(NDC) = mean(dr);  
    
    if (FP(NDC)<FP_desejado)
        break; 
    end
    
end    
    

%encontrar o valor mínimo
[~,ind] = find((FP<FP_desejado).*(FP~=0)); %encontrar o primeiro NT que respeita a restrição

NDC = 0;
if isempty(ind)
    NDC = nan; 
    warning('####### Não encontro Nalpha ##########');
    elseif ind(1)>1
        NDC = interp1(FP([ind(1)-1,ind(1)]),[ind(1)-1,ind(1)],FP_desejado);
    else
        NDC = 1; %para apenas uma detecção 
end
    
if(NNTmax ==2)
    NDC=2;
end

if(NNTmax ==1)
    NDC=1;
end

