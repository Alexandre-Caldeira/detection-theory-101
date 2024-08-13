function [aThresholds,gThresholds] = get_beta_CGST_thresholds(K, M, alpha )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    Alpha_k         = ones(1,K)*(alpha/K);    
    Gamma_k = ((1-alpha)/K).*ones(1,K);
    Resolution      = (1/0.00001); %(1/0.0001);                 
    Xvalues         = 0:1/Resolution:1;            
    Null         	= betapdf(Xvalues, 1, M-1);
    Null            = Null/sum(Null);            
    Chi2_Norm       = Null/sum(Null);             
 
    k               = 1;                           
    aThresholds(k)	= 1 - Alpha_k(k).^(1./(M-1));  
    gThresholds(k)	= 1-(1- Gamma_k(k)).^(1./(M-1));
    TruncInd_Ra      = round(aThresholds(k)*Resolution);
    TruncInd_Rg      = round(gThresholds(k)*Resolution);           
    
    for k = 2:K
        NullTrunc                   = Null;                                                
        NullTrunc(TruncInd_Ra:end)  = zeros(1, length(NullTrunc(TruncInd_Ra:end)));    
        NullTrunc(1:TruncInd_Rg)    = zeros(1, length(NullTrunc(1:TruncInd_Rg)));
        
        Null2                       = conv(Chi2_Norm, NullTrunc);   
        Null2                       = Null2 / (sum(Null2) / (1 - sum(Gamma_k(1:(k-1))) - sum(Alpha_k(1:(k-1)))));

        TruncInd_Ra                 = findIndex(Null2, sum(Null2) - Alpha_k(k)); 
        aThresholds(k)              = TruncInd_Ra/Resolution;  
        TruncInd_Rg                 = findIndex(Null2, Gamma_k(k), 1);
        gThresholds(k)              = TruncInd_Rg/Resolution;
        Null                        = Null2; 
    end   

end