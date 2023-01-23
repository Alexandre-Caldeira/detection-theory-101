function [Thresholds,Stage_FPRs] = chesnayeThresholds(K,alphaVector,gammaVector,Resolution, validate,NumT)
%chesnayeThresholds returns efficacy and futility thresholds for the CGST
%   K:                  number of stages
%   alphaVector:        false positive rate for each stage (efficacy)
%   gammaVector:        true negative rate for each stage (futility)
%   validate:           enables Monte Carlo Simulation for validation
%   trialSamplesSize:   number of trials for validation step

digits(32)
%% Verify number of inputs
if (nargin<=2)
    gammaVector = [];
    Resolution      = (1/0.0001); % the resolution of the approximated distribution of the test statistic
    validate = false;
    NumT = NaN;
    Stage_FPRs = NaN;
end
if (nargin<=3)
    Resolution      = (1/0.0001); % the resolution of the approximated distribution of the test statistic
    validate = false;
    NumT = NaN;
    Stage_FPRs = NaN;
end
if (nargin<=4)
    validate = false;
    NumT = NaN;
    Stage_FPRs = NaN;
end

%% Calculate thresholds for the CGST

Xvalues         = 0:1/Resolution:35;                                                        % the x-axis along which the distribution of the test statistic will be defined
Null         	= chi2pdf(Xvalues,2);                                                       % null distribution for test statistic at stage k
Chi2_Norm       = Null/sum(Null);                                                           % normalise 
k               = 1;                                                                        % stage 1

if (length(gammaVector)<k)
    Thresholds(k)	= -2*log( alphaVector(k) );                                                 % for the first stage, we don't need to do any convolutions
    TruncInd_R      = round(Thresholds(k)*Resolution);                                          % The location where the stage 1 distribution will be truncated (later at stage two)

    for k = 2:K
        NullTrunc                   = Null;                                                     % reset null hip.
        NullTrunc(TruncInd_R:end)   = zeros(1, length(NullTrunc(TruncInd_R:end)));              % truncate distribution: this is the stage k-1 distribution upon entering stage k
        Null2                       = conv(Chi2_Norm, NullTrunc);                               % convolve the truncated distribution with the stage two distribution (given by Chi2_Norm) to give the distibution of the summary statistic for stage k
        Null2                       = Null2 / (sum(Null2) / (1 - sum(alphaVector(1:(k-1)))));   % normalise. Note: area is reduced due to the truncation
        TruncInd_R                  = findIndex(Null2, sum(Null2) - alphaVector(k));            % find the critical threshold. Function "findIndex" is a search function (probably not the most efficient, but it works)
        Thresholds(k)               = TruncInd_R/Resolution;                                    % translate the index to a threshold
        Null                        = Null2; 
    end
else
    aThresholds(k)	= -2*log( alphaVector(k) );                                                 % for the first stage, we don't need to do any convolutions
    gThresholds(k)	= -2*log(1- gammaVector(k) );
    TruncInd_Ra      = round(aThresholds(k)*Resolution);                                          % The location where the stage 1 distribution will be truncated (later at stage two)
    TruncInd_Rg      = round(gThresholds(k)*Resolution);  
    Null            = Null/sum(Null);             	 

    for k = 2:K
        NullTrunc                   = Null;                                                     % reset null hip.
        NullTrunc(TruncInd_Ra:end)  = zeros(1, length(NullTrunc(TruncInd_Ra:end)));              % truncate distribution: this is the stage k-1 distribution upon entering stage k
        NullTrunc(1:TruncInd_Rg)    = zeros(1, length(NullTrunc(1:TruncInd_Rg)));
        Null2                       = conv(Chi2_Norm, NullTrunc);                               % convolve the truncated distribution with the stage two distribution (given by Chi2_Norm) to give the distibution of the summary statistic for stage k
        Null2                       = Null2 / (sum(Null2) / (1 - sum(gammaVector(1:(k-1))) - sum(alphaVector(1:(k-1)))));
        TruncInd_Ra                 = findIndex(Null2, sum(Null2) - alphaVector(k));            % find the critical threshold. Function "findIndex" is a search function (probably not the most efficient, but it works)
        aThresholds(k)              = TruncInd_Ra/Resolution;                                    % translate the index to a threshold
        TruncInd_Rg                 = findIndex(Null2, gammaVector(k), 1);
        gThresholds(k)              = TruncInd_Rg/Resolution;
        Null                        = Null2; 
    end
    
    Thresholds = [aThresholds;gThresholds];
end

%% Validation step
if (validate)
    FP      = zeros(1, K);                          % number of false-positives]
    k = 1;
    for ti=1:NumT

        Ps          = rand(1,K);                    % random p values: all assumptions are met, so the FPR should be exact. 
        Plog        = -2*log(Ps);                   % fisher transform

        % check for rejections
        if sum(Plog(1:k)) > Thresholds(k)           % H0 rejected at stage 1 --> stop
            FP(k) = FP(k) + 1;
        else
            k = k+1;
        end
        
        if k==K+1
            break
        end
    end

    % total FPR
%     disp('Resulting FPR: ')
%     disp(sum(FP) / NumT)
    % stage-wise FPRs
    Stage_FPRs = FP/NumT;
%     disp('Stage_FPRs: ')
%     disp(Stage_FPRs)      
end

end

