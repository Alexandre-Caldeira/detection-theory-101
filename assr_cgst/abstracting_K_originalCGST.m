%% 
clear all
close all
clc

% parameters
K               = 5;                            % analyse data 3 times
TotalAlpha      = 0.05;                         % pre-specified alpha-level for the full sequential test
Alpha_k         = ones(1,K)*(TotalAlpha/K);     % alpha to spend at each stage
% Gamma_k         = [0.1, 0.15, 0.2, 0.25, 0.29]; % ...!
Gamma_k = ((1-TotalAlpha)/K).*ones(1,K);
Resolution      = (1/0.0001);                   % the resolution of the approximated distribution of the test statistic
Xvalues         = 0:1/Resolution:35;            % the x-axis along which the distribution of the test statistic will be defined
Null         	= chi2pdf(Xvalues,2);           % null distribution for test statistic at stage k
Null            = Null/sum(Null);             	% normalise 
Chi2_Norm       = Null/sum(Null);             	% normalise 

% stage 1 threshold
k               = 1;                                % stage 1
aThresholds(k)	= -2*log( Alpha_k(k) );                                                 % for the first stage, we don't need to do any convolutions
gThresholds(k)	= -2*log(1- Gamma_k(k) );
TruncInd_Ra      = round(aThresholds(k)*Resolution);                                          % The location where the stage 1 distribution will be truncated (later at stage two)
TruncInd_Rg      = round(gThresholds(k)*Resolution);           

for k = 2:K
    disp(k)
    NullTrunc                   = Null;                                                     % reset null hip.
    NullTrunc(TruncInd_Ra:end)  = zeros(1, length(NullTrunc(TruncInd_Ra:end)));              % truncate distribution: this is the stage k-1 distribution upon entering stage k
    NullTrunc(1:TruncInd_Rg)    = zeros(1, length(NullTrunc(1:TruncInd_Rg)));
    
    Null2                       = conv(Chi2_Norm, NullTrunc);                               % convolve the truncated distribution with the stage two distribution (given by Chi2_Norm) to give the distibution of the summary statistic for stage k
    Null2                       = Null2 / (sum(Null2) / (1 - sum(Gamma_k(1:(k-1))) - sum(Alpha_k(1:(k-1)))));
    % Null2                       = Null2 / (1 - sum(Gamma_k(1:(k-1))) - sum(Alpha_k(1:(k-1))));   % normalise. Note: area is reduced due to the truncation
    
    TruncInd_Ra                  = findIndex(Null2, sum(Null2) - Alpha_k(k));            % find the critical threshold. Function "findIndex" is a search function (probably not the most efficient, but it works)
    aThresholds(k)              = TruncInd_Ra/Resolution;                                    % translate the index to a threshold
    TruncInd_Rg                 = findIndex(Null2, Gamma_k(k), 1);
    gThresholds(k)              = TruncInd_Rg/Resolution;
    Null                        = Null2; 
end
disp('')
disp('--------------------------------')
disp(aThresholds)
disp(gThresholds)
% Thresholds =[aThresholds;gThresholds]
Thresholds = aThresholds;
%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % verify FPR  % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% K       = 5;                % test the 5-staged design
FP      = zeros(1, K);      % number of false-positives
TN      = zeros(1, K);      % number of true-negatives
NumT    = 1*1e6;          % number of tests to carry out

for ti=1:NumT
    Ps          = rand(1,K);            % random p values: all assumptions are met, so the FPR should be exact. 
    Plog        = -2*log(Ps);           % fisher transform
    Detected    = false;                % 
    
    % check for rejections
    for k=1:K
        if sum(Plog(1:k)) >= aThresholds(k)
            FP(k) = FP(k)+1;
            break
        elseif sum(Plog(1:k)) <= gThresholds(k)
            TN(k) = TN(k)+1;
            break
        end
    end
end

Stage_FPRs  = FP/NumT               % stage-wise FPRs
FPR         = sum(FP) / NumT        % total FPR

Stage_TNRs  = TN/NumT               % stage-wise FPRs
TNR         = sum(TN) / NumT        % total FPR