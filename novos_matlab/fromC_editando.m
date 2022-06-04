%% 
clear all
close all
clc

% parameters
K               = 3;                            % analyse data 3 times
TotalAlpha      = 0.05;                         % pre-specified alpha-level for the full sequential test
Alpha_k         = ones(1,3)*(TotalAlpha/3);     % alpha to spend at each stage
Resolution      = (1/0.0001);                   % the resolution of the approximated distribution of the test statistic
Xvalues         = 0:1/Resolution:35;            % the x-axis along which the distribution of the test statistic will be defined
Null         	= chi2pdf(Xvalues,2);           % null distribution for test statistic at stage k
Chi2_Norm       = Null/sum(Null);             	% normalise 

sum(Null/Resolution)

% stage 1 threshold
k               = 1;                                % stage 1
Thresholds(k)	= -2*log( Alpha_k(k) );             % for the first stage, we don't need to do any convolutions
TruncInd_R      = round(Thresholds(k)*Resolution);  % The location where the stage 1 distribution will be truncated (later at stage two)

% stage 2 threshold
k                           = 2;                                                        % stage 2
NullTrunc                   = Null;                                                     % 
NullTrunc(TruncInd_R:end)   = zeros(1, length(NullTrunc(TruncInd_R:end)));              % truncate distribution: this is the stage 1 distribution upon entering stage 2
Null2                       = conv(Chi2_Norm, NullTrunc);                               % convolve the truncated distribution with the stage two distribution (given by Chi2_Norm) to give the distibution of the summary statistic for stage 2
Null2                       = Null2 / ( sum(Null2) / (1 - sum(Alpha_k(1:(k-1)) ) ) );       % normalise. Note: area is reduced due to the truncation
TruncInd_R                  = findIndex(Null2, sum(Null2) - Alpha_k(k));                % find the critical threshold. Function "findIndex" is a search function (probably not the most efficient, but it works)
Thresholds(k)               = TruncInd_R/Resolution;                                    % translate the index to a threshold
Null                        = Null2;                                                    % 

sum(Null/Resolution)

% stage 3 threshold
k                           = 3;                                                % stage 3
NullTrunc                   = Null;                                             %
NullTrunc(TruncInd_R:end)   = zeros(1, length(NullTrunc(TruncInd_R:end)));      % truncation
Null2                       = conv(Chi2_Norm, NullTrunc);                       % convolve the stage two truncated distribution with the stage three distribution (given again by Chi2_Norm) to give the distibution of the summary statistic for stage 3
Null2                       = Null2 / (sum(Null2) / sum(NullTrunc));            % normalise 
TruncInd_R                  = findIndex(Null2, sum(Null2) - Alpha_k(k));        % find truncation index
Thresholds(k)               = TruncInd_R/Resolution;                            % and translate the index to a threshold   


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % verify FPR  % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
K       = 3;                % test the 3-staged design
FP      = zeros(1, K);      % number of false-positives
NumT    = 1000000;          % number of tests to carry out
for ti=1:NumT
    
    Ps          = rand(1,K);            % random p values: all assumptions are met, so the FPR should be exact. 
    Plog        = -2*log(Ps);           % fisher transform
    Detected    = false;                % 
    
    % check for rejections
    k = 1;
    if sum(Plog(1:k)) > Thresholds(k)       % H0 rejected at stage 1 --> stop
        FP(k) = FP(k) + 1;
    else
        k = 2;
        if sum(Plog(1:k)) > Thresholds(k)       % H0 rejected at stage 2 --> stop
            FP(k) = FP(k) + 1;
        else
            k = 3;
            if sum(Plog(1:k)) > Thresholds(k)   % H0 rejected at stage 3 --> stop
                FP(k) = FP(k) + 1;
            end
        end
    end
end

Stage_FPRs  = FP/NumT               % stage-wise FPRs
FPR         = sum(FP) / NumT        % total FPR

