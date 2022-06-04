%% cleanup
close all;clearvars;clc;

%% loading data

% ! Ensure data ins on Matlab's path
load Subject_ABR_Data.mat

sub1 = DB00{1}; % subject 1 in exam without stimulus
[N,sampleSize] = size(sub1);

% SUB1 = fft(sub1);
% T2 in time-domain out-performed freq.-domain T2

% V = NxL;
% matrix V = NxQ => vij = ith epoch, jth feature
% where N is the number of epochs 
% and Q are features extracted from each epoch
Q = 25;
V = nan(N,Q); 
for i = 1:1:N
    for j = 6:6:sampleSize
        V(i,j/6) = mean(sub1(i,j-5:j));
        % xi => mean feature value down the ith column
        % "In particular, the 15 ms windows following stimuli onset
        % are split into 25 segments of equal duration (0.6 ms per 
        % segment), and the mean is taken across each segment, 
        % giving an Nx25-dimensional feature matrix V"
    end
end

if isempty(find(isnan(V)~=0, 1))
    disp('Matrix properly constructed, no NaNs left.')
end

%% HOTTERLING'S T-Squared statistic
%     Inputs:
%          X - multivariate data matrix. 
%      alpha - significance level (default = 0.05).
%
%     Output:
%          N - sample-size.
%          Q - variables.
%          T2 - Hotelling's T-Squared statistic.
%          Chi-sqr. or F - the approximation statistic test.
%          df's - degrees' of freedom of the approximation statistic test.
%          P - probability that null Ho: is true.
%
% CODE ADAPTED FROM:
% Hotelling's T-Squared test for one multivariate sample. 
%  Created by A. Trujillo-Ortiz and R. Hernandez-Walls
%             Facultad de Ciencias Marinas
%             Universidad Autonoma de Baja California
%             Apdo. Postal 453
%             Ensenada, Baja California
%             Mexico.
%             atrujo@uabc.mx
%             And the special collaboration of the post-graduate students of the 2002:2
%             Multivariate Statistics Course: Karel Castro-Morales, Alejandro Espinoza-Tenorio,
%             Andrea Guia-Ramirez.
%  Copyright (C) December 2002

K       = 3;                % test the 3-staged design
FP      = zeros(1, K);      % number of false-positives
NumT    = N;                % number of tests to carry out

% from Generate_Thresholds_ExampleCode.m:
Thresholds = [8.1887,  11.0080,   13.4886]; 
Ps = zeros(1,K);

for ti=1:NumT
    X = V(1:K,:); % FIXME
    mu=zeros([1,Q]); % EXPECTED MEAN
    m=mean(X); %Mean vector from data matrix X.
    S=cov(X);  %Covariance matrix from data matrix X.
    T2=N*(m-mu)*inv(S)*(m-mu)'; %Hotelling's T-Squared statistic.

    % "The T2 statistic can then be transformed into an F 
    % statistic using:

    F = T2*(N-Q)/(Q*(N-1));

    % which follows an F-distribution with Q and N – Q degrees of
    % freedom (DOF) under H0. It is worth noting that the number 
    % of epochs N should be larger than the number of features Q, 
    % else S1 cannot be calculated.

    v1=Q;  %Numerator degrees of freedom.
    v2=N-Q;  %Denominator degrees of freedom.
    Ps(ti)=1-fcdf(F,v1,v2);  %Probability that null Ho: is true.
    
end

for ti=1:NumT
%     Ps          = rand(1,K);            % random p values: all assumptions are met, so the FPR should be exact. 
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

% disp('Detection rate:')
% disp(sum(numDetec)/numTests)



