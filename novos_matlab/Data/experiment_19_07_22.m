% Referencias:
% Sci-Hub | The Convolutional Group Sequential Test: reducing test time for evoked potentials. IEEE Transactions on Biomedical Engineering, 1–1 | 10.1109/tbme.2019.2919696
% https://sci-hub.se/10.1109/tbme.2019.2919696
% 
% Sci-Hub | A group sequential test for ABR detection. International Journal of Audiology, 1–10 | 10.1080/14992027.2019.1625486
% https://sci-hub.se/10.1080/14992027.2019.1625486

%% cleanup
close all;clearvars;clc;

%% loading data

% ! Ensure data is on Matlab's path
load Subject_ABR_Data.mat
load 5Thresholds.mat

for SS = 1:12
subN = SS;
sub1 = DB50{subN};
% sub1 = DB40{1};
% DB00{1};  subject 1 in exam without stimulus
[N,sampleSize] = size(sub1(:,1:75));
sub1 = sub1(:,1:75);
% SUB1 = fft(sub1);
% T2 in time-domain out-performed freq.-domain T2

% V = NxL;
% matrix V = NxQ => vij = ith epoch, jth feature
% where N is the number of epochs 
% and Q are features extracted from each epoch
Q = 25;
N = 3000;
V = nan(N,Q); 
for i = 1:1:N
    for j = 3:3:sampleSize
        V(i,j/3) = mean(sub1(i,j-2:j));
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


%% 

% DESIRED_ALPHA = 0.15; % 0.01;
% gamma = [0.2, 0.4,0.25];
% p = hotT2(X);



K       = 5;                % test the 3-staged design
FP      = zeros(1, K);      % number of false-positives
NumT    = 1000000;        % number of tests to carry out
% X1 = V(1:1000,:);
% X2 = V(1001:2000,:);
% X3 = V(2001:3000,:);
% Ps = [hotT2(X1),hotT2(X2),hotT2(X3)];

X1 = V(1:600,:);
X2 = V(601:1200,:);
X3 = V(1201:1800,:);
X4 = V(1801:2400,:);
X5 = V(2401:3000,:);

Ps = [hotT2(X1),hotT2(X2),hotT2(X3),hotT2(X4),hotT2(X5)];
Plog        = -2*log(Ps);           % fisher transform
aThresholds = Thresholds(1,:);
gThresholds = Thresholds(2,:);
for k=1:K
    Psum(k) = sum(Plog(1:k));
end

disp('subN = ')
disp(subN)
disp('Ps = ')
disp(Ps)
disp('Plog = ')
disp(Plog)
disp('Summary = ')
disp(Psum)
disp('alpha thresholds = ')
disp(aThresholds)
disp('gamma thresholds = ')
disp(gThresholds)

%%
for ti=1:NumT
%     Ps = 
%     Ps          = rand(1,K);            % random p values: all assumptions are met, so the FPR should be exact. 
    Detected    = false;                % 
    
    % check for rejections
    k = 1;
    if sum(Plog(1:k)) > aThresholds(k)       % H0 rejected at stage 1 --> stop
        FP(k) = FP(k) + 1;
    else
        k = 2;
        if sum(Plog(1:k)) > aThresholds(k)       % H0 rejected at stage 2 --> stop
            FP(k) = FP(k) + 1;
        else
            k = 3;
            if sum(Plog(1:k)) > aThresholds(k)   % H0 rejected at stage 3 --> stop
                FP(k) = FP(k) + 1;
            else
                k = 4;
                if sum(Plog(1:k)) > aThresholds(k)   % H0 rejected at stage 3 --> stop
                    FP(k) = FP(k) + 1;
                else
                    k = 5;
                    if sum(Plog(1:k)) > aThresholds(k)   % H0 rejected at stage 3 --> stop
                        FP(k) = FP(k) + 1;
                    end
                end
            end
        end
    end
end

Stage_FPRs  = FP/NumT               % stage-wise FPRs
% FPR         = sum(FP) / NumT        % total FPR
% end
%% adaptations I made
% for ti=1:NumT
%     
% %     Ps          = rand(1,K);            % random p values: all assumptions are met, so the FPR should be exact. 
%     Plog        = -2*log(Ps);           % fisher transform
%     Detected    = false;                % 
%     
%     % check for rejections
%     for k=1:K
%         if sum(Plog(1:k)) > Thresholds(k)       % H0 rejected at stage 2 --> stop
%             FP(k) = FP(k) + 1;
%         end
%     end
% end
TN      = zeros(1, K);      % number of false-positives
for ti=1:NumT
%     Ps = 
%     Ps          = rand(1,K);            % random p values: all assumptions are met, so the FPR should be exact. 
    Detected    = false;                % 
    
    % check for rejections
    k = 1;
    if sum(Plog(1:k)) < gThresholds(k)       % H0 rejected at stage 1 --> stop
        TN(k) = TN(k) + 1;
    else
        k = 2;
        if sum(Plog(1:k)) < gThresholds(k)       % H0 rejected at stage 2 --> stop
            TN(k) = TN(k) + 1;
        else
            k = 3;
            if sum(Plog(1:k)) < gThresholds(k)   % H0 rejected at stage 3 --> stop
                TN(k) = TN(k) + 1;
            else
                k = 4;
                if sum(Plog(1:k)) < gThresholds(k)   % H0 rejected at stage 3 --> stop
                    TN(k) = TN(k) + 1;
                else
                    k = 5;
                    if sum(Plog(1:k)) < gThresholds(k)   % H0 rejected at stage 3 --> stop
                        TN(k) = TN(k) + 1;
                    end
                end
            end
        end
    end
end

Stage_TNRs  = TN/NumT               % stage-wise FPRs

end