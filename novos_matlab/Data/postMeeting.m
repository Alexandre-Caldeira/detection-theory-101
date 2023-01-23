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
load Thresholds.mat
digits(32);

% Get thresholds
validate = 1;
K = 5;
% alphaVector = [0.002,0.002,0.002,0.002,0.002];
alphaVector = (0.01/K).*ones(1,K);
gammaVector = (0.09/K).*ones(1,K);
Resolution  = (1/0.0001);
NumT    = 100000;
[Thresholds,Stage_FPRs] = chesnayeThresholds(K,alphaVector,gammaVector,Resolution, validate, NumT);
alpha = Thresholds(1,:);

%%
% 8?   9?  
for SS = 1:12   
disp('-------------------------------------------------------')
sub = DB00{SS};
% sub1 = DB40{1};
% DB00{1};  subject 1 in exam without stimulus
start = 5;
samples = 75;
scaling = 0.001;
sub = scaling*sub(:,start:start+samples-1);
sub = sub - mean(sub,'all');
[Nsub,sampleSize] = size(sub);


% SUB1 = fft(sub1);
% T2 in time-domain out-performed freq.-domain T2

% V = NxQ;
% matrix V = NxQ => vij = ith epoch, jth feature
% where N is the number of epochs 
% and Q are features extracted from each epoch
Q = 25;
N = 3000;
starti = 1;
% starti = abs(Nsub-N);
V = nan(N,Q); 
for i = starti:1:starti+N-1
    for j = 3:3:sampleSize
        V(i-starti+1,j/3) = mean(sub(i,j-2:j));
        % xi => mean feature value down the ith column
        % "In particular, the 15 ms windows following stimuli onset
        % are split into 25 segments of equal duration (0.6 ms per 
        % segment), and the mean is taken across each segment, 
        % giving an Nx25-dimensional feature matrix V"
    end
end

% Q = 25;
% Nmax = 3000;
% V = nan(Nmax,Q); 
% diff = Nsub-Nmax;
% for i = 1:1:Nmax
%     for j = 6:6:sampleSize
%         V(i,j/6) = mean(sub(i+diff,j-5:j));
%         % xi => mean feature value down the ith column
%         % "In particular, the 15 ms win dows following stimuli onset
%         % are split into 25 segments of equal duration (0.6 ms per 
%         % segment), and the mean is taken across each segment, 
%         % giving an Nx25-dimensional feature matrix V"
%     end
% end

if isempty(find(isnan(V)~=0, 1))
    disp('Matrix properly constructed, no NaNs left.')
else 
    error('Matrix has NaN!') 
end


%%


% DESIRED_ALPHA = 0.15; % 0.01;
% gamma = [0.2, 0.4,0.25];
% p = hotT2(X);

% alpha = [12.4292,   16.0489,   19.1948,   22.0851,   24.7741];
% gamma = [0.2107,    1.6728,    4.4597,    8.9531,   24.7741];

K       = 5;                % test the 5-staged design
FP      = zeros(1, K);      % number of false-positives
% NumT    = 1000000;          % number of tests to carry out

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

%--------

% for idx = 1:size(Ps,2)
%     if Ps(idx)<0.001
%         Ps(idx) = 0.001;
%     end
% end
% Plog  = -2*log(Ps);
Plog = zeros(size(Ps));
SPlog = zeros(size(Ps));

% for ti=1:NumT
    
%     Ps          = rand(1,K);            % random p values: all assumptions are met, so the FPR should be exact. 
%     Plog        = -2*log(Ps);           % fisher transform
%     Detected    = false;                % 
    
    % check for rejections
    for k=1:K
%         Plog(k) = chi2inv(1-Ps(k),2);
        Plog(k) = -2*log(Ps(k));    
        SPlog(k) = sum(Plog(1:k));
%         SPlog(k) = sum(Plog(1:k));
%         SPlog(k) = SPlog(k-1)+Plog(k);
        if SPlog(k) > alpha(k)       % H0 rejected at stage 2 --> stop
            FP(k) = FP(k) + 1;
        end
    end
% end

% for ti=1:NumT
%     
% %     Ps          = rand(1,K);            % random p values: all assumptions are met, so the FPR should be exact. 
%     Plog        = -2*log(Ps);           % fisher transform
%     Detected    = false;                % 
%     
%     % check for rejections
%     k = 1;
%     if sum(Plog(1:k)) > Thresholds(k)       % H0 rejected at stage 1 --> stop
%         FP(k) = FP(k) + 1;
%     else
%         k = 2;
%         if sum(Plog(1:k)) > Thresholds(k)       % H0 rejected at stage 2 --> stop
%             FP(k) = FP(k) + 1;
%         else
%             k = 3;
%             if sum(Plog(1:k)) > Thresholds(k)   % H0 rejected at stage 3 --> stop
%                 FP(k) = FP(k) + 1;
%             end
%         end
%     end
% end

% Stage_FPRs  = FP/NumT               % stage-wise FPRs
% FPR         = sum(FP) / NumT        % total FPR

disp('Subject:')
disp(SS)
disp('Table 1:')
disp(Ps)
disp('Table 2:')
disp(SPlog)
disp('-------------------------------------------------------')
end