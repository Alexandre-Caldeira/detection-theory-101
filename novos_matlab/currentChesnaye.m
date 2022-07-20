% Referencias:
% Sci-Hub | The Convolutional Group Sequential Test: reducing test time for evoked potentials. IEEE Transactions on Biomedical Engineering, 1–1 | 10.1109/tbme.2019.2919696
% https://sci-hub.se/10.1109/tbme.2019.2919696
% 
% Sci-Hub | A group sequential test for ABR detection. International Journal of Audiology, 1–10 | 10.1080/14992027.2019.1625486
% https://sci-hub.se/10.1080/14992027.2019.1625486

%% cleanup
close all;clearvars;clc;

%% loading data

% ! Ensure data ins on Matlab's path
load Subject_ABR_Data.mat
% for SS = 1:12
sub1 = DB00{1};
% sub1 = DB40{1};
% DB50{1};  subject 1 in exam without stimulus
[N,sampleSize] = size(sub1);

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

%% Hotelling'S T-Squared statistic
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



nK = 3; %numero de etapas
max = 600; % floor(N/nK);

step = floor(max/nK);
numDetec = zeros(size(1:26:max-25));

% pHist = zeros(1,max);
% summary = zeros(1,max);
% phiSum = zeros(1,max); 
% phiSum = cell(1,nK);
numTests = 0;
pHistSem = [];
pHist = [];
summary = [];


% K = 0:max:N-max;
K = 1:max:3000;
DESIRED_ALPHA = 0.15; % 0.01;
gamma = [0.2, 0.4,0.25];
TAXA_POR_K = [];
phiHist = [];

for i = 1:numel(K)
    % each stage gets disjoint epochs
    X = V(K(i):K(i)+max-1,:);
    [N,Q] = size(X);

    alpha = DESIRED_ALPHA/nK;

    mu=zeros([1,Q]); % EXPECTED MEAN
    m=mean(X); %Mean vector from data matrix X.
    S=cov(X);  %Covariance matrix from data matrix X.
    
    T2=N*(m-mu)*pinv(S)*(m-mu)'; %Hotelling's T-Squared statistic.
%     T2=N*(m-mu)/S*(m-mu)'; % exact method

    % "The T2 statistic can then be transformed into an F 
    % statistic using:

    F = T2*(N-Q)/(Q*(N-1));
    
    % which follows an F-distribution with Q and N – Q degrees of
    % freedom (DOF) under H0. It is worth noting that the number 
    %of epochs N should be larger than the number of features Q, 
    % else S1 cannot be calculated.

    v1=Q;  %Numerator degrees of freedom.
    v2=N-Q;  %Denominator degrees of freedom.
    P=1-fcdf(F,v1,v2);  %Probability that null Ho: is true.

    pHist = [pHist,P];
    
    phiHist = [chi2inv(1-P,2),phiHist];
    
    
    


    summary = [summary,sum(phiHist)];
    
    if i==1
        xx = 0:0.0001:20;
        phi = chi2pdf(xx,i+1);
        cummPhi = cumtrapz(xx,phi);
        
        idxA = find(1-cummPhi<alpha,1);
        A = xx(idxA)
        
        idxC = find(cummPhi>gamma(i),1);
        C = xx(idxC)
    else
        oldPhi = phi;
        
%         xx = xx(1:idxA);
%         phi = conv(oldPhi(1:idxA),chi2pdf(xx,i+1),'same');
%         xx = xx(idxC:idxA);
        phi = conv(chi2pdf(xx,i+1),phi,'same');
        cummPhi = cumtrapz(xx,phi);
        
        idxA = find(1-cummPhi<alpha,1);
        A = xx(idxA)
        idxC = find(cummPhi>gamma(i),1);
        C = xx(idxC)

        
%         phi = chi2pdf(xx,i+1);

        
        
%         aux = find(1-cummPhi<alpha,1);
        
        % 1-cummPhi(aux) == alpha
        
%         C = xx(find(cummPhi>gamma,1));

    end
    figure(i)
    plot(xx,phi)
    drawnow
    
    
    
    if P >= alpha
%         disp('H0 cannot be rejected')
        % No signal
        numDetec(i) = 0;
%         numDetec = [numDetec,0];
    else
        disp('HIT! H0 rejected, signal detected.')
        % Mean vectors results significant.
        numDetec(i) = 1;
%         numDetec = [numDetec,1];
    end
    
%     quantile(phiSum{2},alpha)
    
    numTests = numTests+1;
    TAXA_POR_K = [TAXA_POR_K, sum(numDetec)/numTests];
end
% end

disp('Detection rate:')
disp(sum(numDetec)/numTests)
disp('')
disp('Rejection rate: (FPR ou TPR)')
disp(1-sum(numDetec)/numTests)
summary
pHist

%% too afraid to delete these ideias:
%     fcdf(F,v1,v2,'upper')

%     pHistSem = [pHistSem,P];
%     pHist = [pHist,P2];

%     summary(K) = sum(pHist);
%     if i>1
%         phi = [chi2inv(1-P,2+i), phi];
% %         phiSum(i) = {conv(phiSum{i-1},fpdf(X(:,randi(Q,1,Q)),v1,v2))};
%     else
%         phi = chi2inv(1-P,2+i);
% %         phiSum(1) = {fpdf),v1,v2)};
%     end

