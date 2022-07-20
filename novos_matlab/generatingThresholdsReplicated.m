% clearvars;clc;close all
K               = 5;                            % analyse data K times
TotalAlpha      = 0.01;                         % alpha for the full test
alphaVector     = ones(1,K)*(TotalAlpha/K); 
gammaVector     = [0.1, 0.15, 0.2, 0.25, 0.29];
Thresholds      = chesnayeThresholds(K,alphaVector,gammaVector);
disp(Thresholds)
save('5Thresholds.mat','Thresholds')