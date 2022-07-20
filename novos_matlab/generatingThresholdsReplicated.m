K               = 5;                            % analyse data 3 times
TotalAlpha      = 0.05;                         % pre-specified alpha-level for the full sequential test
alphaVector     = ones(1,K)*(TotalAlpha/K);     % alpha to spend at each stage
gammaVector     = [0.1, 0.15, 0.2, 0.25, 0.29];
Thresholds      = chesnayeThresholds(K,alphaVector);
disp(Thresholds)
save('5Thresholds.mat','Thresholds')