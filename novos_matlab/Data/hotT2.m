function p = hotT2(X)
digits(32);
%hotT2 Summary  Hotelling'S T-Squared statistic
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

    [N,Q] = size(X);
    mu=zeros([1,Q]); % EXPECTED MEAN
    m=mean(X); %Mean vector from data matrix X.
    S=cov(X);  %Covariance matrix from data matrix X.
    
    T2=N*(m-mu)*pinv(S)*(m-mu)'; %Hotelling's T-Squared statistic.
%     T2=N*(m-mu)/S*(m-mu)'; % exact method
%     T2=N*((m-mu)/S)*ctranspose(m-mu); % exact method

    % "The T2 statistic can then be transformed into an F 
    % statistic using:

    F = T2*(N-Q)/(Q*(N-1));
    
    % which follows an F-distribution with Q and N – Q degrees of
    % freedom (DOF) under H0. It is worth noting that the number 
    % of epochs N should be larger than the number of features Q, 
    % else S1 cannot be calculated.

    v1=Q-1;  %Numerator degrees of freedom.
    v2=N-Q;  %Denominator degrees of freedom.
    p=fcdf(F,v1,v2);  %Probability that null Ho: is true
end

