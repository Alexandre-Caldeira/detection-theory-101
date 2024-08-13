%% Setup
clearvars; close all; clc

%% From
% file:///C:/PPGEE/Assessing%20CGST%20on%20ASSR/main_papers/chesnaye2019_cgst_theory.pdf#page=8.50
% The Convolutional Group Sequential Test: reducing
% test time for evoked potentials
% Michael A. Chesnaye, Steven L. Bell, James M. Harte, & David M. Simpson
% 

% teste by Alexandre Caldeira PPGEE UFMG 12/08/24

%% Parameters (pg. 6)
% The total α-level per dB SL condition was set to 0.01, which was also 
% split equally across the 5 stages, i.e. α1 = α2 = α3 = α4 = α5 = 0.002. 
% For each dB SL condition, the fraction of tests rejected for futility was
% set to 0.1, 0.15, 0.2, 0.25, and 0.29 for stages 1, 2, 3, 4, and 5 
% respectively, i.e. γ1 = 0.1, γ2 = 0.15, γ3 = 0.2, γ4 = 0.25, and 
% γ5 = 0.29 ...
% ... the following thresholds for efficacy:
A1 = 12.429; A2 = 16.049; A3 = 19.195; A4 = 22.085; A5 = 24.774;
% along with the following thresholds for futility:
C1 = 0.211; C2 = 1.673; C3 = 4.46; C4 = 8.953; C5 = A5; % = 24.774. 
% At each stage of the analysis, the sub-sample of 600 epochs was analysed
% using the Hotellings T 2test.

A = [A1,A2,A3,A4,A5];
C = [C1,C2,C3,C4,C5];

