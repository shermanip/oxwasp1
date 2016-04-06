clc;
clear all;
close all;

%STEP 1: initalize variables and save it

%number of factors
k_max = 10;
%number of EM steps
n_EM = 100;
%array of log likelihoods
lnL_array = -inf*ones(k_max,1);
%BIC_optimal
BIC_optimal = inf;
%store the total number of runs
n_totalRun = 0;

%save it
save('fa_bic.mat');