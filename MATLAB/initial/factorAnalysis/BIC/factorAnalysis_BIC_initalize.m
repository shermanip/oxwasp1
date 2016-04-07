%PROCEDURE: FACTOR ANALYSIS BIC - TO BE RUN LOCALLY
    %STEP 1: Initalize variables
%save variables to a binary file
    %k_max: the number of maximum factors to look at
    %n_EM: number of EM steps
    %lnL_array: array of maximum log likelihood for each k
    %BIC_optimal: the minimum BIC found
function factorAnalysis_BIC_initalize()

    k_max = 15;
    n_EM = 100;
    lnL_array = -inf*ones(k_max,1);
    BIC_optimal = inf;
    n_totalRun = 0;

    %save the variables
    save('/data/tinamou/sip/results/fa_bic.mat');
end