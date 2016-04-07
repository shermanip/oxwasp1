%FUNCTION FACTOR ANALYSIS BOOTSTRAP BIC - RUN ON GREYPARTRIDGE
    %STEP 1: initalize variables and save them in a binary file
        %k_max: maximum number of factors to investigate
        %n_EM: number of EM steps
        %n_run: number of initial values to run
        %lnL_array: n_bootstrap x k_max matrix of log likelihoods
        %n_bootstrap: number of bootstraps done
function factorAnalysis_bootstrapBIC_initalize()

    %number of factors
    k_max = 20;
    %number of EM steps
    n_EM = 100;
    %number of initial values
    n_run = 1;
    %array of log likelihoods
    lnL_array = zeros(0);
    %number of bootstraps done
    n_bootstrap = 0;

    %save the variables
    save('/data/greypartridge/oxwasp/oxwasp15/sip/fa_bootstrap.mat');

end