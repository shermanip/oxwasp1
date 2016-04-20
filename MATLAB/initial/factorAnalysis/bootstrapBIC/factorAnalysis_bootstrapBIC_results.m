%FUNCTION: FACTOR ANALYSIS BOOTSTRAP BIC
    %STEP 3: RESULTS
    %displays the BIC for all k and bootstraps
    %display the mean and std of the optimial k
function factorAnalysis_bootstrapBIC_results

    %load the variables from the binary file
        %load the variables and past results
        %k_max: maximum number of factors to investigate
        %n_EM: number of EM steps
        %n_run: number of initial values to run
        %lnL_array: n_bootstrap x k_max matrix of log likelihoods
        %n_bootstrap: number of bootstraps done
    load('/data/greypartridge/oxwasp/oxwasp15/sip/fa_bootstrap.mat');
    
    %dataset properties
    length = 100;
    area = length^2;
    n = 100;

    %work out the BIC
    BIC = -2*lnL_array + repmat(log(n)*area*((1:k_max)+1),n_bootstrap,1);

    %box plot the BIC for each k
    figure;
    boxplot(BIC);
    xlabel('Number of factors');
    ylabel('BIC (nat)');

    %find the optimal k for each bootstrap
    [~,k_optimal] = min(BIC,[],2);
    %set display to scientific notation
    format longE;
    %display the 95% confidence interval of optimial k
    disp(strcat(num2str(n_bootstrap),' bootstraps done'));
    disp('mean and std of optimial k');
    disp(mean(k_optimal));
    disp('+/-');
    disp(std(k_optimal));
end