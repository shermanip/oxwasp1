%FUNCTION: FACTOR ANALYSIS INITALIZE VARIABLES FOR EM
%PARAMETERS:
    %p: number of features
    %k: number of factors
%RETURN:
    %loading: pxk loading matrix
    %noise_vector: p-vector of instrinic noise
function [loading,noise_vector] = factorAnalysis_EM_initalize(p,k)
    loading = normrnd(0,1E3,p,k); %loading
    noise_vector = gamrnd(1,1E3,p,1); %vector of noise
end

