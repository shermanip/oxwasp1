%FUNCTION: FACTOR ANALYSIS E STEP
%PARAMETERS:
    %loading: p x k design matrix
    %noise_vector: p vector of the noise
    %X: n x p design matrix
    %k: number of factors
%RETURN:
    %Y: n x k design matrix of the factors
    %Y: k x k covariance matrix of the factors
function [Y,Y_cov] = factorAnalysis_EStep(loading,noise_vector,X,k)
    
    noise_inv = diag(1./noise_vector);
    Y = ((loading'*noise_inv*loading+eye(k)) \ (loading'*noise_inv*X'))';
    Y_cov = inv(eye(k)+loading'*noise_inv*loading);

end