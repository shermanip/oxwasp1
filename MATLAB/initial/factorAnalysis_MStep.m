%FUNCTION: FACTOR ANALYSIS M STEP
%PARAMETERS:
    %X: n x p design matrix
    %Y: n x k design matrix of the factors
    %Y_cov: k x k covariance matrix of the factors
    %n: sample size
    %p: number of features
%RETURN:
    %loading: p x k matrix
    %noise_vector: p vector of noise
function [loading,noise_vector] = factorAnalysis_MStep(X,Y,Y_cov,n,p)
    
    %declare empty vector
    noise_vector = zeros(p,1);

    %estimate the loading matirx
    loading = X'* (Y / (n*Y_cov+Y'*Y));
    
    %estimate the noise
    residual = X' - loading*Y';
    for i = 1:p
        lambda = loading(i,:)';
        noise_vector(i) = lambda'*Y_cov*lambda + sumsqr(residual(i,:))/n;
    end

end

