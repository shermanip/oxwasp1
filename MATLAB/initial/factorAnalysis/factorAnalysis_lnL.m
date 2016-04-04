%FUNCTION FACTOR ANALYSIS LOG LIKELIHOOD
%PARAMETERS:
    %loading: loading matrix lambda
    %noise: diagional matrix phi
    %cov_x: estimated covariance matrix
    %n: sample size
%RETURN:
    %lnL: log likelihood
function lnL = factorAnalysis_lnL(loading,noise_vector,cov_x,n,p)

    %covariance = loading*loading' + diag(noise)
    covariance = loading*loading';
    for i = 1:p
        covariance(i,i) = covariance(i,i) + noise_vector(i);
    end

    %work out the ln_det using the diagional of the Cholesky decomposition
    L = chol(covariance);
    ln_det = 2*sum(log(abs(diag(L))));

    %return the ln likelihood
    lnL = -0.5*n*ln_det - 0.5*n*p*log(2*pi) - 0.5*trace(covariance\cov_x);
    
end