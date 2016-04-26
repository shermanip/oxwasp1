%FUNCTION: LIKELIHOOD OF LINEAR REGRESSION
%PARAMETERS:
    %Y: column vector of responses
    %X: design matrix
    %beta: parameter vector
    %variance: the variance of residual%
%RETURN:
    %n vector of likelihoods (one for each data)
function l = linearReg_likelihood(Y,X,beta,variance)

    l = -0.5*log(2*pi) - 0.5*log(variance) - 0.5*(Y-X*beta).^2/variance;
    l = exp(l);
    
end

