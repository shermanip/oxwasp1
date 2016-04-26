%FUNCTION: MIXTURE OF LINEAR REGRESSION LOG LIKELIHOOD
%PARAMETERS:
    %Y: column vector of responses
    %X: design matrix
    %beta: column matrix of regression parameters
    %variance_array: column vector of residual variance
    %mixing_array: column vector of mixing coefficients
function lnL = mixture_logLikelihood(Y,X,beta_array,variance_array,mixing_array)

    %get m (number of models) and n (number of mean-variance pairs)
    m = numel(mixing_array);
    n = numel(Y);
    %define matrix of likelihoods for each data (dim1) and each model(dim2)
    lnL = zeros(n,m);
    %for each model
    for i = 1:m
        %get a n vector of likelihoods for model i
        lnL(:,i) = linearReg_likelihood(Y,X,beta_array(:,i),variance_array(i));
    end
    %multiply the likelihood by the prior
    lnL = lnL .* repmat(mixing_array',n,1);
    %for each data, marginalise over all models, then sum the logs
    lnL = sum(log(sum(lnL,2)));

end