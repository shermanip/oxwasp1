%FUNCTION: MIXTURE OF LINEAR REGRESSION LOG LIKELIHOOD
%PARAMETERS:
    %Y: column vector of responses (n)
    %X: design matrix (nx2)
    %beta: column matrix of regression parameters (2xm)
    %variance_array: column vector of residual variance (m)
    %mixing_array: column vector of mixing coefficients (m)
%RETURN:
    %log likelihood (scalar) of the mixture of lienar regression
function lnL = mixture_logLikelihood(Y,X,beta_array,variance_array,mixing_array)

    %check parameters are correct
    %get m (number of models) and n (number of mean-variance pairs)
    m = 0;
    n = 0;
    checkParameters(Y,X,beta_array,variance_array,mixing_array);

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

    %NESTED FUNCTION TO CHECK PARAMETERS
    %check the parameters, display error message otherwise
    function checkParameters(Y,X,beta_array,variance_array,mixing_array)
        
        %CHECK DATA TYPE
        if ~iscolumn(Y) %check if response vector is a column vector
            error('Error in mixture_logLikelihood(Y,X,beta_array,variance_array,mixing_array), Y is not a column vector');
        end
        if ~iscolumn(variance_array) %check if variance_array is a column vector
            error('Error in mixture_logLikelihood(Y,X,beta_array,variance_array,mixing_array), variance_array is not a column vector');
        end
        if ~iscolumn(mixing_array) %check if mixing_array is a column vector
            error('Error in mixture_logLikelihood(Y,X,beta_array,variance_array,mixing_array), mixing_array is not a column vector');
        end
        
        %CHECK MATRIX SIZE
        [n,p_x] = size(X); %get size of the design matrix
        [p_b,m] = size(beta_array); %get size of the column matrix
        %check if there are 2 columns in X
        if p_x~=2
            error('Error in mixture_logLikelihood(Y,X,beta_array,variance_array,mixing_array), number of columns in X is not 2');
        end
        %check if there are 2 rows in beta_array
        if p_b~=2
            error('Error in mixture_logLikelihood(Y,X,beta_array,variance_array,mixing_array), number of columns in X is not 2');
        end
        %check if the number of rows in X is the same as the size of Y
        if n~=numel(Y)
            error('Error in mixture_logLikelihood(Y,X,beta_array,variance_array,mixing_array), number of rows in X is not the same as the number of elements in Y');
        end
        %check if the number of columns in beta_array is the same as the size of variance_array
        if m~=numel(variance_array)
           error('Error in mixture_logLikelihood(Y,X,beta_array,variance_array,mixing_array), number of columns in beta_array is not the same as the number of elements in variance_array');
        end
        %check if the number of columns in beta_array is the same as the size of mixing_array
        if m~=numel(mixing_array)
            error('Error in mixture_logLikelihood(Y,X,beta_array,variance_array,mixing_array), number of columns in beta_array is not the same as the number of elements in mixing_array');
        end
    end
    
end