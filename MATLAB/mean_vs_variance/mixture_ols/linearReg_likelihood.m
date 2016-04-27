%FUNCTION: LIKELIHOOD OF LINEAR REGRESSION
%PARAMETERS:
    %Y: column vector of responses (n)
    %X: design matrix (nx2)
    %beta: parameter vector (2)
    %variance: the variance of residual (scalar)
%RETURN:
    %n vector of likelihoods (one for each data)
function l = linearReg_likelihood(Y,X,beta,variance)

    %check parameters
    checkParameters(Y,X,beta,variance);

    %for each data, calculate the likelihood given the parameters
    l = -0.5*log(2*pi) - 0.5*log(variance) - 0.5*(Y-X*beta).^2/variance;
    l = exp(l);
    
    %NESTED FUNCTION TO CHECK PARAMETERS
    %check the parameters, display error message otherwise
    function checkParameters(Y,X,beta,variance)
        
        %CHECK DATA TYPE
        if ~iscolumn(Y) %check if column vector
            error('Error in linearReg_likelihood(Y,X,beta,variance), Y is not a column vector');
        end
        if ~iscolumn(beta) %check if column vector
            error('Error in linearReg_likelihood(Y,X,beta,variance), beta is not a column vector');
        end
        if ~isscalar(variance) %check if scalar
            error('Error in linearReg_likelihood(Y,X,beta,variance), variance is not a scalar');
        end
        
        %CHECK MATRIX SIZE
        if numel(beta)~=2 %check if beta is a 2-vector
            error('Error in linearReg_likelihood(Y,X,beta,variance), beta does not have 2 elements');
        end
        [n,p] = size(X);
        %check if there are 2 columns in X
        if p~=2
            error('Error in linearReg_likelihood(Y,X,beta,variance), number of columns in X is not 2');
        end
        %check if the number of rows in X is the same as the size of Y
        if n~=numel(Y)
            error('Error in linearReg_likelihood(Y,X,beta,variance), number of rows in X is not the same as the number of elements in Y');
        end
    end
    
end

