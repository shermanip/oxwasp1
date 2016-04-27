%FUNCTION: GAUSSAN WEIGHTED LEAST SQUARES REGRESSION
%PARAMETER:
    %X: design matrix (nx2)
    %Y: response vector (n)
    %mu: mean (scalar)
    %width: gaussian width (scalar)
%RETURN:
    %beta: regression parameter vector (2)
    %cov: covariance matrix of beta (2x2)
    %MSE: weighted mean squared error (scalar)
function [beta,cov,MSE] = gaussianWeighted_ols(X,Y,mu,width)

    %check the parameters
    checkParameters(X,Y,mu,width);
    
    %calculate the weights
    w = exp(-0.5 * ( (X(:,2)-mu) / width ).^2 );
    %scale the weights so that the largest weight is 1
    w = (w/max(w));
    
    %work out diag(w)*X
    wX = X;
    for feature = 1:2
        wX(:,feature) = X(:,feature).*w;
    end

    %do weighted least squares
    beta = (X'*wX)\(X'*(w.*Y));
    %work out the weighted mean squared error
    MSE = sum(w.*((Y-X*beta).^2)) / sum(w);
    %estimate the covariance matrix
    cov = wX/(X'*wX);
    cov = (cov'*cov) * MSE;
    
    %NESTED FUNCTION: TO CHECK PARAMETERS WITH
    function checkParameters(X,Y,mu,width)
        %check the parameters, display error message otherwise
        if ~iscolumn(Y) %check if column vector
            error('Error in gaussianWeighted_ols(X,Y,mu,width), Y is not a column vector');
        end
        if ~isscalar(mu) %check if scalar
            error('Error in gaussianWeighted_ols(X,Y,mu,width), mu is not a scalar');
        end
        if ~isscalar(width) %check if scalar
            error('Error in gaussianWeighted_ols(X,Y,mu,width), width is not a scalar');
        end
        [n,p] = size(X);
        %check if there are 2 columns in X
        if p~=2
            error('Error in gaussianWeighted_ols(X,Y,mu,width), number of columns in X is not 2');
        end
        %check if the number of rows in X is the same as the size of Y
        if n~=numel(Y)
            error('Error in gaussianWeighted_ols(X,Y,mu,width), number of rows in X is not the same as the number of elements in Y');
        end
    end
end

