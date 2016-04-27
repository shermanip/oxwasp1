%FUNCTION: GAUSSAN WEIGHTED LEAST SQUARES REGRESSION
%PARAMETER:
    %X: design matrix
    %Y: response vector
    %mu: mean
    %width: gaussian width
%RETURN:
    %beta: regression parameter
    %cov: covariance matrix of beta
    %MSE: weighted mean squared error
function [beta,cov,MSE] = gaussianWeighted_ols(X,Y,mu,width)

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
    
end

