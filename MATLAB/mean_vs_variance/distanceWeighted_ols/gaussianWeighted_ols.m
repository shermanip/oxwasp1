function [beta,cov,MSE] = gaussianWeighted_ols(X,Y,mu,width)

    w = exp(-0.5 * ( (X(:,2)-mu) / width ).^2 );
    w = (w/max(w));
    
    wX = X;
    for feature = 1:2
        wX(:,feature) = X(:,feature).*w;
    end

    beta = (X'*wX)\(X'*(w.*Y));
    MSE = sum(w.*((Y-X*beta).^2)) / sum(w);
    cov = wX/(X'*wX);
    cov = (cov'*cov) * MSE;
    
end

