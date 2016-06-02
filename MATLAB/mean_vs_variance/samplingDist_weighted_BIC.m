%FUNCTION: SAMPLING DISTRIBUTION WEIGHTED LEAST SQUARES REGRESSION
%This fits a polynomial on the mean vs variance data and plot the fit and
%returns the BIC
%PARAMETERS:
    %sample_mean: column vector of sample means
    %sample_var: column vector of sample variances
    %wantPlot: binary 1 if want plots of the fitted regression
    %p_max: the maxmimum polynomial order to fit
%RETURN:
    %BIC: row vector the BIC for each polynomial order
function BIC = samplingDist_weighted_BIC(sample_mean,sample_var,wantPlot,p_max,n_bin)

    %check the parameters are of the correct type
    checkParameters(sample_mean,sample_var,wantPlot,p_max);

    %get the number of pixels
    area = numel(sample_mean);
    
    %define array for log likelihood for each polynomial order
    lnL_array = zeros(p_max,1);

    %get the response vector
    Y = sample_var;

    %set the weights
    w = 1./(Y.^2);
    w = (w/max(w)); %normalise the weight so the maximum is 1
    
    %normalize the response vector to have zero mean and std 1
    y_mean = mean(sample_var);
    y_std = std(sample_var);
    Y = (Y-y_mean)/y_std;
    
    %for each polynomial order
    for p = 1:p_max
        
        %m = number of polynomial orders
        m = p+1;
        
        %define design matrix with polynomial features
            %rows: for each data
            %columns: for each polynomial feature
        X = repmat(sample_mean,1,p).^repmat(1:p,area,1);
        %normalize the design matrix so that the columns have mean zero andstd 1
        x_mean = mean(X,1);
        x_std = std(X,[],1);
        
        %append constant of 1
        X = [ones(area,1), (X-repmat(x_mean,area,1))./repmat(x_std,area,1)];

        %calculate diag(w)*X without creating diag(w) itself
        wX = X;
        for feature = 1:m
            wX(:,feature) = X(:,feature).*w;
        end

        %estimate the weighted regression parameter
        beta = (X'*wX)\(X'*(w.*Y));
        %calculate the weighted mean squared error
        MSE = sum(w.*((Y-X*beta).^2)) / sum(w);
        %estimate the covariance matrix
        cov = wX/(X'*wX);
        cov = (cov'*cov) * MSE;
        
        %un-normalise the mean squared error
        MSE = MSE * y_std^2;
        %un-normalise the covariance matrix
        cov_units = y_std ./ [1,x_std];
        cov_units = cov_units'*cov_units;
        cov = cov.*cov_units;
        %un-normalise the regression parameter
        beta(2:end) = y_std*beta(2:end)./x_std';
        beta(1) = y_mean + beta(1)*y_std - sum(x_mean'.*beta(2:end));

        %if want the plot...
        if wantPlot
            
            %plot the histogram
            ax = plotHistogramHeatmap(sample_mean,sample_var,n_bin);
            %resize the histogram
            set(ax, 'Position', [600,800,400,300]);
            %set the label of the histogram
            xlabel('Sample mean (arb. unit)');
            ylabel('Sample variance {(arb. unit^2)}');
            %hold on to plot more things on it
            hold on;
            
            %get range of x to plot the regression
            x_min = min(sample_mean);
            x_max = max(sample_mean);
            x_length = x_max-x_min;
            x_plot = (x_min:(x_length/1000):x_max)';
            %given the range of x, create the design matrix X_plot
            X_plot = [ones(numel(x_plot),1),repmat(x_plot,1,p).^repmat(1:p,numel(x_plot),1)];

            %predict the response given the x
            y_hat = X_plot*beta;
            %plot it
            plot(X_plot(:,2),y_hat,'r--');
            %work out the 90% prediction interval
            y_se = sqrt(diag(X_plot*cov*X_plot')+MSE);
            plot(x_plot,y_hat,'r--');
            plot(x_plot,y_hat+norminv(0.975)*y_se,'r--');
            plot(x_plot,y_hat-norminv(0.975)*y_se,'r--');
            hold off;
        end

        %work out the log likelihood
        lnL_array(p) = -(area/2) * (log(2*pi) + log(MSE) + 1);
    end

    %work out the BIC for each parameter
    BIC = (-2*lnL_array + ((1:p_max)+1)'*log(area))';
    
    %NESTED FUNCTION: CHECK PARAMETERS are of the correct type
    function checkParameters(sample_mean,sample_var,wantPlot,p_max)
        
        %check if sample_mean is a column vector, if not throw
        if ~iscolumn(sample_mean)
            error('sample_mean is not a column vector');
        end
        
        %check if sample_var is a column vector, if not throw
        if ~iscolumn(sample_var)
            error(' sample_var is not a column vector');
        end
        
        %check if sample_mean and sample_var has the same length, else throw
        n1 = numel(sample_mean);
        n2 = numel(sample_var);
        if n1~=n2
            error('sample_mean and sample_var are not the same length');
        end
        
        %check if wantPlot is a scalar
        if (~isscalar(wantPlot))
            error('wantPlot is not a scalar');
        %else check if wantPolt is either 0 or 1
        elseif ~ismember(wantPlot,[0,1])
            error('wantPlot has to be either 0 or 1');  
        end
        
        %check if p_max is a positive integer
        if ( (~isscalar(p_max)) || (p_max <= 0) || (floor(p_max)~=p_max) )
            error('p_max must be a positive integer');
        end
    end

end

