%load the data
folder_location = '/data/tinamou/sip/block_images/orginial'; %location of the data
[sample_mean,sample_var,area] = load_meanVariance(folder_location);

p_max = 6;
lnL_array = zeros(p_max,1);

for p = 1:p_max
    %number of polynomial orders
    m = p+1;

    %get the response vector
    Y = sample_var;

    w = 1./(Y.^2);
    w = (w/max(w));

    X = repmat(sample_mean,1,p).^repmat(1:p,area,1);
    x_mean = mean(X,1);
    x_std = std(X,[],1);
    y_mean = mean(sample_var);
    y_std = std(sample_var);

    Y = (Y-y_mean)/y_std;
    X = [ones(area,1), (X-repmat(x_mean,area,1))./repmat(x_std,area,1)];

    wX = X;
    for feature = 1:m
        wX(:,feature) = X(:,feature).*w;
    end

    beta = (X'*wX)\(X'*(w.*Y));
    MSE = sum(w.*((Y-X*beta).^2)) / sum(w);
    %estimate the covariance matrix
    cov = wX/(X'*wX);
    cov = (cov'*cov) * MSE;

    MSE = MSE * y_std^2;
    %un-normalise the covariance matrix

    cov_units = y_std ./ [1,x_std];
    cov_units = cov_units'*cov_units;
    cov = cov.*cov_units;
    beta(2:end) = y_std*beta(2:end)./x_std';
    beta(1) = y_mean + beta(1)*y_std - sum(x_mean'.*beta(2:end));

    plotHistogramHeatmap(sample_mean,sample_var);
    hold on;
    x_min = min(sample_mean);
    x_max = max(sample_mean);
    x_length = x_max-x_min;
    x_plot = (x_min:(x_length/1000):x_max)';
    X_plot = [ones(numel(x_plot),1),repmat(x_plot,1,p).^repmat(1:p,numel(x_plot),1)];

    y_hat = X_plot*beta;
    plot(X_plot(:,2),y_hat,'r--');

    y_se = sqrt(diag(X_plot*cov*X_plot')+MSE);
    plot(x_plot,y_hat,'r--');
    plot(x_plot,y_hat+norminv(0.975)*y_se,'r--');
    plot(x_plot,y_hat-norminv(0.975)*y_se,'r--');
    hold off;

    lnL_array(p) = -(area/2) * (log(2*pi) + log(MSE) + 1);
end

BIC = -2*lnL_array + (1:p_max)'*log(100);
figure;
plot(BIC);
xlabel('Polynomial order');
ylabel('BIC (nat)');