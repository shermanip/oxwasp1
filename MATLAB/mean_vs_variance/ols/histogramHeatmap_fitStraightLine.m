%FUNCTION: PLOT HISTOGRAM HEAT MAP AND FIT LINEAR REGRESSION
%PARAMETERS:
    %stack: an array of images in matrix form (3 dimensions)
    %nbin: number of bins to bin the data in (scalar)
%RETURN:
    %coefficients: [mdl.Coefficients.Estimate,mdl.Coefficients.SE]
    %n: sample size
    %mdl: model object
    %mean_range: [min(sample_mean),max(sample_mean)]
function [coefficients,n,mdl,mean_range] = histogramHeatmap_fitStraightLine(stack,nbin)
    
    %check the parameters
    checkParameters(stack,nbin);

    %work out the sample mean
    sample_mean = reshape(mean(stack,3),[],1);
    %work out the sample std
    sample_var = reshape(var(stack,0,3),[],1);

    %get the number of data points
    area = numel(sample_mean);
    %fit linear model, mean vs variance
    mdl = LinearModel.fit(sample_mean,sample_var);

    %get the range of the means
    x_min = min(sample_mean);
    x_max = max(sample_mean);
    x_length = x_max-x_min;
    x_plot = x_min:(x_length/nbin):x_max;
    %given the range of means, plot the predicted variance
    [y_hat,y_error] = predict(mdl,x_plot','Prediction','observation');
    
    %bin the data
    [N,c] = hist3([sample_var,sample_mean],[nbin,nbin]);
    
    %plot the histogram as a heat map
    figure;
    %normalize N so that the colormap is the frequency density
    imagesc(cell2mat(c(2)),cell2mat(c(1)),N/( (c{2}(2)-c{2}(1))*(c{1}(2)-c{1}(1)) ) );
    hold on;
    axis xy; %switch the y axis
    colorbar; %display the colour bar
    xlabel('Sample grey value mean (AU)');
    ylabel('Sample grey value variance (AU^{2})');
    %plot the linear regression
    plot(x_plot,y_hat,'r--');
    plot(x_plot,y_error,'r--');
    hold off;

    %return the regression coefficient
    coefficients = [mdl.Coefficients.Estimate,mdl.Coefficients.SE];
    %return the sample size
    n = area;
    %return the min and max of the means
    mean_range = [min(sample_mean),max(sample_mean)];
    
    %NESTED FUNCTION: check the parameters are of the correct tye
    function checkParameters(stack,nbin)
        %check stack is an array of matrices
        if ndims(stack)~=3
            error('Error in histogramHeatmap_fitStraightLine(stack,nbin), stack is not an array of matrices');
        end
        %check nbin is a scalar
        if ~isscalar(nbin)
            error('Error in histogramHeatmap_fitStraightLine(stack,nbin), nbin is not a scalar');
        end
    end
end

