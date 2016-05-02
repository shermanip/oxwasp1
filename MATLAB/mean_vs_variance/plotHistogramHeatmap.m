%PROCEDURE: PLOT HISTOGRAM HEAT MAP OF THE MEAN VS VARIANCE RELATIONSHIP
    %note: this function is for ALL the mean-variance pairs#
%PARAMETERS:
    %sample_mean: n vector of the sample mean grey values
    %sample_var: n vector of the sample variance grey values
function plotHistogramHeatmap(sample_mean,sample_var,nbin)

    checkParameters(sample_mean,sample_var);

    %number of marginal bins to bin the data into
    if nargin==2
        nbin = 3000;
    end
    %bin the data
    [N,c] = hist3([sample_var,sample_mean],[nbin,nbin]);
    
    %plot the heatmap
    figure;
    %normalize N so that the colormap is the frequency density
    imagesc(cell2mat(c(2)),cell2mat(c(1)),N/( (c{2}(2)-c{2}(1))*(c{1}(2)-c{1}(1)) ) );
    axis xy; %switch the y axis
    ylim([0.05E6,4E5]); %set the y limit
    colorbar; %display the colour bar
    xlabel('Sample grey value mean (AU)'); %label the axis
    ylabel('Sample grey value variance (AU^{2})'); %label the axis
    
    %NESTED FUNCTION: CHECK PARAMETERS
    function checkParameters(sample_mean,sample_var)
        %check if sample_mean is a column vector, if not throw
        if ~iscolumn(sample_mean)
            error('Error in plotHistogramHeatmap(sample_mean,sample_var), sample_mean is not a column vector');
        end
        %check if sample_var is a column vector, if not throw
        if ~iscolumn(sample_var)
            error('Error in plotHistogramHeatmap(sample_mean,sample_var), sample_var is not a column vector');
        end
        %check if sample_mean and sample_var has the same length, else throw
        n1 = numel(sample_mean);
        n2 = numel(sample_var);
        if n1~=n2
            error('Error in plotHistogramHeatmap(sample_mean,sample_var), sample_mean and sample_var are not the same length');
        end
    end

end

