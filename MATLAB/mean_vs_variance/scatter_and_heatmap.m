clc;
clear all;
close all

%get the mean-variance data
folder_location = '/data/tinamou/sip/block_images/orginial'; %location of the data
[sample_mean,sample_var,area] = load_meanVariance(folder_location);

%plot the histogram as a heatmap
plotHistogramHeatmap(sample_mean,sample_var);

%set random seed
rng(21236878);
%partition the data
c = cvpartition(numel(sample_mean),'KFold',500);
%scatter plot selected partitions
for i = 37:38
    figure('Position', [600,800,400,300]);
    scatter(sample_mean(c.test(i)),sample_var(c.test(i)));
    xlabel('Sample mean (arb. unit)'); %label the axis
    ylabel('Sample variance {(arb. unit^2)}'); %label the axis
end
%display the number of samples in the partitions
disp('Number of samples');
disp(c.TestSize([37,38]));
