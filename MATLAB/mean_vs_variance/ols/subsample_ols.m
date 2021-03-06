clc;
clear all;
close all;

%SCRIPT: linear regression on subsample of SAMPLE, BACKGROUND and FOAM
%plots the lienar regressions

folder_location = '/data/tinamou/sip/block_images/orginial'; %location of the data

%array of images
[stack,length,area,n] = load_stack(folder_location);
%work out the sample mean
sample_mean = reshape(mean(stack,3),[],1);
%work out the sample std
sample_var = reshape(var(stack,0,3),[],1);

%for SAMPLE, BACKGROUND, FOAM subsamples, plot histogram heatmap and plot linear regression
%fit linear regression of subsample of the SAMPLE
    %beta_1 = regression vector parameter
    %n_1 = number of sub-sample
    %model_1 = linear model object
    %range_1 = [min,max] of the mean
[beta_1,n_1,model_1,range_1] = histogramHeatmap_fitStraightLine(stack(236:(236+1506),328:(328+1398),:),1000);
ylim([0.05E6,2.5E5]); %adjust the y axis of this histogram
%fit linear regression of subsample of the BACKGROUND
[beta_2,n_2,model_2,range_2] = histogramHeatmap_fitStraightLine(stack(32:(32+964),22:(22+190),:),200);
%fit linear regression of subsample of the FOAM
[beta_3,n_3,model_3,range_3] = histogramHeatmap_fitStraightLine(stack(1942:(1942+38),12:(12+1972),:),200);

% %ANOVA F STATISTICS
% beta_bar = (n_1*beta_1(2,1)+n_2*beta_2(2,1)+n_3*beta_3(2,1))/(n_1+n_2+n_3);
% within_variance = ((n_1-2)*beta_1(2,2)^2+(n_2-2)*beta_2(2,2)^2+(n_3-2)*beta_3(2,2)^2) / (n_1+n_2+n_3-6);
% between_variance = (n_1*(beta_1(2,1)-beta_bar)^2 + n_2*(beta_2(2,1)-beta_bar)^2 + n_3*(beta_3(2,1)-beta_bar)^2)/(n_1+n_2+n_3-1);
% F = between_variance / within_variance;
% fcdf(F,n_1+n_2+n_3-1,n_1+n_2+n_3-6,'upper');

%for all data, bin the data
plotHistogramHeatmap(sample_mean,sample_var);
hold on;

%plot linear regression of SAMPLE
x_min = range_1(1);
x_max = range_1(2);
x_length = x_max-x_min;
x_plot = x_min:(x_length/1000):x_max; %the range of means to plot over
%predict the variance given x_plot
[y_hat,y_error] = predict(model_1,x_plot','Prediction','observation');
%plot the mean and 95% confidence interval
plot(x_plot,y_hat,'r--');
plot(x_plot,y_error,'r--');

%plot linear regression of BACKGROUND
x_min = range_2(1);
x_max = range_2(2);
x_length = x_max-x_min;
x_plot = x_min:(x_length/1000):x_max; %the range of means to plot over
%predict the variance given x_plot
[y_hat,y_error] = predict(model_2,x_plot','Prediction','observation');
%predict the variance given x_plot
plot(x_plot,y_hat,'r--');
plot(x_plot,y_error,'r--');

%plot linear regression of FOAM
x_min = range_3(1);
x_max = range_3(2);
x_length = x_max-x_min;
x_plot = x_min:(x_length/1000):x_max; %the range of means to plot over
%predict the variance given x_plot
[y_hat,y_error] = predict(model_3,x_plot','Prediction','observation');
%predict the variance given x_plot
plot(x_plot,y_hat,'r--');
plot(x_plot,y_error,'r--');

hold off;