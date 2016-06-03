clc;
clear all;
close all;

%SCRIPT: (multiple) gaussian weighted regression
    %the regressions are centered according to kmeans on the means
    %plot gradient vs gaussian width
    %plot weighted mse vs gaussian width
    %plot the weighted regression for different widths

%USER DEFINED VARIBLES
rng(118055149); %set random seed
K = 3; %number of regressions
folder_location = '/data/tinamou/sip/block_images/orginial'; %location of the data
width_array = 2:0.1:6; %gaussian widths to investigate
width_select = 10.^[2,3,4,5]; %gaussian widths to plot

%PROGRAM STARTS
%load the mean-variance pair data
[sample_mean,sample_var,area] = load_meanVariance(folder_location);

%normalize the data, estimate the mean and std of the mean and variance
mean_Y = mean(sample_var); %mean of the variances
std_Y = std(sample_var); %std of the variances
mean_X = mean(sample_mean); %mean of the means
std_X = std(sample_mean); %std of the means

%normalize the data to have 0 mean, std 1
Y = (sample_var-mean_Y)/std_Y;
X = [ones(area,1),(sample_mean-mean_X)/std_X]; %X is a design matrix

%get the centroids of the sample_means using K-means
[~,k_means] = kmeans(sample_mean,K);
%sort the centroids in accending order
k_means = sort(k_means);
%normalise the centroids
k_means = (k_means - mean_X)/std_X;

%declare matrix containing gradient and its standard error
    %dim 1: for each gaussian width
    %dim 2: estimate, standard error
gradient_array = zeros(numel(width_array),2);
%declare array to store the mse for each guassian width
mse_array = zeros(numel(width_array),1);
%for each entroid
for k = 1:K
    %for each width to be investigated
    for i = 1:numel(width_array)
        %get the normalized width
        width = (10^(width_array(i)))/std_X;
        %do weighted least squares, get the estimate and covariance matrix
        [beta,beta_cov,mse] = gaussianWeighted_ols(X,Y,k_means(k),width);
        %save the estimate of the gradient
        gradient_array(i,1) = beta(2);
        %save the estimate of the standard error of the gradient
        gradient_array(i,2) = sqrt(beta_cov(2,2));
        %save the estimate of the mse
        mse_array(i) = mse;
    end
    
    %un-normalise the gradient and the standard error
    gradient_array = gradient_array*std_Y/std_X;
    
    %plot width vs gradient with the standard error as the error bar
    figure;
    errorbar(width_array,gradient_array(:,1),gradient_array(:,2));
    xlim([min(width_array)-0.05,max(width_array)+0.05]);
    xlabel('log_{10} Gaussian width (log_{10}AU)');
    ylabel('Gradient (AU)');
    
    %plot width vs mse
    figure;
    plot(width_array,mse_array);
    xlabel('log_{10} Gaussian width (log_{10}AU)');
    ylabel('Weighted MSE (AU^{2})');
end

%declare column matrix for storing regression parameters for each cluster
beta_array = zeros(2,K);
%declare array of matrices storing the covariance matrix
beta_cov_array = zeros(2,2,K);
%declare K vector storing the mse
mse_array = zeros(K,1);

%for each guassian to plot
for j = 1:numel(width_select)
    %get the normalised width
    width = width_select(j)/std_X;
    
    %for each centroid
    for i = 1:K
        %get the estimate, covariance and mse of the weighted least squares
        [beta,beta_cov,mse] = gaussianWeighted_ols(X,Y,k_means(i),width);
        %save the regression parameter
        beta_array(:,i) = beta;
        %save the covariance matrix
        beta_cov_array(:,:,i) = beta_cov;
        %save the mse
        mse_array(i) = mse;
    end
    
    %un-normalise the gaussian width
    width = width*std_X;
    %un-normalise the k_means
    k_means = k_means*std_X + mean_X;
    %un-normalise the mse
    mse_array = mse_array * std_Y^2;
    %un-normalise the covariance matrix
    beta_cov_array = beta_cov_array.*repmat([std_Y^2,std_Y^2/std_X;std_Y^2/std_X,std_Y^2/std_X^2],1,1,K);
    %un-normalise the regression parameters
    beta_array(1,:) = beta_array(1,:)*std_Y + mean_Y - mean_X*beta_array(2,:)*std_Y/std_X;
    beta_array(2,:) = beta_array(2,:)*std_Y/std_X;

    %plot heatmap of the histogram
    ax = plotHistogramHeatmap(sample_mean,sample_var);
    set(ax, 'Position', [600,800,400,300]);
    xlabel('Sample mean (arb. unit)');
    ylabel('Sample variance {(arb. unit^2)}');
    hold on;

    %for each model
    for i = 1:K
        %set the range of x to plot for a given centroid
        x_min = k_means(i)-0.5E4;
        x_max = k_means(i)+0.5E4;
        x_length = x_max-x_min;
        %vector of means to plot
        x_plot = x_min:(x_length/1000):x_max;
        %convert it to design matrix
        X_plot = [ones(numel(x_plot),1),x_plot'];
        
        %using the values in x_plot, get the estimate of the variance
        y_hat = X_plot*beta_array(:,i);
        %get the standard error
        y_se = sqrt(diag(X_plot*beta_cov_array(:,:,i)*X_plot')+mse_array(i));
        
        %plot the mean and 95% prediction interval
        plot(x_plot,y_hat,'r--');
        plot(x_plot,y_hat+norminv(0.975)*y_se,'r--');
        plot(x_plot,y_hat-norminv(0.975)*y_se,'r--');
        %plot the centroid of the sample mean on the prediction mean
        scatter(k_means(i),[1,k_means(i)]*beta_array(:,i),'r');
    end
    hold off;
    
    %normalise the k_means so it can be used next for loop
    k_means = (k_means - mean_X)/std_X;
end