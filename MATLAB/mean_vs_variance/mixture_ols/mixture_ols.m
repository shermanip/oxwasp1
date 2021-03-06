clc;
clear all;
close all;

%SCRIPT: FIT MIXTURE OF LINEAR REGRESSION
%USER DEFINED VARIABLES
rng(24255965); %set random seed
m = 3; %number of mixtures
n_EM = 50; %number of EM steps
n_samples = 10000; %number of samples of mean-variance pair to be used
folder_location = '/data/tinamou/sip/block_images/orginial'; %location of the data
n_repeat = 10; %number of times to repeat the analysis

%PROGRAM STARTS HERE

%load the mean-variance pair data
[sample_mean,sample_var,area] = load_meanVariance(folder_location);

%take samples of mean-variance pair
index = randperm(area,n_samples);
area = n_samples;
X = sample_mean(index); %assign the samples of means
Y = sample_var(index); %assign the samples of variances

%normalize the data, estimate the mean and std of the mean and variance
mean_Y = mean(Y); %mean of the variances
std_Y = std(Y); %std of the variances
mean_X = mean(sample_mean); %mean of the means
std_X = std(sample_mean); %std of the means

%normalize the data to have 0 mean, std 1
Y = (Y-mean_Y)/std_Y;
X = [ones(area,1),(X-mean_X)/std_X]; %X is a design matrix

%define mixture variables
mixing_array = zeros(m,1); %column vector of mixing ratios
beta_array = zeros(2,m); %column matrix of regression parameters
variance_array = zeros(m,1); %column vector of mse
lnL_array = zeros(n_EM,n_repeat); %column vector of log likelihood at each EM step

%repeat n_repeat times
for j = 1:n_repeat
    %initial E STEP
    %assign random responsbility to each data point
    %responsbility is a matrix of posterior probability of the data being in the ith model
        %dim 1: for each data
        %dim 2: for each model
    responsbility = rand(area,m); %assign it randomly
    %normalize it so that the row sum is one
    responsbility = responsbility./repmat(sum(responsbility,2),1,m);
    
    %do EM algorithm for n_EM times
    for steps = 1:n_EM
        %M STEP, for each model
        for i = 1:m
            %get the weights, aka the responsbility of model i
            w = responsbility(:,i);
            %wX = diag(w)*X
            wX = X;
            for feature = 1:2
                wX(:,feature) = X(:,feature).*w;
            end
            %calculate the weighted least squares
            beta_array(:,i) = (X'*wX)\(X'*(w.*Y));
            %calculate the weighted MSE
            variance_array(i) = sum(w.*((Y-X*beta_array(:,i)).^2)) / sum(w);
        end
        %calculate the mean responsbility
        mixing_array = (sum(responsbility,1)/area)';

        %E STEP, for each model
        for i = 1:m
            %calculate the posterior probability (to a constant) for each data belonging to model i
            responsbility(:,i) = linearReg_likelihood(Y,X,beta_array(:,i),variance_array(i))*mixing_array(i);
        end
        %normalize the responsbility so that the sum of the rows = 1
        responsbility = responsbility./repmat(sum(responsbility,2),1,m);

        %save the log likelihood
        lnL_array(steps,j) = mixture_logLikelihood(Y,X,beta_array,variance_array,mixing_array);
    end
end

%un-normalize the regression parameters
beta_array(1,:) = beta_array(1,:)*std_Y + mean_Y - mean_X*beta_array(2,:)*std_Y/std_X;
beta_array(2,:) = beta_array(2,:)*std_Y/std_X;
prediction_error = sqrt(variance_array.*std_Y^2);

%plot heatmap of the histogram
plotHistogramHeatmap(sample_mean,sample_var);
hold on;

%get a range of means to be used to plot the mixture of linear regression
x_min = min(sample_mean);
x_max = max(sample_mean);
x_length = x_max-x_min;
x_plot = x_min:(x_length/1000):x_max;
%for each model
for i = 1:m
    %using the values in x_plot, plot the model
    y_hat = [ones(numel(x_plot),1),x_plot']*beta_array(:,i);
    plot(x_plot,y_hat,'r-');
    plot(x_plot,y_hat+norminv(0.975)*prediction_error(i),'r:');
    plot(x_plot,y_hat-norminv(0.975)*prediction_error(i),'r:');
end
hold off;

%plot the log likelihood at each EM step
figure;
plot(lnL_array,'b');
xlabel('Number of EM steps');
xlim([1,n_EM]);
ylabel('Log likelihood (nat)');