%SCRIPT: WEIGHTED LEAST SQUARES
    %do weighted least squares regression on the mean and variance data,
    %the weights are prop. to 1/variance^2. The BIC was calculated for each
    %polynomial order and each bootstrap samples

%set the number of bootstrap samples to do
n_repeat = 40;
%maximum polynomial order to investigate
p_max = 6;

%set the random seed
rng(27689355);
%set seeds for each repeat of the experiment
seeds = randi([0,intmax],n_repeat-1,1);

%load the data
folder_location = '/data/tinamou/sip/block_images/orginial'; %location of the data
[sample_mean,sample_var,area] = load_meanVariance(folder_location);

%define array of BIC for each polynomial and each bootstrap
    %column: for each polynomial order
    %row: for each bootstrap
BIC_array = zeros(n_repeat,p_max);

%do weighted least squares fit on the orginial data and get the BIC and plot the regression
BIC = samplingDist_weighted_BIC(sample_mean,sample_var,1,p_max);
%save the BIC in the array
BIC_array(1,:) = BIC;

%for each bootstrap
parfor i = 2:n_repeat
    %set the random seed
    rng(seeds(i-1));
    %get the index of the boostrap samples
    bootstrap_index = randi([1,area],area,1);
    %do weighted least squares fit on the boostrap sample and get the BIC
    BIC = samplingDist_weighted_BIC(sample_mean(bootstrap_index),sample_var(bootstrap_index),0,p_max);
    %save the BIC in the BIC array
    BIC_array(i,:) = BIC;
end

%box plot the BIC for each polynomial
figure;
boxplot(BIC_array);
xlabel('Polynomial order');
ylabel('BIC (nat)');