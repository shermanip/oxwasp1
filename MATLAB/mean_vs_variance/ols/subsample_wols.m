clc;
clear all;
close all;

%SCRIPT: weighted linear regression on subsample of SAMPLE, BACKGROUND and FOAM
%plots the lienar regressions

rng(93016144);
n_repeat = 40;
seeds = randi([0,93016144],3,40);

n_bins = [1000,100,100];

subsample_index = [
    236,1506,328,1398;
    32, 964,22,190;
    1942,32,12,1972
];

folder_location = '/data/tinamou/sip/block_images/orginial'; %location of the data

%array of images
[stack,length,area,n] = load_stack(folder_location);

%for SAMPLE, BACKGROUND, FOAM subsamples, plot histogram heatmap and plot linear regression
%fit linear regression of subsample of the SAMPLE
    %beta_1 = regression vector parameter
    %n_1 = number of sub-sample
    %model_1 = linear model object
    %range_1 = [min,max] of the mean

for i = 1:3
    sub_sample = stack(subsample_index(i,1):(subsample_index(i,1)+subsample_index(i,2)),subsample_index(i,3):(subsample_index(i,3)+subsample_index(i,4)),:);
    sample_mean = reshape(mean(sub_sample,3),[],1);
    sample_var = reshape(var(sub_sample,0,3),[],1);
    samplingDist_weighted_BIC(sample_mean,sample_var,1,2,n_bins(i));
end

format longE
disp('BIC');
for i = 1:3
    disp(cell2mat(strcat({'For material '},num2str(i))));
    sub_sample = stack(subsample_index(i,1):(subsample_index(i,1)+subsample_index(i,2)),subsample_index(i,3):(subsample_index(i,3)+subsample_index(i,4)),:);
    sample_mean = reshape(mean(sub_sample,3),[],1);
    sample_var = reshape(var(sub_sample,0,3),[],1);
    BIC_array = zeros(n_repeat,2);
    parfor j = 1:n_repeat
        rng(seeds(i,j));
        bootstrap_index = randi([1,numel(sample_mean)],numel(sample_mean),1);
        BIC = samplingDist_weighted_BIC(sample_mean(bootstrap_index),sample_var(bootstrap_index),1,2,n_bins(i));
        BIC_array(j,:) = BIC;
    end
    disp('Order 1     Order 2');
    disp('Mean');
    disp(mean(BIC_array));
    disp('Std');
    disp(std(BIC_array));
end
