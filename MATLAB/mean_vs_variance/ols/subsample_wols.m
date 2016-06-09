clc;
clear all;
close all;

%SCRIPT: weighted linear regression on subsample of SAMPLE, BACKGROUND and FOAM
%plots the lienar regressions

%set random seed and for each core
rng(93016144);
n_repeat = 40;
seeds = randi([0,93016144],3,40);

%set histogram bins for materials
n_bins = [1000,100,100];

%set the position of the resampled
    %rows: for each material
    %column: [x position of top left, width, y position of top left, height]
subsample_index = [
    236,1506,328,1398;
    32, 964,22,190;
    1942,32,12,1972
];

folder_location = '/data/tinamou/sip/block_images/orginial'; %location of the data

%array of images
[stack,length,area,n] = load_stack(folder_location);

%for SAMPLE, BACKGROUND, FOAM subsamples, plot histogram heatmap and plot linear regression
for i = 1:3
    %get the sample of pixels
    sub_sample = stack(subsample_index(i,1):(subsample_index(i,1)+subsample_index(i,2)),subsample_index(i,3):(subsample_index(i,3)+subsample_index(i,4)),:);
    %estimate the mean and reshape it into a vector
    sample_mean = reshape(mean(sub_sample,3),[],1);
    %estimate the variance and reshape it into a vector
    sample_var = reshape(var(sub_sample,0,3),[],1);
    %plot the weighted regression
    samplingDist_weighted_BIC(sample_mean,sample_var,1,2,n_bins(i));
end

format longE
disp('BIC');
%for each mateiral
for i = 1:3
    
    %display the material number
    disp(cell2mat(strcat({'For material '},num2str(i))));
    %get the resampled
    sub_sample = stack(subsample_index(i,1):(subsample_index(i,1)+subsample_index(i,2)),subsample_index(i,3):(subsample_index(i,3)+subsample_index(i,4)),:);
    %estimate the mean and reshape it into a vector
    sample_mean = reshape(mean(sub_sample,3),[],1);
    %estimate the variance and reshape it into a vector
    sample_var = reshape(var(sub_sample,0,3),[],1);
    %declare array of BIC for 1st and 2nd order for each n_repeat
    BIC_array = zeros(n_repeat,2);
    
    %repeat n_repeat times
    parfor j = 1:n_repeat
        %get the core random seed
        rng(seeds(i,j));
        %bootstrap the sample and get the BIC for 1st and 2nd order
        bootstrap_index = randi([1,numel(sample_mean)],numel(sample_mean),1);
        BIC = samplingDist_weighted_BIC(sample_mean(bootstrap_index),sample_var(bootstrap_index),1,2,n_bins(i));
        %save the BIC
        BIC_array(j,:) = BIC;
    end
    
    %print the mean and std BIC for both orders
    disp('Order 1     Order 2');
    disp('Mean');
    disp(mean(BIC_array));
    disp('Std');
    disp(std(BIC_array));
end
