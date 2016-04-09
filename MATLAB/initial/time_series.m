clc;
clear all;
close all;

%TO BE RUN LOCALLY

%array of images
length = 1996; %length of the data
area = length^2; %area of the data
n = 100; %sample size

%array of mean and standard deviation of each image
greyValue_mean_std = zeros(n,2);

%for each image, save the mean and std
for i = 1:n
    stack = double(imread(strcat('/data/tinamou/sip/block_images/orginial/block_',num2str(i),'.tif')));
    stack = reshape(stack,[],1);
    greyValue_mean_std(i,1) = mean(stack);
    greyValue_mean_std(i,2) = std(stack);
end

%error bar mean and standard error
figure;
errorbar(greyValue_mean_std(:,1),greyValue_mean_std(:,2)/sqrt(area));
xlabel('Observation number');
ylabel('Mean grey value');
xlim([1,n]);

%plot autocorrelation
ax = figure;
autocorr(greyValue_mean_std(:,1));
title('');
set(ax, 'Position', [600,800,400,300]);

%plot partial autocorrelation
ax1 = figure;
parcorr(greyValue_mean_std(:,1));
title('');
set(ax1, 'Position', [600,800,400,300]);