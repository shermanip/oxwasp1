clc;
clear all;
close all;

%STEP 3: plot results

%load the results
load('fa_bic.mat');
length = 100;
area = length^2;
n = 100;

%work out the BIC and plot it
BIC = -2*lnL_array'+log(n)*area*((1:k_max)+1);
figure;
plot(BIC);

%plot the individual factor noises
loading_max = max(max(loading_optimal.^2)); %find the maximum element squared
loading_min = min(min(loading_optimal.^2)); %find the minimum element squared
figure;
%for each factor
for i = 1:8
    subplot(2,4,i);
    %heat map plot the noise
    imagesc(reshape((loading_optimal(:,i)).^2,length,length),[loading_min,loading_max]);
    set(gca,'xtick',[]); %turn the ticks off
    set(gca,'ytick',[]); %turn the ticks off
end
%place the colorbar
hp4 = get(subplot(2,4,4),'Position');
colorbar('Position', [0.05  0.04  0.05  hp4(2)+hp4(3)*2.1]);

%plot the instrinic noise
figure;
imagesc(reshape(noise_optimal,length,length));
colorbar;
set(gca,'xtick',[]); %turn the ticks off
set(gca,'ytick',[]); %turn the ticks off