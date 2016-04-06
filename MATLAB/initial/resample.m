clc;
clear all;
close all;

%TO BE RUN ON A SERVER, NUMBER OF CORES TO BE PRE-DEFINED

%array of images
length = 1996; %length of the data
area = length^2; %area of the data
n = 100; %sample size
stack = zeros(length,length,n);

%for each image, save the pixel values
for i = 1:n
    stack(:,:,i) = imread(strcat('/data/tinamou/sip/block_images/orginial/block_',num2str(i),'.tif'));
end

%set size of the resize image
shrink_size = 100;
%define stack of resize images
resample_stack = zeros(shrink_size,shrink_size,n);
%select which pixel to sample from
index_grid = round(length*(1:shrink_size)/(shrink_size+1));
%for each image, resample
for i = 1:n
    resample_stack(:,:,i) = stack(index_grid,index_grid,i);
end

%save it
save('/data/tinamou/sip/block_images/resample.mat','resample_stack');