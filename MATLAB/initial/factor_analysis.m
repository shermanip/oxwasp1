clc;
clear all;
close all;

%number of factors
k = 3;
%number of EM steps
n_EM = 3;

%array of images
length = 100;
area = length^2;
stack = zeros(length,length,100);

%LOAD VARIABLES

%for each image, save the pixel values
n = 100;
for i = 1:n
    if i<=10
        stack(:,:,i) = imread(strcat('/data/tinamou/sip/block_images/stack_100/compress_block_images_stack00',num2str(i-1),'.tif'));
    else
        stack(:,:,i) = imread(strcat('/data/tinamou/sip/block_images/stack_100/compress_block_images_stack0',num2str(i-1),'.tif'));
    end
end

%define design matrix
X = zeros(n,area);
%reshape the image to a row vector
for i = 1:n
    X(i,:) = reshape(stack(:,:,i),1,area);
end

%center the data at 0
mean_x = mean(X,1);
X = X - repmat(mean_x,n,1);
%estimate the covariance matrix
cov_x = cov(X);

%declare the loading matrix and the noise
loading = normrnd(0,30,area,k);
noise_vector = gamrnd(1,1E3,area,1); %vector of noise

%declare vector of lnL
lnL = zeros(n_EM,1);

%start online figure of the factor noise and intrinsic noise
ax = figure;
colormap gray; %colourmap
subplot(2,1,1) %plot factor noise
imagesc(reshape(diag(loading*loading'),length,length));
subplot(2,1,2); %plot intrinsic noise
imagesc(reshape(noise_vector,length,length))
ax.NextPlot = 'replaceChildren';

%INITIAL E-STEP
[Y,Y_cov] = factorAnalysis_EStep(loading,noise_vector,X,k);

%repeat n_EM times
for j = 1:n_EM

    %M STEP
    [loading,noise_vector] = factorAnalysis_MStep(X,Y,Y_cov,n,area);

    %plot the factor and intrinsic noise
    subplot(2,1,1) %factor
    imagesc(reshape(diag(loading*loading'),length,length));
    subplot(2,1,2); %intrinsic
    imagesc(reshape(noise_vector,length,length))
    drawnow;
    
    %E STEP
    [Y,Y_cov] = factorAnalysis_EStep(loading,noise_vector,X,k);

    %get log likelihood
    lnL(j) = factorAnalysis_lnL(loading,noise_vector,cov_x,n,area);
end

figure;
plot(1:n_EM,lnL);
xlim([1,n_EM]);

% x_bar = reshape(mean(stack,3),1,area);
% x_fa = (loading'*noise_inv*loading+eye(k))\(loading'*noise_inv)*x_bar';
% x_fa = (loading*x_fa)';
% figure;
% colormap gray;
% imagesc(reshape(x_fa+mean_x,length,length));