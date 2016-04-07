clc;
clear all;
close all;

%TO BE RUN LOCALLY

%USER DEFINE VARIABLES
%number of factors
k_max = 4;
%number of EM steps
n_EM = 50;

%array of images
length = 100;
area = length^2;
n = 100;
stack = zeros(length,length,n);

%LOAD VARIABLES

%for each image, save the pixel values
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

%declare array of noise images
factor_noise_array = zeros(length,length,4);
intrinsic_noise_array = zeros(length,length,4);

%EM ALGORITHM

%for factors 1,2,3,...,k_max
parfor k = 1:k_max

    %declare the loading matrix and the noise
    loading = normrnd(0,1E3,area,k);
    noise_vector = gamrnd(1,1E3,area,1); %vector of noise

    %repeat n_EM times
    for j = 1:n_EM
        %E STEP
        [Y,Y_cov] = factorAnalysis_EStep(loading,noise_vector,X,k);

        %M STEP
        [loading,noise_vector] = factorAnalysis_MStep(X,Y,Y_cov,n,area);   
    end
    
    %save the factor and intrinsic noise
    factor_noise_array(:,:,k) = reshape(diag(loading*loading'),length,length);
    intrinsic_noise_array(:,:,k) = reshape(noise_vector,length,length);
end

%PLOTTING factor and intrinsic noise
figure;
for j = 1:k_max %for each factor
    
    %plot factor noise
    subplot(k_max,2,2*j-1) %define subplot
    %heatmap
    imagesc(factor_noise_array(:,:,j));
    colorbar;
    %put the type of noise at the top
    if (j==1)
        title('Factor noise');
    end
    ylabel(strcat('k=',num2str(j)));
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    
    %plot intrinsic noise
    subplot(k_max,2,2*j);
    %heatmap
    imagesc(intrinsic_noise_array(:,:,j));
    colorbar;
    %put the type of noise at the top
    if (j==1)
        title('Intrinsic noise');
    end
    ylabel(strcat('k=',num2str(j)));
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
end