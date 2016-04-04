clc;
clear all;
close all;

%number of factors
k_max = 4;
%number of EM steps
n_EM = 8;
n_run = 12;

%array of images
length = 100;
area = length^2;
stack = zeros(length,length,100);

%LOAD VARIABLES

%for each image, save the pixel values
n = 100;
for i = 1:n
%     if i<=10
%         stack(:,:,i) = imread(strcat('/data/tinamou/sip/block_images/stack_100/compress_block_images_stack00',num2str(i-1),'.tif'));
%     else
%         stack(:,:,i) = imread(strcat('/data/tinamou/sip/block_images/stack_100/compress_block_images_stack0',num2str(i-1),'.tif'));
%     end
    if i<=10
        stack(:,:,i) = imread(strcat('/data/greypartridge/not-backed-up/oxwasp/oxwasp15/sip/stack_100/compress_block_images_stack00',num2str(i-1),'.tif'));
    else
        stack(:,:,i) = imread(strcat('/data/greypartridge/not-backed-up/oxwasp/oxwasp15/sip/stack_100/compress_block_images_stack0',num2str(i-1),'.tif'));
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

%declare array of log likelihoods
lnL_array = zeros(n_EM,n_run,k_max);

%for factors 1,2,3,...,k_max
for k = 1:k_max
    %n_run times in parallel
    parfor i = 1:n_run
        %declare the loading matrix and the noise
        loading = normrnd(0,30,area,k);
        noise_vector = gamrnd(1,1E3,area,1); %vector of noise

        %repeat n_EM times
        for j = 1:n_EM
            %E STEP
            [Y,Y_cov] = factorAnalysis_EStep(loading,noise_vector,X,k);

            %M STEP
            [loading,noise_vector] = factorAnalysis_MStep(X,Y,Y_cov,n,area);
            
            %save the log likelihood
            lnL_array(j,i,k) = factorAnalysis_lnL(loading,noise_vector,cov_x,n,area);
        end
    end
end

%plot it
plot_array = [];
figure;
colour = hsv;
for j = 1:k_max
    plot_array(j,:) = plot(lnL_array(:,:,j),'Color',colour(floor(j*64/k_max),:));
    hold on;
end
xlim([1,n_EM]);
xlabel('Number of EM steps');
ylabel('lnL (nat)');
hold off;
legend(plot_array(:,1),'k=1','k=2','k=3','k=4');

%save the variables
save('EMlnL.mat','lnL_array');