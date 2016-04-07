clc;
clear all;
close all;

%TO BE RUN ON A SERVER: n_run cores

%USER DEFINED
%number of factors
k_max = 4;
%number of EM steps
n_EM = 8;
%number of initial values
n_run = 12;

%array of images
length = 100;
area = length^2;
n = 100;
stack = zeros(length,length,n);

%LOAD VARIABLES

%for each image, save the pixel values
for i = 1:n
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

%EM ALGORITHM

%for factors 1,2,3,...,k_max
for k = 1:k_max
    %n_run times in parallel
    parfor i = 1:n_run
        %declare the loading matrix and the noise
        loading = normrnd(0,1E3,area,k);
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

%PLOT

%plot the log likelihoods
plot_array = []; %array of plots
figure;
colour = hsv; %select colour of the lines from hsv
%for each k
for j = 1:k_max
    %plot the lnL and save the plot in plot_array
    plot_array(j,:) = plot(lnL_array(:,:,j),'Color',colour(floor(j*64/k_max),:));
    hold on;
end
%set limits and labels
xlim([1,n_EM]);
xlabel('Number of EM steps');
ylabel('lnL (nat)');
hold off;
%plot the legend using plot_array
legend(plot_array(:,1),'k=1','k=2','k=3','k=4');

%save the variables
save('/data/greypartridge/not-backed-up/oxwasp/oxwasp15/sip/EMlnL.mat','lnL_array');