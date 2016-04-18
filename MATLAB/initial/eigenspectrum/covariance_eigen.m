clc;
clear all;
close all;

%TO BE RUN ON A SERVER, SET NUMBER OF CORES = 16

%user define
n_bootstrap = 100;
parallel_profile = gcp;
n_core = parallel_profile.NumWorkers;

%array of images
length = 100; %length of image
area = length^2; %area of image
n = 100; %number of images
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

%DATA ANALYSIS

%estimate the correlation coefficient
R = corrcoef(X);

%declare array of eigenvalues for each bootstrap and core
    %dim1: for each bootstrap
    %dim2: for each PC
    %dim3: for each core
eigenvalue_array_core = zeros(n_bootstrap,6,n_core);
%declare array of PC1 vectors for each bootstrap and core
    %dim1: for each pixel
    %dim2: for each bootstrap
    %dim3: for each core
PC1_array_core = zeros(area,n_bootstrap,n_core);

%declare array of eigenvalues and PC1 for each bootstrap + 1 (combining the cores together)
eigenvalue_array = zeros(n_bootstrap*n_core+1,6); %dim1: boostrap. dim2: eigen
PC1_array = zeros(area,n_bootstrap*n_core+1); %dim1: pixel, dim2: boostrap

%estimate the covariance matrix
Sigma = cov(X,1);
%without bootstrap
%get the few biggest eigenvector and eigenvalues and save it
[eigenvector,eigenvalue] = eigs(Sigma,6);
S = eigenvector;
%save the non-bootstrap eigenvalue and PC1 at the end of eignevalue_array and PC1_array
eigenvalue_array(end,:) = (diag(eigenvalue))';
PC1_array(:,end) = S(:,1);

%for each core
parfor i = 1:n_core
    %n_bootstrap times
    for j = 1:n_bootstrap
        %get a bootstrap sample of X and estimate the covariance
        Sigma = cov(X(randi(n,n,1),:),1);
        %get the few biggest eigenvalues and save it
        [eigenvector,eigenvalue] = eigs(Sigma,6);
        eigenvalue_array_core(j,:,i) = (diag(eigenvalue))';
        PC1_array_core(:,j,i) = eigenvector(:,1);
    end
end

%combine the matrices from each core
for i = 1:n_core
    eigenvalue_array( ((i-1)*n_bootstrap+1) : i*n_bootstrap, :) = eigenvalue_array_core(:,:,i);
end
for i = 1:n_core
    PC1_array(:, ((i-1)*n_bootstrap+1) : i*n_bootstrap) = PC1_array_core(:,:,i);
end


%PLOT%%%

%plot the correlation coefficient as a heat map
figure;
imagesc(R,[-1,1]);
colorbar;

%box plot the bootstrap eigenvalues
figure;
boxplot(eigenvalue_array);
xlabel('Eigenvalue number');
ylabel('Eigenvalue (AU)');

%for each eigenvector, plot the non-bootstrap eigenimage
figure('Position', [100, 100, 800, 400]);
for i = 1:6
    subplot(2,3,i);
    imagesc(((reshape((S(:,i))',length,length))).^2);
    colorbar;
    set(gca,'xtick',[]); %remove x tick
    set(gca,'ytick',[]); %remove y tick
    title(strcat('PC',num2str(i))); %put title
end

%express the mean image using PC basis %LOAD STACK
x_mean = reshape(mean(stack,3),1,area);
figure;
for i = 1:6
    %get the transformation matrix
    S_i = (S(:,1:i));
    %transform the data
    x_bar = S_i*(S_i'*x_mean');
    %plot it
    subplot(2,3,i);
    colormap gray;
    imagesc((reshape(x_bar',length,length)));
    set(gca,'xtick',[]); %remove x tick
    set(gca,'ytick',[]); %remove y tick
    title(strcat(num2str(i),' PC basis'));
end

%get a design matrix of the coefficient of the PC squared
%LOAD X
pca_coefficient = (X*S).^2;
%plot it as a heat map
figure;
imagesc(pca_coefficient);
colorbar;
xlabel('PCA number');
ylabel('Data number');

%plot the time series of the 1st PC coefficient bootstrapping
figure;
plot(pca_coefficient(:,1));
ylabel('PCA1 coefficient^{2} (AU^{2})');
xlabel('Observation number');
xlim([1,n]);

%using the bootstrap samples, get the 1st PC coefficient
%LOAD X
pca1_timeseries = (X*PC1_array).^2;
%plot it with std
figure;
errorbar(mean(pca1_timeseries,2),std(pca1_timeseries,0,2));
ylabel('PCA1 coefficient^{2} (AU^{2})');
xlabel('Observation number');
xlim([1,n]);
%plot it with standard error
figure;
errorbar(mean(pca1_timeseries,2),std(pca1_timeseries,0,2)/sqrt(n_core*n_bootstrap));
ylabel('PCA1 coefficient^{2} (AU^{2})');
xlabel('Observation number');
xlim([1,n]);

%save the variables
save('/data/greypartridge/not-backed-up/oxwasp/oxwasp15/sip/eigenspectrum.mat','R','eigenvalue_array','S','PC1_array','n_core','n_bootstrap');