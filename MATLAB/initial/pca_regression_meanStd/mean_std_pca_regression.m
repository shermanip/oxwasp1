clc;
clear all;
close all;

%TO BE RUN LOCALLY (can be extended to server by increasing n_bootstrap)

%array of images
length = 1996; %length of the images
area = length^2; %area of the images
n = 100; %number of images
stack = zeros(length,length,n);

%for each image, save the pixel values
for i = 1:n
    stack(:,:,i) = imread(strcat('/data/tinamou/sip/block_images/orginial/block_',num2str(i),'.tif'));
end

%work out the sample mean
sample_mean = reshape(mean(stack,3),[],1);
%work out the sample std
sample_std = reshape(std(stack,0,3),[],1);

%number of times to shake the data
n_bootstrap = 100;
%define arrays to store vectors
com_array = zeros(n_bootstrap,2); %center of mass
u_array = zeros(n_bootstrap,2); %PC1
y_error_array = zeros(n_bootstrap,1); %error

%for n_boostrap times
parfor i = 1:n_bootstrap
    
    %shake the data according to the sampling distribution
    sample_mean_shake = sample_mean + trnd(n-1,area,1).*sample_std/sqrt(n);
    sample_std_shake = sqrt(((n-1)*(sample_std).^2)./chi2rnd(n-1,area,1));
    
    %bootstrapping
    %sample_mean_shake = sample_mean(randi(area,area,1));
    %sample_std_shake = sample_std(randi(area,area,1));
    
    %do PCA_regression
        %r_com is the center of mass
        %u_1 is the PC1 vector
        %y_error is the mean squared error in the PC2 direction, projected in the y direction
    [r_com,u_1,y_error] = PCA_regression(sample_mean_shake,sample_std_shake);
    %save u_1 in u_array
    u_array(i,:) = u_1';
    %save r_com in com_array
    com_array(i,:) = r_com';
    %save y_error in  y_error_array
    y_error_array(i) = y_error;
end

%get the mean center of mass
r_com = (mean(com_array,1))';
x_bar = r_com(1);
y_bar = r_com(2);
%get the std of the center of mass
r_com_std = (std(com_array,0,1))';

%get the mean y_error
y_error = mean(y_error_array);

%work out the mean and standard deviation of the gradient
m_array = u_array(:,2)./u_array(:,1); %work out the graident
m = mean(m_array); %estimate the mean
m_std = std(m_array); %estimate the standard deviation

%bin the data (1000 x 1000 bins)
    %N is the count in each bin
    %c are the centers of each bin
[N,c] = hist3([sample_std,sample_mean],[1000,1000]);
%get the x values of the plots
x_plot = cell2mat(c(2));
%get the y values of the plots
y_plot = m*(x_plot-x_bar)+y_bar;

%plot the histogram as a heat map
figure;
%normalize N so that the colormap is the frequency density
imagesc(cell2mat(c(2)),cell2mat(c(1)),N/( (c{2}(2)-c{2}(1))*(c{1}(2)-c{1}(1)) ) );
axis xy; %switch the y axis
ylim([150,700]); %set the y limit
colorbar; %display the colour bar
hold on;
%plot the PCA regression
plot(x_plot,y_plot,'r--'); %mean
plot(x_plot,y_plot+y_error,'r--'); %mean +/- std
plot(x_plot,y_plot-y_error,'r--'); %mean +/- std
xlabel('Sample grey value mean (AU)');
ylabel('Sample grey value std (AU)');
hold off;

%set the format to be in long scientific notation
format longE;

%print the gradient mean and standard deviation
disp('Gradient mean');
disp(m);
disp('Gradient std');
disp(m_std);

%print the COM mean and standard deviation
disp('COM mean (mean_greyValue, std_greyValue)');
disp(r_com);
disp('COM std (mean_greyValue, std_greyValue)');
disp(r_com_std);