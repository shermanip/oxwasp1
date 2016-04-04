clc;
clear all;
close all;

%array of images
stack = zeros(1996,1996,100);

%for each image, save the pixel values
n = 100;
for i = 1:n
    stack(:,:,i) = imread(strcat('/data/tinamou/sip/block_images/block_',num2str(i),'.tif'));
end

%work out the sample mean
sample_mean = reshape(mean(stack,3),[],1);
%work out the sample std
sample_std = reshape(std(stack,0,3),[],1);

%get the number of times to shake the data
n_bootstrap = 100;
%define arrays to store vectors
com_array = zeros(n_bootstrap,2); %center of mass
u_array = zeros(n_bootstrap,2); %PC1
y_error_array = zeros(n_bootstrap,1); %error

for i = 1:n_bootstrap
    
    %shake the data
    sample_mean_shake = sample_mean + trnd(99,1996^2,1).*sample_std/sqrt(100);
    sample_std_shake = sqrt((99*(sample_std).^2)./chi2rnd(99,1996^2,1));
    
%     sample_mean_shake = sample_mean(randi(1996^2,1996^2,1));
%     sample_std_shake = sample_std(randi(1996^2,1996^2,1));
    
    %do PCA_regression
    [r_com,u_1,y_error] = PCA_regression(sample_mean_shake,sample_std_shake);
    %get the gradient
    u_array(i,:) = u_1';
    %get the center of mass
    com_array(i,:) = r_com';
    y_error_array(i) = y_error;
end

%get the mean center of mass
r_com = (mean(com_array))';
x_bar = r_com(1);
y_bar = r_com(2);
%get the std of the center of mass
r_com_std = (std(com_array))';

%get the mean y_error
y_error = mean(y_error_array);

%work out the mean and stdgradient
m_array = u_array(:,2)./u_array(:,1);
m = mean(m_array);
m_std = std(m_array);

%bin the data
[N,c] = hist3([sample_std,sample_mean],[1000,1000]);
%get the x values of the plots
x_plot = cell2mat(c(2));
%get the y values of the plots
y_plot = m*(x_plot-x_bar)+y_bar;

%plot the histogram as a heat map
figure;
imagesc(cell2mat(c(2)),cell2mat(c(1)),N/( (c{2}(2)-c{2}(1))*(c{1}(2)-c{1}(1)) ) );
axis xy; %switch the y axis
ylim([150,700]); %set the y limit
colorbar; %display the colour bar
hold on;
%plot the PCA regression
plot(x_plot,y_plot,'r--');
plot(x_plot,y_plot+y_error,'r--');
plot(x_plot,y_plot-y_error,'r--');
xlabel('Sample grey value mean (AU)');
ylabel('Sample grey value std (AU)');
hold off;

format longE;
disp('Gradient mean');
disp(m);
disp('Gradient std');
disp(m_std);

disp('COM mean');
disp(r_com);
disp('COM std');
disp(r_com_std);