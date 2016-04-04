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
    stack(:,:,i) = imread(strcat('/data/greypartridge/not-backed-up/oxwasp/oxwasp15/sip/block_images/block_',num2str(i),'.tif'));
end

%declare of p values for each pixel in the image
p_value_array = zeros(length,length);

%for each column
parfor j = 1:length
    %for each row
    for i = 1:length
        %get the (i,j) pixel and shape it into a column vector
        x = reshape(stack(i,j,:),[],1);
        %fit Normal distribution
        distribution_fit = fitdist(x,'Normal');
        %get p value for chi squared goodness of fit
        [~,p] = chi2gof(x,'CDF',distribution_fit);
        %save the value
        p_value_array(i,j) = p;
    end
end

%save it
save('/data/greypartridge/not-backed-up/oxwasp/oxwasp15/sip/p_values.mat','p_value_array');

%plot the p values
figure;
imagesc(p_value_array,[0,1]); %set scale from 0 to 1
colorbar; %add color bar

%plot the log 10 scale
log_p = log10(p_value_array); %log10
figure;
%plot heat map of log 10 scale
imagesc(log_p);
colorbar; %add color bar
hold on; %in addition, find [row,col] significant pixels at 10% level
[row_array, column_array] = find(log_p < (-log10(area)-1));
%plot circles on significant pixels (x,y)
scatter(column_array,row_array);
hold off;

%print the quantile of the p values
p_quantile = quantile(reshape(p_value_array,[],1),[0.05,0.5,0.95]);
disp('90% two tailed quantile');
disp(strcat('medium:',num2str(p_quantile(2))));
disp(strcat('upper error:',num2str(p_quantile(3)-p_quantile(2))));
disp(strcat('lower error:',num2str(p_quantile(2)-p_quantile(1))));