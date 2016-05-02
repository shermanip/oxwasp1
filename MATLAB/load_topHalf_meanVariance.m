%FUNCTION: LOAD THE MEAN-VARIANCE DATA
%PARAMETERS:
    %folder_location: string containing the location of the file
function [sample_mean,sample_var,area] = load_topHalf_meanVariance(folder_location)
    %load the data
        %stack is an array of matrices containing the grey values
        %area is the number of pixels in each image
    [stack,~,area,~] = load_topHalf_stack(folder_location);
    
    %work out the sample mean
    sample_mean = reshape(mean(stack,3),[],1);
    %work out the sample std
    sample_var = reshape(var(stack,0,3),[],1);
end

