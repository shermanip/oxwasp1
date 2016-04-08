%FUNCTION: FACTOR ANALYSIS BOOTSTRAP BIC - RUN ON GREYPARTRIDGE
    %STEP 2: get a bootstrap sample and fit FA onto it for a range of k
    %append the results to the binary file
function factorAnalysis_bootstrapBIC_learn()

    %load the variables and past results
        %k_max: maximum number of factors to investigate
        %n_EM: number of EM steps
        %n_run: number of initial values to run
        %lnL_array: n_bootstrap x k_max matrix of log likelihoods
        %n_bootstrap: number of bootstraps done
    load('/data/greypartridge/oxwasp/oxwasp15/sip/fa_bootstrap.mat');
    %make the variables required for each worker visible
    n_run = n_run;
    n_EM = n_EM;
    
    %vector of log likelihoods for this single bootstrap
    lnL_array_bootstrap = -inf*ones(1,k_max);

    %dataset properties
    length = 100;
    area = length^2;
    n = 100;
    %array of images
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

    %get a bootstrap sample
    X = X(randi(n,n,1),:);

    %center the data at 0
    mean_x = mean(X,1);
    X = X - repmat(mean_x,n,1);
    %estimate the covariance matrix
    cov_x = cov(X);

    %for each k
    parfor k = 1:k_max
        %get the log likelihood to be -infinity
        lnL_k = -inf;
        %for n_run times
        for j = 1:n_run
            %declare the loading matrix and the noise randomly
            [loading,noise_vector] = factorAnalysis_EM_initalize(area,k);
            %for n_EM times
            for i = 1:n_EM
                %E STEP
                [Y,Y_cov] = factorAnalysis_EStep(loading,noise_vector,X,k);
                %M STEP
                [loading,noise_vector] = factorAnalysis_MStep(X,Y,Y_cov,n,area);  
            end
            %get the log liklihood
            lnL = factorAnalysis_lnL(loading,noise_vector,cov_x,n,area);
            %if the log likelihood is bigger than lnL_k, update it
            if lnL>lnL_k
                lnL_k = lnL;
            end
        end
        %save the log liklihood
        lnL_array_bootstrap(k) = lnL_k;
    end

    %append results of the log likelihood
    lnL_array = [lnL_array;lnL_array_bootstrap];
    %increment the number of bootstraps
    n_bootstrap = n_bootstrap+1;
    %display the number of bootstraps done
    disp(strcat(num2str(n_bootstrap),' bootstraps done'));

    %save the variables to the binary file
    save('/data/greypartridge/oxwasp/oxwasp15/sip/fa_bootstrap.mat','k_max','n_EM','n_run','lnL_array','n_bootstrap');
end