%PROCEDURE: FACTOR ANALYSIS BIC - TO BE RUN LOCALLY
    %STEP 2: fit FA onto the data for a range of k, update:
        %lnL_array - array of loglikehoods for each k
        %BIC_optimal - the optimial BIC found
        %loading_optimal - the optimial loading matrix
        %noise_optimal - the optimial vector of instrinic noise
        %n_totalRun - number of initial values run for EM
    %and save it to a binary file
function factorAnalysis_BIC_learn(n_run)
    %load initalize variables and past results
    load('/data/tinamou/sip/results/fa_bic.mat'); 
        %k_max
        %n_EM
        %lnL_array
        %BIC_optimial
        %n_totalRun
        %loading_optimal
        %noise_optimal

    %define properties of the dataset
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

    %for each k
    for k = 1:k_max
        %get the maximum log likelihood found from the array
        lnL_k = lnL_array(k);
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
            %calculate the BIC
            BIC = -2*lnL+log(n)*area*(k+1);
            %if BIC is smaller than BIC_optimal, update the BIC_optimal, loading and noise optimal
            if BIC<BIC_optimal
                BIC_optimal = BIC;
                loading_optimal = loading;
                noise_optimal = noise_vector;
            end
        end
        %save the log liklihood in the array
        lnL_array(k) = lnL_k;
        %display the progress - number of factors done
        disp(strcat(num2str(k),' factors done'));
    end

    %update the number of runs
    n_totalRun = n_totalRun + n_run;
    %display the total number of runs
    disp(strcat(num2str(n_totalRun),' runs done so far'));
    
    %save the variables
    save('/data/tinamou/sip/results/fa_bic.mat','k_max','n_EM','lnL_array','BIC_optimal','n_totalRun','loading_optimal','noise_optimal');
end