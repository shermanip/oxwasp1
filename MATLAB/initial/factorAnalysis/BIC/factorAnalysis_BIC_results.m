%PROCEDURE: FACTOR ANALYSIS BIC - TO BE RUN LOCALLY
    %STEP 3: plot results from the saved binary file
        %individual factor noise
        %instrinic noise
function factorAnalysis_BIC_results()
    %load the results
        %k_max
        %n_EM
        %lnL_array
        %BIC_optimial
        %n_totalRun
        %loading_optimal
        %noise_optimal
    load('/data/tinamou/sip/results/fa_bic.mat');
    
    %display the total number of runs
    disp(strcat(num2str(n_totalRun),' runs done so far'));
    
    %define the dataset properties
    length = 100;
    area = length^2;
    n = 100;

    %work out the BIC and plot it
    BIC = -2*lnL_array'+log(n)*area*((1:k_max)+1);
    disp(find(min(BIC)==BIC));
    figure;
    plot(BIC);
    xlim([1,k_max]);
    xlabel('Number of factors');
    ylabel('BIC (nat)');

    %plot the individual factor noises
    loading_max = max(max(loading_optimal.^2)); %find the maximum element squared
    loading_min = min(min(loading_optimal.^2)); %find the minimum element squared
    figure;
    %for each factor
    for i = 1:8
        subplot(2,4,i);
        %heat map plot the noise
        imagesc(reshape((loading_optimal(:,i)).^2,length,length),[loading_min,loading_max]);
        set(gca,'xtick',[]); %turn the ticks off
        set(gca,'ytick',[]); %turn the ticks off
    end
    %place the colorbar
    hp4 = get(subplot(2,4,4),'Position');
    colorbar('Position', [0.07  0.04  0.05  hp4(2)+hp4(3)*2.1]);

    %plot the instrinic noise
    figure;
    imagesc(reshape(noise_optimal,length,length));
    colorbar;
    set(gca,'xtick',[]); %turn the ticks off
    set(gca,'ytick',[]); %turn the ticks off
end