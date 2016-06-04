%SCRIPT: FIT COMPOUND POISSON USING EM ALGORITHM FOR DIFFERENT INITIAL PARAMETERS
%DETAILS:
    %The data are simulated with known parameters.
    %EM algorithm was initalised by assigning random Poisson distributed
    %latent variables with 3 different parameters (1,5,9). This was done 50
    %times for each parameter.
    %The parameters and log likelihood are plotted for each EM step.

%set random seed
rng(137721447);

%parameters
n_sample = 100; %sample size of the simulated data
alpha = 1; %normal scale parameter
time_exposure = 1; %rate scale parameter

mu_actual = 3; %mean parameter
var_actual = 0.1; %variance parameter
rate_actual = 5; %rate parameter

%simulate the data
[X,Y_actual] = simulateData(n_sample,alpha*mu_actual/time_exposure,alpha^2*var_actual/time_exposure^2,rate_actual*time_exposure);

n_EM = 300; %number of EM steps
n_repeat = 50; %number of initial starting values for each Poisson parameter

%set seeds for each core
seeds = randi([0,intmax],n_repeat,3,'int32');

%declare array which stores results for every EM step and each Poisson parameter
lnL_array = zeros(n_EM,n_repeat,3); %log likelihood
m_array = zeros(n_EM,n_repeat,3); %mean parameter
v_array = zeros(n_EM,n_repeat,3); %variance parameter
rate_array = zeros(n_EM,n_repeat,3); %rate parameter

%array of Poisson parameters to investigate
poissonParameters = [1,5,9];
%for each Poisson parameter
for j = 1:3
    %repeat n_repeat times
    parfor i = 1:n_repeat
        %reset the random seed
        rng(seeds(i,j));
        %run EM and get array of mean, variance, rate parameters and log likelihood
        [m,v,rate,lnL] = EM_compoundPoisson(X,poissonParameters(j),n_EM,[1,100],alpha,time_exposure);
        m_array(:,i,j) = m'; %save the mean parameters
        v_array(:,i,j) = v'; %save the variance parameters
        rate_array(:,i,j) = rate'; %save the rate parameters
        lnL_array(:,i,j) = lnL'; %save the log likelihood
    end
end

%plot the mean parameter for each EM step
figure('Position', [400, 400, 400, 300]);
plot(m_array(:,:,1),'b');
hold on;
plot(m_array(:,:,2),'r');
plot(m_array(:,:,3),'g');
plot([0,n_EM],[mu_actual,mu_actual],'k--','LineWidth',2);
hold off;
xlabel('Number of EM Steps');
ylabel('Mean estimate (arb. unit)');

%plot the variance parameter for each EM step
figure('Position', [400, 400, 400, 300]);
semilogy(v_array(:,:,1),'b');
hold on;
semilogy(v_array(:,:,2),'r');
semilogy(v_array(:,:,3),'g');
semilogy([0,n_EM],[var_actual,var_actual],'k--','LineWidth',2);
hold off;
xlabel('Number of EM Steps');
ylabel('Variance estimate (arb. unit)');

%plot the rate parameter for each EM step
figure('Position', [400, 400, 400, 300]);
plot(rate_array(:,:,1),'b');
hold on;
plot(rate_array(:,:,2),'r');
plot(rate_array(:,:,3),'g');
plot([0,n_EM],[rate_actual,rate_actual],'k--','LineWidth',2);
hold off;
xlabel('Number of EM Steps');
ylabel('Rate estimate (arb. unit)');

%plot the likelihood for each EM step
figure('Position', [400, 400, 400, 300]);
plot(lnL_array(:,:,1),'b');
hold on;
plot(lnL_array(:,:,2),'r');
plot(lnL_array(:,:,3),'g');
hold off;
xlabel('Number of EM Steps');
ylabel('Log likelihood (nats)');