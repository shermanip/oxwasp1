%SCRIPT: FIT COMPOUND POISSON USING EM ALGORITHM FOR DIFFERENT INITIAL PARAMETERS
%DETAILS:
    %The data are simulated 50 times with known parameter. The EM algorithm
    %was initalised by assigning random Poisson distributed and was used to
    %fit
    
%set random seed
rng(64623562);

%parameters
n_sample = 100; %sample size of the simulated data
alpha = 1; %normal scale parameter
time_exposure = 1; %rate scale parameter

mu_actual = 3; %mean parameter
var_actual = 0.1; %variance parameter
rate_actual = 5; %rate parameter

n_EM = 300; %number of EM steps
n_repeat = 50; %number of initial starting values for each Poisson parameter

%set seeds for each core
seeds = randi([0,intmax],n_repeat,1,'int32');

%declare array which stores results for every EM step and each Poisson parameter
lnL_array = zeros(n_EM,n_repeat); %log likelihood
m_array = zeros(n_EM,n_repeat); %mean parameter
v_array = zeros(n_EM,n_repeat); %variance parameter
rate_array = zeros(n_EM,n_repeat); %rate parameter

%array of Poisson parameters to investigate
poissonParameter = 4;

%repeat n_repeat times
parfor i = 1:n_repeat
    %reset the random seed
    rng(seeds(i));
    
    %simulate the data
    [X,~] = simulateData(n_sample,alpha*mu_actual/time_exposure,alpha^2*var_actual/time_exposure^2,rate_actual*time_exposure);
    
    %run EM and get array of mean, variance, rate parameters and log likelihood
    [m,v,rate,lnL] = EM_compoundPoisson(X,poissonParameter,n_EM,[1,100],alpha,time_exposure);
    m_array(:,i) = m'; %save the mean parameters
    v_array(:,i) = v'; %save the variance parameters
    rate_array(:,i) = rate'; %save the rate parameters
    lnL_array(:,i) = lnL'; %save the log likelihood
end

%plot the mean parameter for each EM step
figure;
plot(m_array,'r');
hold on;
plot([0,n_EM],[mu_actual,mu_actual],'k--','LineWidth',2);
hold off;
xlabel('Number of EM Steps');
ylabel('Mean estimate (arb. unit)');

%plot the variance parameter for each EM step
figure;
semilogy(v_array,'r');
hold on;
semilogy([0,n_EM],[var_actual,var_actual],'k--','LineWidth',2);
hold off;
xlabel('Number of EM Steps');
ylabel('Variance estimate (arb. unit)');

%plot the rate parameter for each EM step
figure;
plot(rate_array,'r');
hold on;
plot([0,n_EM],[rate_actual,rate_actual],'k--','LineWidth',2);
hold off;
xlabel('Number of EM Steps');
ylabel('Rate estimate (arb. unit)');

%plot the likelihood for each EM step
figure;
plot(lnL_array,'b');
xlabel('Number of EM Steps');
ylabel('Log likelihood (nats)');