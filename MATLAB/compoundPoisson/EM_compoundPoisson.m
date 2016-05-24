%FUNCTION: EM FOR COMPOUND POISSON
%PARAMETERS
    %X (column vector): vector of observables
    %initialPoisson_parameter (positive scalar): initial latent variables are Poisson distributed with this parameter
    %n_EM (positive integer): number of EM steps
    %alpha (positive scalar): scaling of the Normal variables
    %time_exposure (positive scalar): scaling of the estimated Poisson parameter
%RETURN
    %mean_array (column vector): estimated mean Normal parameter at each EM step
    %var_array (column vector): estimated variance Normal parameter at each EM step
    %rate_array (column vector): estimated Poisson rate at each EM step
    %lnL_array (column vector): log likelihood at each EM step
function [mean_array,var_array,rate_array,lnL_array] = EM_compoundPoisson(X,initialPoisson_parameter,n_EM,alpha,time_exposure)

    %check parameters are of the correct type
    checkParameters(X,initialPoisson_parameter,n_EM,alpha,time_exposure);

    %get the sample size
    n_sample = numel(X);

    %INITIAL E STEP
        %assign the latent variables Poisson randomly
    %Y is a column vector of the latent variables
    Y = poissrnd(initialPoisson_parameter,n_sample,1);
    %Zeta is a column vector of E[1/Y|X=x], initially set it 1./Y
    Zeta = Y;
    Zeta(Y==0) = nan; %for Y=0, set nan;
    Zeta(Y~=0) = 1./(Y(Y~=0));

    %declare array storing variables at each EM step
    mean_array = zeros(n_EM,1); %mean Normal parameter
    var_array = zeros(n_EM,1); %variance Normal parameter
    rate_array = zeros(n_EM,1); %rate Poisson parameter
    lnL_array = zeros(n_EM,1); %log likelihood
    
    %temporary variable, likelihood of a single data point
    L = zeros(n_sample,1);
    
    %for n_EM times
    for i = 1:n_EM
        
        %M STEP
        rate = sum(Y)/n_sample/time_exposure; %estimate the rate
        m = time_exposure * sum(X) / alpha / sum(Y); %estimate the mean
        %estimate the variance
        v = ( time_exposure^2/alpha^2*sum(((X(Y~=0)).^2).*Zeta(Y~=0)) - 2*time_exposure*m/alpha*sum(X(Y~=0)) + m^2*sum(Y) ) /sum(Y~=0);

        %save the estimated parameters for this EM step
        mean_array(i) = m;
        var_array(i) = v;
        rate_array(i) = rate;

        %E STEP
        %foe each data, estimate the expected sufficient statistic
        for j = 1:n_sample
            %get the expected sufficient statistic [y = E[Y|X=x], zeta = E[1/Y|X=x], l = p_X(x)
            [y,zeta,l] = expectation_compoundPoisson(100,X(j),alpha*m/time_exposure,alpha^2*v/time_exposure^2,rate*time_exposure);
            %save the sufficient statistic and likelihood
            Y(j) = y;
            Zeta(j) = zeta;
            if isnan(l) %if the likelihood is not aviable (this happens when x=0)...
                l = 1; %set it to 1
            end
            L(j) = l; %save the liklihood
        end
        
        %save the log likelihood of the full dataset
        lnL_array(i) = sum(log(L));
    end
    
    %FUNCTION: CHECK THE PARAMETERS ARE OF THE CORRECT TYPE
    function checkParameters(X,initialPoisson_parameter,n_EM,alpha,time_exposure)
        %X is a column vector
        if ( ~iscolumn(X) )
            error('X needs to be a column vector');
        end
        %initialPoisson_parameter is a positive scalar
        if ( (~isscalar(initialPoisson_parameter)) || (initialPoisson_parameter <= 0) );
            error('initialPoisson_parameter needs to be a positive scalar');
        end
        %n_EM is a positive integer
        if ( (~isscalar(n_EM)) || (n_EM <= 0) || (n_EM ~= floor(n_EM)) )
            error('n_EM needs to be a positive integer scalar');
        end
        %alpha is a positive scalar
        if ( (~isscalar(alpha)) || (alpha <= 0) );
            error('alpha needs to be a positive scalar');
        end
        %time_exposure is a positive scalar
        if ( (~isscalar(time_exposure)) || (time_exposure <= 0) );
            error('time_exposure needs to be a positive scalar');
        end
    end

end

