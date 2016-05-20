%FUNCTION: REJECTION SAMPLES FOR THE CONDITIONAL COMPOUND POSSION
%DESCRIPTION:
    %Y ~ Poisson(psy) and is latent
    %X|Y ~ N(Y*m,Y*v), X is observed
    %Rejection sampling draw samples of Y|X=x
    %Proposal is Poisson(Psy)
%PARAMETERS:
    %n_samples: positive integer, number of samples
    %x: scalar, the value of the conditional
    %psy: positive scalar, parameter for the Poisson random variable
    %m: scalar, mean parameter for observed given latent Normal distribution
    %v: positive scalar, variance parameter for observed given latent Normal distribution
%RETURN:
    %Y_conditional_array: column vector (size n_samples) containing the samples
function [y_sample, reject_array] = rejection_compoundPoisson(n_sample,x,psy,m,v)

    %check the parameters are the correct type
    checkParameters(n_sample,x,psy,m,v);
    
    %DECLARE VARIABLES
    %work out the scalar factor for the proposal
    if x^2<(v^2/4/m^2)
        gamma = 1/(1-2*m^2*x^2/v^2);
    else
        gamma = sqrt(v^2+4*m^2*x^2)-v;
    end
    scalarFactor = sqrt(2*m^2/gamma) * exp( -(2*x*m-gamma)^2/(4*v*gamma) ) / (1-exp(-psy));
    %array of samples of Y|X=x
    y_sample = zeros(n_sample,1);
    %array of number of rejected samples per accepted sample
    reject_array = zeros(n_sample,1);

    %Plot the distribution and the proposal
    y_max = 10;
    f = ones(11,1);
    g = ones(11,1);
    for y = 1:y_max
        f(y+1) = psy^y/factorial(y)/sqrt(v*y)*exp(-0.5*(x-y*m)^2/(y*v));
        g(y+1) = psy^y/factorial(y) * scalarFactor;
    end
    figure;
    bar([f,g]);

    %for n_sample times
    for i = 1:n_sample
        
        %SET INITIAL VARIABLES
        %set the boolean flags
        isAccept = 0; %has the sample been accepted
        %count the number of rejections
        rejects = 0;
        
        %while the sample has not been accepted
        while ~isAccept
            
            %sample a non-zero poisson
            acceptPoisson = 0;
            while ~acceptPoisson
                y = poissrnd(psy);
                if y ~=0
                    acceptPoisson = 1;
                end
            end
            
            %target distribution
            f = exp(y*log(psy)-gammaln(y+1)-0.5*log(v*y)-0.5*(x-y*m)^2/(y*v));
            %proposal
            g = scalarFactor * exp(y*log(psy)-gammaln(y+1));

            %check if the proposal is bigger than the target
            if (f>g)
                disp(f);
                disp(g);
                error('Proposal is wrong');
            end
            %do rejection sampling
            if ( (g*rand(1)) < f)
                isAccept = 1;
            else
                %if reject, increment rejection count
                rejects = rejects + 1;
            end
        end
        %save the sample and the rejection count
        y_sample(i) = y;
        reject_array(i) = rejects;
    end

    %FUNCTION: CHECK THE PARAMETERS ARE OF THE CORRECT TYPE
    function checkParameters(n_sample,x,psy,m,v)
        %n_samples is a positive integer scalar
        if ( (~isscalar(n_sample)) || (n_sample <= 0) || (n_sample ~= floor(n_sample)) )
            error('n_samples needs to be a positive integer scalar');
        end
        %x is a scalar
        if ~isscalar(x);
            error('x needs to be a scalar');
        end
        %poissonParameter is a positive scalar
        if ( (~isscalar(psy)) || (psy <= 0) );
            error('psy needs to be a positive scalar');
        end
        %mean is a scalar
        if ~isscalar(m);
            error('m needs to be a scalar');
        end
        %variance is a positive scalar
        if ( (~isscalar(v)) || (v <= 0) );
            error('v needs to be a positive scalar');
        end
    end

end

