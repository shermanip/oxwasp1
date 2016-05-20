%FUNCTION: ABC SAMPLES FOR THE CONDITIONAL COMPOUND POSSION
%DESCRIPTION:
    %Y ~ Poisson(poissonParameter) and is latent
    %X|Y ~ N(Y*mean,Y*variance), X is observed
    %ABC draw samples of Y|X=x
    %This is done by drawing samples of Y, then X|Y, accept Y if X=x to some uncertainity
%PARAMETERS:
    %n_samples: positive integer, number of samples
    %x: scalar, the value of the conditional
    %poissonParameter: positive scalar, parameter for the Poisson random variable
    %mean: scalar, mean parameter for observed given latent Normal distribution
    %variance: positive scalar, variance parameter for observed given latent Normal distribution
    %uncertainity: possitive scalar, accept if x-uncertainity < X < x+uncertainity
%RETURN:
    %Y_conditional_array: column vector (size n_samples) containing the samples
function [Y_conditional_array,rejection_array] = ABC_compoundPoisson(n_samples,x,poissonParameter,mean,variance,uncertainity)
    
    %check the parameters are the correct type
    checkParameters(n_samples,x,poissonParameter,mean,variance,uncertainity);

    %declare array of samples
    Y_conditional_array = zeros(n_samples,1);
    %declare array of number of rejections
    rejection_array = zeros(n_samples,1);
    %declare variable for counting the number of accepted samples
    n_accept = 0;
    %declare variable for counting the number of rejections
    n_reject = 0;
    
    %while n_accept is less than required
    while n_accept < n_samples
        %simulated the latent variable
        Y = poissrnd(poissonParameter);
        %simulated the observed given the latent
        X = normrnd(Y*mean,sqrt(Y*variance));
        %if the observed is the same as the required conditional...
        if ( (X>(x-uncertainity)) && (X<(x+uncertainity)) )
            %accept, updated n_accept and save Y to Y_conditional_array
            n_accept = n_accept + 1;
            Y_conditional_array(n_accept) = Y;
            %save the rejection count and reset it
            rejection_array(n_accept) = n_reject;
            n_reject = 0;
        else
            n_reject = n_reject + 1;
        end
    end
    
    %FUNCTION: CHECK THE PARAMETERS ARE OF THE CORRECT TYPE
    function checkParameters(n_samples,x,poissonParameter,mean,variance,uncertainity)
        %n_samples is a positive integer scalar
        if ( (~isscalar(n_samples)) || (n_samples <= 0) || (n_samples ~= floor(n_samples)) )
            error('n_samples needs to be a positive integer scalar');
        end
        %x is a scalar
        if ~isscalar(x);
            error('x needs to be a scalar');
        end
        %poissonParameter is a positive scalar
        if ( (~isscalar(poissonParameter)) || (poissonParameter <= 0) );
            error('poissonParameter needs to be a positive scalar');
        end
        %mean is a scalar
        if ~isscalar(mean);
            error('mean needs to be a scalar');
        end
        %variance is a positive scalar
        if ( (~isscalar(variance)) || (variance <= 0) );
            error('variance needs to be a positive scalar');
        end
        %uncertainity is a positive scalar
        if ( (~isscalar(uncertainity)) || (uncertainity <= 0) );
            error('uncertainity needs to be a positive scalar');
        end
    end

end

