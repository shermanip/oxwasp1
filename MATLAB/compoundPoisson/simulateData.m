%FUNCTION: SIMULATE COMPOUND POISSON DATA
%PARAMETERS:
    %n_sample (positive integer): number of simulated data
    %m (scalar): mean of Normal
    %v (positive scalar): variance of Normal
    %poissonParameter (positive scalar): parameter for the Poisson distribution
%RETURN:
    %X (column vector): n_sample observable compound poisson variables
    %Y (column vector): n_sample latent poisson variables
function [X,Y] = simulateData(n_sample,m,v,poissonParameter)

    %check the parameters are of the correct type
    checkParameters(n_sample,m,v,poissonParameter);

    %simulate the latent Poisson variables
    Y = poissrnd(poissonParameter,n_sample,1);
    %given the latent Poisson variables, simulate the Normal observables
    X = normrnd(Y*m,sqrt(Y*v));

    %FUNCTION: CHECK THE PARAMETERS ARE OF THE CORRECT TYPE
    function checkParameters(n_sample,m,v,poissonParameter)
        %n_samples is a positive integer scalar
        if ( (~isscalar(n_sample)) || (n_sample <= 0) || (n_sample ~= floor(n_sample)) )
            error('n_sample needs to be a positive integer scalar');
        end
        %m is a scalar
        if ~isscalar(m);
            error('m needs to be a scalar');
        end
        %v is a positive scalar
        if ( (~isscalar(v)) || (v <= 0) );
            error('v needs to be a positive scalar');
        end
        %poissonParameter is a positive scalar
        if ( (~isscalar(poissonParameter)) || (poissonParameter <= 0) );
            error('poissonParameter needs to be a positive scalar');
        end
    end

end