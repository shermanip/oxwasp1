%FUNCTION: CONDITIONAL EXPECTED SUFFICIENT STATISTICS FOR COMPOUND POISSON
%PARAMETER:
    %truncation (positive integer): number of terms in the sum
    %x (scalar): conditional observable
    %m (scalar): mean Normal parameter
    %v (positive scalar): variance Normal parameter
    %psy (positive scalar): Poisson parameter
%RETURN:
    %expectation: E(Y|X=x)
    %zeta: E(1/Y|X=x)
    %Z: p_X(x)
function [expectation,zeta,Z] = expectation_compoundPoisson(truncation,x,m,v,psy)

    %CHECK THE PARAMETERS ARE OF THE CORRECT TYPE
    checkParameters(truncation,x,m,v,psy);

    %if condition that the observed compound poisson varaible is 0...
    if x==0
        %expected poisson variable is 0
        expectation = 0;
        %expected 1/poisson varaible is not defined
        zeta = nan;
        %p_X(x) = infinity
        Z = nan;
    %else for positive x...
    else
        %define the y's to be sumed
        y = (1:truncation);
        
        %define a matrix of terms of the log likelihood
            %columns: for each y = 1,2,3,...
            %rows: each term
        lnz = [
            -psy*ones(truncation,1)';
            y*log(psy);
            -gammaln(y+1);
            -0.5*log(2*pi)*ones(truncation,1)';
            -0.5*log(v)*ones(truncation,1)';
            -0.5*(x-y*m).^2./(v*y)
        ];
        
        %Z = p_X(x) by marginalising the joint distribution
        Z = sum( exp( sum([lnz;-0.5*log(y)],1) ) );
        
        %work out the expectation of Y by truncating the expectation sum
        expectation = sum( exp( sum([lnz;0.5*log(y)],1) ) ) / Z;
        
        %work out the expectation of 1/Y by truncating the expectation sum
        zeta = sum( exp( sum([lnz;-(3/2)*log(y)],1) ) ) / Z;
    end

    %FUNCTION: CHECK THE PARAMETERS ARE OF THE CORRECT TYPE
    function checkParameters(truncation,x,m,v,psy)
        %n_samples is a positive integer scalar
        if ( (~isscalar(truncation)) || (truncation <= 0) || (truncation ~= floor(truncation)) )
            error('truncation needs to be a positive integer scalar');
        end
        %x is a scalar
        if ~isscalar(x);
            error('x needs to be a scalar');
        end
        %m is a scalar
        if ~isscalar(m);
            error('m needs to be a scalar');
        end
        %v is a positive scalar
        if ( (~isscalar(v)) || (v <= 0) );
            error('v needs to be a positive scalar');
        end
        %psy is a positive scalar
        if ( (~isscalar(psy)) || (psy <= 0) );
            error('psy needs to be a positive scalar');
        end
    end
    
end

