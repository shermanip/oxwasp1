function expectation = expectation_compoundPoisson(truncation,x,m,v,psy)


    y = 1:truncation;
    Z = partitionFunction(truncation,x,m,v,psy);

    expectation = sum( exp( y*log(psy)-gammaln(y)-0.5*log(y)-0.5*(x-m*y).^2./(v*y) ) ) / Z;


end

