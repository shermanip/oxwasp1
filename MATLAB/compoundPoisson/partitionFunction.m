function Z = partitionFunction(truncation,x,m,v,psy)

    y = 1:truncation;
    
    Z = sum( exp( y*log(psy)-gammaln(y+1)-0.5*log(y)-0.5*(x-m*y).^2./(v*y) ) );

end

