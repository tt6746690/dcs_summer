% the blackwell model for European options
% using a Gaussian distribution for evaluating an expectation
function Price = BlsMC1(SO,K,r,T,sigma,NRepl)
    % NRepl : number of replications, i.e. num_sample
    nuT = (r - 0.5*sigma^2) * T;
    siT = sigma * sqrt(T);
    DiscPayoff = exp(-r*T)*max(0, SO*exp(nuT + siT*randn(NRepl,1))-K);
    % monte carlo integration to estimate expected value
    Price = mean(DiscPayoff); 
end