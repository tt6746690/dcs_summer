
% added confidence interval for the point estimate
function [Price, CI] = BlsMC2(SO, K, r, T, sigma, NRep1)
    nuT = (r - 0.5*sigma^2)*T;
    siT = sigma * sqrt(T);
    DiscPayoff = exp(-r*T)*max(0, SO*exp(nuT+siT*rand(NRepl, 1))-K);
    [Price, VarPrice, CI] = normfit(DisPayof)
    % normfit: normal paraemter estimates
    % [muHat,sigmaHat,muCI,sigmaCI] = normfit(x, alpha=0.05)
    % returns mean/variance estimator and the corresponding CI given 1-0.05 confidence
end