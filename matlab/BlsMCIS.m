% Importance sampling mehtods to price out-of-money vanilla call
function [Price, CI] = BlsMCIS(SO, K, r, T, sigma, NRep1)
    nuT = (r - 0.5*sigma^2)*T;
    siT = sigma * sqrt(T);
    % the mean of Gaussian
    ISnuT = log(K/SO) - 0.5*sigma^2*T;
    % Gaussian noise
    Veps = randn(NRepl, 1);
    VY = ISnuT + isT*Veps;
    % likelihood ratio, f,g are both Gaussian
    ISRatios = exp( (2*(nuT - ISnuT)*VY - nuT^2 + ISnuT^2)/2/siT^2 );
    DiscPayoff = exp(-r*T)*max(0, (SO*exp(VY)-K));
    [Price, VarPrice, CI] = normfit(DiscPayoff.*ISRatios)
end