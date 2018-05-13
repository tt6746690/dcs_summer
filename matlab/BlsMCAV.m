% added confidence interval for the point estimate
function [Price, CI] = BlsMCAV(SO, K, r, T, sigma, NRep1)
    nuT = (r - 0.5*sigma^2)*T;
    siT = sigma * sqrt(T);
    % (Veps, 1-Veps) is the negatively correlated input stream
    Veps = rand(NRepl, 1);
    Payoff1 = max(0, SO*exp(nuT+siT*Veps) - K);
    Payoff2 = max(0, SO*exp(nuT+siT*(-Veps)) - K);
    DiscPayoff = exp(-r*T) * 0.5 * (Payoff1+Payoff2);
    [Price, VarPrice, CI] = normfit(DiscPayoff);
end