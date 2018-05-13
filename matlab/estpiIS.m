% estimating pi with importance sampling 
% m: number of samples
% L: number of subintervals
% try 
%   estpiIS(1000, 10)
%   estpiIS(1000, 100)
function z = estpiIS(m, L)
    % define left end-points of sub-intervals 
    % x = j:i:k = [j,j+i,j+2*i,...,j+m*i] 
    %   where  m = fix((k-j)/i)
    % i.e. if L=10, s = [0, 1/10, ..., 9/10]
    s = (0:(1/L):(1-1/L));
    % h(x) = \sqrt{1 + x^2} is the integrands
    hvals = sqrt(1 - s.^2);
    % get cumulative probabilities
    % i.e. given A, cumsum(A) = (A(1), A(1)+A(2), ...) 
    cs = cumsum(hvals);
    % sample m times
    for j=1:m 
        % locate the sub-intervals 
        % rand represents x sampled uniformly from [0,1], i.e. f(x)=1
        % want to know which sub-interval x falls in, i.e. which k
        loc = sum(rand * cs(L) > cs) + 1;
        % sample uniformly within sub-interval to compute h(x)
        x = (loc - 1)/L + rand/L;
        % compute q_k = \frac{h(s_k)}{\sum_j h(s_j)}
        %   note q_k is same for any sample x in k-th inteval
        p = hvals(loc) / cs(L);
        % compute estimate
        est(j) = sqrt(1 - x.^2) / (p*L);
    end 
    z = 4*sum(est)/m;
end