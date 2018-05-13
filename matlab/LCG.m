% linear congruential generator LCG
% generates a deterministic sequence of samples from a Unif(0,1)
% USeq : a vector of uniform samples
% ZSeq : a vector of integers generated, seed = ZSeq(-1)
function [USeq, ZSeq] = LCG(a,c,m,seed,N)
    ZSeq = zeros(N, 1);
    USeq = zeros(N, 1);
    for i = 1:N
        seed = mod(a*seed+c, m);
        ZSeq(i) = seed;
        USeq(i) = seed/m;
    end
end

% testing
% a=5; c=3; m=16; seed=7; N=20;
% [USeq, ZSeq] = LCG(a,c,m,seed,N);
% fprintf(1, '%2d %2d %6.4f \n', [(1:N)', ZSeq, USeq]')