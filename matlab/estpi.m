% function for estimating pi
% \pi = \int_0^1 \sqrt{1 - x^&2} dx = \frac{\pi}P{4}
% m: number of samples
function out=estpi(m)
    z = sqrt(1-rand(1,m).^2);
    out = 4*sum(z)/m;
end