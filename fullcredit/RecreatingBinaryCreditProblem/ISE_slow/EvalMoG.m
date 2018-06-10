%mu : Sxk 
%sigma S*S*k
function [ density ] = EvalMoG(weights,mu,sigma,z)
[~,k] = size(mu);
density = 0;
for i=1:k
    m = mu(:,i);
    s = sigma(:,:,i);
    density = density + weights(i)*mvnpdf(z,m,s);
end

end

