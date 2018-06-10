function [ MoGSamples ] = SampleMoG(weights,mu,sigma,NMoG)
        [S,k,~] = size(mu);

        idx = datasample(1:k,NMoG,'Weights',weights);
        MoGSamples = zeros(NMoG,S);
        for i=1:NMoG
            m = mu(:,idx(i));
            s = sigma(:,:,idx(i));
            MoGSamples(i,:) = mvnrnd(m,s,1); %Note: Can speed up by batching on idx   
        end
end

