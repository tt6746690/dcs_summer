function [mu] = GlassermanMu(Z0,H,BETA,tail,EAD,LGC)
    options = optimset('LargeScale','off','MaxFunEvals',15000,'Display','final-detailed','FinDiffType','central','TolX',1.0e-15,'TolFun',1.0e-15,'MaxIter',15000);    
    [mu,E,exitflag,output] = fminsearch(@(z) energy(z), Z0, options)
    
    function E = energy(z)
        weights = EAD.*LGC;
        pncz = pnc(z); % for one z sample
        [~,theta] = GlassermanPTheta(pncz,weights,tail);
        psi = @(theta,pnc) sum(log(sum(pnc.*exp(weights.*theta),2)),1);
        E = theta*tail - psi(theta,pncz) + 0.5*(z'*z); %multiply by -1 to transform into minimization problem
    end

    function pncz = pnc(z)
        [N,C] = size(H);
        NZ = 1;
        denom = (1-sum(BETA.^2,2)).^(1/2);
        BZ = BETA*z;
        CH = H;
        CHZ = repmat(CH,1,1,NZ);
        BZ = reshape(BZ,N,1,NZ);
        CBZ = repelem(BZ,1,C);
        PINV = (CHZ - CBZ) ./ denom;
        PHI = normcdf(PINV);
        PHI = [zeros(N,1,NZ) PHI];
        pncz = diff(PHI,1,2); %column wise diff
    end
end