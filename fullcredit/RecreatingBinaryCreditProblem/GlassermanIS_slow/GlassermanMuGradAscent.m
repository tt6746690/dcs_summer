function [ mu ] = GlassermanMuGradAscent( mu, H, BETA, tail, EAD, LGC, tol, step, maxIter, iterDecay, decay )

    disp('Initial guess')
    vpa(mu)
    %S = size(mu,1);
    updateSize = 1e5;
    iter = 1;
    %weights = EAD.*LGC;
    while(updateSize > tol && iter < maxIter)
        pnc = ComputePNC(H,BETA,mu);
        gradZ = ComputeGradZ(pnc, H, BETA, tail, EAD, LGC, mu);
  
        updateSize = norm(gradZ);
        update = step*gradZ./updateSize;
        mu = mu + update;
        updateSize
        if(mod(iter,iterDecay) == 0)
            step = decay*step;
        end
        iter = iter + 1;
    end
    disp('iteration count')
    iter
    disp('Final guess')
    vpa(mu)
   
    function pncz = ComputePNC(H,BETA,z)
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

