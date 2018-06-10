function [ logProb ] = LogDensityAtZ(z,H,BETA,tail,EAD,LGC)
    S = size(z,2);
    p = creditEvent(BETA,z,H);
    mu = computeMu(LGC,p,EAD);
    sigma = computeSigma(LGC,p,EAD);
    logProb = log((1-normcdf((tail - mu) / sigma))) - 0.5*(z*z') - 0.5*S*log(2*pi);
   
    function [p] = creditEvent(BETA,z,H)
        NZ = 1;
        N = size(BETA,1);
        C = size(H,2);
        denom = (1-sum(BETA.^2,2)).^(1/2); 
        BZ = BETA*z';
        CH = H;
        CHZ = repmat(CH,1,1,NZ);
        BZ = reshape(BZ,N,1,NZ);
        CBZ = repelem(BZ,1,C);
        PINV = (CHZ - CBZ) ./ denom;
        PHI = normcdf(PINV);
        PHI = [zeros(N,1,NZ) PHI];
        p = diff(PHI,1,2); %column wise diff
    end

    function [mu] = computeMu(LGC,p,EAD)
        mu = sum(EAD.*sum(LGC.*p,2));
    end

    function [sigma] = computeSigma(LGC,p,EAD)
        [N,C] = size(LGC);
        A = zeros(N,(C-1)*C/2);
        index = 1;
        for a=1:C
            for b=1:(a-1)
                A(:,index) = ((LGC(:,a) - LGC(:,b)).^2).*p(:,a).*p(:,b);
                index = index + 1;
            end
        end
        B = sum(A,2);
        sigma = sqrt(sum((EAD.^2).*B));
    end
end

