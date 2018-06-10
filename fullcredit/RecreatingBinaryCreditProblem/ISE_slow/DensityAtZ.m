%pi(z)
%ICS = initial credit states
function [prob] = DensityAtZ(z,H,BETA,tail,EAD,LGC) 
    p = creditEvent(BETA,z,H);
    mu = computeMu(LGC,p,EAD);
    sigma = computeSigma(LGC,p,EAD);
    prob = (1-normcdf((tail - mu) / sigma)) * mvnpdf(z);
    prob = max(prob, 1e-100);
   
    function [p] = creditEvent(BETA,z,H)
        NMC = 1;
        N = size(BETA,1);
        C = size(H,2);
        z = z'; % for some reason rand(S,1) is giving a 1xS vector
        denom = (1-sum(BETA.^2,2)).^(1/2); 
        BZ = BETA*z;
        CH = H;
        CH = [repmat(-Inf,N,1),CH];
        CHZ = repmat(CH,1,1,NMC);
        BZ = reshape(BZ,N,1,NMC);
        CBZ = repelem(BZ,1,C+1);
        PHI = normcdf((CHZ - CBZ) ./ (denom));
        p = diff(PHI,1,2); %column wise diff
        clear sampleE;
        clear sampleZ;
        clear BETA;
        clear BZ;
        clear CH;
        clear CHZ;
        clear CBZ;
        clear PHI;
    end

    function [mu] = computeMu(LGC,p,EAD)
        mu = sum(sum(EAD.*(LGC.*p)));
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
        %pProd = prod(p,1,2);
        sigma = sqrt(((EAD.^2)')*B);
    end
end