clear all;
N = 2500; % Number of creditors
NZ = 1; % Number of samples from MC
nE = 10000; % Number of epsilion samples to take PER z sample
NRuns = 1; % Number of times to recompute integral before averaging results
S = 1; % Dimension of Z
C = 4;
  
a = zeros(1,NRuns);
v = zeros(1,NRuns);
  
az = zeros(1,NZ);
vz = zeros(1,NZ);

try
    for r=1:NRuns
        totalT = cputime;
        disp(strcat('RUN NUMBER',num2str(r)))
        %Initialize data
        [H, BETA, tail, EAD, CN, LGC, CMM] = ProblemParams(N, S, true);

        disp('BEGIN FINDING SHIFTED MEAN')
        t = cputime;
        [mu,~] = GlassermanMuCon(zeros(S,1),0, H, BETA, tail, EAD, LGC);
        disp(strcat('FINISH FINDING SHIFTED MEAN...',num2str(cputime - t),'s'))
        
        zIndex = 1;
        l = zeros(NZ,1);
        %for zIndex=1:NZ
        while true
            %disp('BEGIN SAMPLING')
            %t = cputime;
            sampleZ = mvnrnd(mu,eye(S),1)';
            %disp(strcat('FINISH SAMPLING...',num2str(cputime - t),'s'))

%             disp('BEGIN COMPUTING PNCZ')
%             t = cputime;
            denom = (1-sum(BETA.^2,2)).^(1/2);
            BZ = BETA*sampleZ;
            CH = H;
            CHZ = repmat(CH,1,1,1);
            BZ = reshape(BZ,N,1,1);
            CBZ = repelem(BZ,1,C);
            PINV = (CHZ - CBZ) ./ denom;
            PHI = normcdf(PINV);
            PHI = [zeros(N,1,1) PHI];
            pncz = diff(PHI,1,2); %column wise diff
%             clear BETA;
%             clear BZ;
%             clear CH;
%             clear CHZ;
%             clear CBZ;
%             clear PHI;
%             clear PInv;
%             disp(strcat('FINISH COMPUTING PNCZ...',num2str(cputime - t),'s'))

%             disp('BEGIN COMPUTING THETA')
%             t = cputime;
            weights = EAD.*LGC;
            [pTheta,theta] = GlassermanPTheta(pncz,weights,tail);
%             disp(strcat('FINISH COMPUTING THETA...',num2str(cputime - t),'s'))

%             disp('BEGIN SAMPLING PNCZ')
%             t = cputime;
            cdf = cumsum(pTheta,2);
            cdf = repelem(cdf,1,1,nE);
            u = rand([N,1,nE*1]);
            isOne = (cdf >= u) == 1;
            ind = (cumsum(isOne,2) == 1);
%             clear isOne;
%             clear u;
%             clear cdf;
%             disp(strcat('FINISH SAMPLING PNCZ...',num2str(cputime - t),'s'))

%             disp('BEGIN COMPUTING LOSS')
%             t = cputime;
            LossMat = repelem(weights,1,1,1*nE).*ind;
            Loss = sum(sum(LossMat,2),1);
            theta = reshape(theta,[1,1,1]);
            B = zeros([N C 1]);
            for j=1:1
                B(:,:,j) = theta(:,:,j)*weights;
            end
            psi = sum(log(sum(pncz.*exp(B),2)),1);
            LRE = reshape(exp(-repelem(theta,1,1,nE).*Loss + repelem(psi,1,1,nE)),1,nE*1,1);
            Loss = reshape(Loss,1,nE*1);
            LRZ = repelem(arrayfun(@(i) exp(-mu'*sampleZ(:,i) + 0.5*(mu'*mu)),1:1),1,nE);
            LR = LRE.*LRZ;
            %l(zIndex) = double(Loss > tail).*LR;
            l = double(Loss > tail).*LR;
            az(zIndex) = mean(vpa(l));
            %vz(zIndex) = vpa(var(l));
            
            if (mod(zIndex,10000) == 0)
              vpa(mean(az))
            end
            zIndex = zIndex + 1;
%             clear C;
%             clear CMM;
%             clear CN;
%             clear denom;
%             clear EAD;
%             clear f;
%             clear H;
%             clear ind;
%             clear l;
%             clear LGC;
%             clear Loss;
%             clear LossMat;
%             clear LR;
%             clear LRE;
%             clear LRZ;
%             clear MoGDen;
%             clear psi;
%             clear pTheta;
%             clear sampleZ;
%             clear theta;
%             clear weights;
%             clear ZDen;
%             clear pncz;
%             disp(strcat('FINISH COMPUTING LOSS...',num2str(cputime - t),'s'))
%             disp(strcat('TOTAL RUNTIME...',num2str(cputime - totalT),'s'))
        end
        a(r) = vpa(mean(az));
        %v(r) = vpa(mean(vz));
        disp(strcat('TOTAL RUNTIME...',num2str(cputime - totalT),'s'))
    end
catch ex
    disp(ex)
end
%[vpa(a); vpa(v)]'
disp('mean')
vpa(a)
vpa(mean(a))
%disp('var')
%vpa(v)
%vpa(mean(v))
