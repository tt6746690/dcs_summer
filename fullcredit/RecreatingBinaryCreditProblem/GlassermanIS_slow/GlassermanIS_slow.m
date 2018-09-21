clear all;
N = 2500; % Number of creditors
NZ = 1; % Number of samples from MC
nE = 10000; % Number of epsilion samples to take PER z sample
NRuns = 1; % Number of times to recompute integral before averaging results
S = 5; % Dimension of Z
C = 4;


[H, BETA, tail, EAD, CN, LGC, CMM] = ProblemParams(N, S, false);
disp('BEGIN FINDING SHIFTED MEAN')
t = cputime;

% both ran without hessian, but with fmincon and fminunc
[Mu1,obj1,exit1] = GlassermanMuCon(zeros(S,1),0, H, BETA, tail, EAD, LGC, true, false);
[Mu2,obj2,exit2] = GlassermanMuCon(zeros(S,1),0, H, BETA, tail, EAD, LGC, true, true);

% find the mu with the smaller objective, i.e. a better local minima
[~,idx1] = min(obj1);
[~,idx2] = min(obj2);
if min(obj1)<min(obj2)
    mu=Mu1(:,idx1);
else
    mu=Mu2(:,idx2);
end
mu
disp(strcat('FINISH FINDING SHIFTED MEAN...',num2str(cputime - t),'s'))

% create file for outputs ...
FID = fopen(sprintf('compare_methods_S%d_l%0.2f.txt', S,tail), 'w');  
fprintf(FID, 'MU: %s\n', sprintf('%d ', mu))
fprintf(FID, 'algo,mean,variance,S,tail\n')

a = zeros(1,NRuns);
v = zeros(1,NRuns);

muci = zeros(2,NRuns);     % CI for MC estimate 
sigmaci = zeros(2, NRuns); % CI for MC estimate variance
  
az1 = zeros(1,NZ);
az2 = zeros(1,NZ);
vz = zeros(1,NZ);

try
    for r=1:NRuns
        totalT = cputime;
        disp(strcat('RUN NUMBER',num2str(r)))
        
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
            az1(zIndex) = mean(vpa(l));
            % vz(zIndex) = vpa(var(l));

            % with normfit
            % alpha_ = 0.05;
            % [m, s, mci, sci] = normfit(l, alpha_);
            % a(r) = vpa(m);
            % v(r) = vpa(s);
            % muci(:,r) = vpa(mci);
            % sigmaci(:,r) = vpa(sci);

            %%%%%%%%%%%%%%%%% naive way
            % disp('BEGIN NAIVE SAMPLING')
            sampleZ = randn(S,1);
            denom = (1-sum(BETA.^2,2)).^(1/2); 
            BZ = BETA*sampleZ;
            CH = H;
            CHZ = repmat(CH,1,1,1);
            BZ = reshape(BZ,N,1,1);
            CBZ = repelem(BZ,1,C);
            A = (CHZ - CBZ);
            PINV = zeros(N,C);
            for j=1:C
                PINV(:,j) = A(:,j) ./ denom;
            end
            PHI = normcdf(PINV);
            PHI = [zeros(N,1,1) PHI];
            pncz = diff(PHI,1,2);
            cdf = cumsum(pncz,2);
            cdf = repelem(cdf,1,1,nE);
            u = rand([N,1,nE*1]);
            isOne = (cdf >= u) == 1;
            ind = (cumsum(isOne,2) == 1);
            weights = EAD.*LGC;
            weights = repelem(weights,1,1,nE*1);
            LossMat = weights.*ind;
            Loss = sum(sum(LossMat,2),1);
            Loss = reshape(Loss,1,nE*1);
            l = double(Loss > tail);
            az2(zIndex) = mean(vpa(l));
            % vz(zIndex) = var(vpa(l));

            if (mod(zIndex,10) == 0)
                % mc estimate: mean of per Z sampled loss
                fprintf(FID, 'gl,%d,%d,%d,%d\n', mean(az1), var(az1), S, tail);
                fprintf('gl,%d,%d\n', mean(az1), var(az1));
                fprintf(FID, 'nv,%d,%d,%d,%d\n', mean(az2), var(az2), S, tail);
                fprintf('nv,%d,%d\n', mean(az2), var(az2));
                % az
                % vpa(mean(az))
                % vpa(var(az))
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
        % a(r) = vpa(mean(az));
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


