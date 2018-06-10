clear all;
N = 2500; % Number of creditors
NZ = 300; % Number of Z samples 
nE = 300; % Number of epsilion samples to take PER z sample
NRuns = 1; % Number of times to recompute integral before averaging results
S = 5; % Dimension of Z

a = zeros(1,NRuns);
v = zeros(1,NRuns);

for r=1:NRuns 
    totalT = cputime;
    disp(strcat('RUN NUMBER',num2str(r)))
    %Initialize data
    [H, BETA, tail, EAD, CN, LGC, CMM, C] = ProblemParams(N, S, true);
    
    disp('BEGIN FINDING SHIFTED MEAN')
    t = cputime;
    [mu,~] = GlassermanMuCon(zeros(S,1),0, H, BETA, tail, EAD, LGC);
    disp(strcat('FINISH FINDING SHIFTED MEAN...',num2str(cputime - t),'s'))
    
    disp('BEGIN SAMPLING')
    t = cputime;
    sampleZ = mvnrnd(mu,eye(S),NZ)';
    disp(strcat('FINISH SAMPLING...',num2str(cputime - t),'s'))

    disp('BEGIN COMPUTING PNCZ')
    t = cputime;
    denom = (1-sum(BETA.^2,2)).^(1/2);
    BZ = BETA*sampleZ;
    CH = H;
    CHZ = repmat(CH,1,1,NZ);
    BZ = reshape(BZ,N,1,NZ);
    CBZ = repelem(BZ,1,C);
    PINV = (CHZ - CBZ) ./ denom;
    PHI = normcdf(PINV);
    PHI = [zeros(N,1,NZ) PHI];
    pncz = diff(PHI,1,2); %column wise diff
    clear BETA;
    clear BZ;
    clear CH;
    clear CHZ;
    clear CBZ;
    clear PHI;
    clear PInv;
    disp(strcat('FINISH COMPUTING PNCZ...',num2str(cputime - t),'s'))
    
    disp('BEGIN COMPUTING THETA')
    t = cputime;
    weights = EAD.*LGC;
    [pTheta,theta] = GlassermanPTheta(pncz,weights,tail);
    disp(strcat('FINISH COMPUTING THETA...',num2str(cputime - t),'s'))

    disp('BEGIN SAMPLING PNCZ')
    t = cputime;
    cdf = cumsum(pTheta,2);
    cdf = repelem(cdf,1,1,nE);
    u = rand([N,1,nE*NZ]);
    isOne = (cdf >= u) == 1;
    ind = (cumsum(isOne,2) == 1);
    clear isOne;
    clear u;
    clear cdf;
    disp(strcat('FINISH SAMPLING PNCZ...',num2str(cputime - t),'s'))

    disp('BEGIN COMPUTING LOSS')
    t = cputime;
    LossMat = repelem(weights,1,1,NZ*nE).*ind;
    Loss = sum(sum(LossMat,2),1);
    theta = reshape(theta,[1,1,NZ]);
    B = zeros([N C NZ]);
    for j=1:NZ
        B(:,:,j) = theta(:,:,j)*weights;
    end
    psi = sum(log(sum(pncz.*exp(B),2)),1);
    LRE = reshape(exp(-repelem(theta,1,1,nE).*Loss + repelem(psi,1,1,nE)),1,nE*NZ,1);
    LRZ = repelem(arrayfun(@(i) exp(-mu'*sampleZ(:,i) + 0.5*(mu'*mu)),1:NZ),1,nE);
    LR = LRE.*LRZ;
    Loss = reshape(Loss,1,nE*NZ);
    l = double(Loss > tail).*LR;
    a(r) = mean(vpa(l));
    v(r) = var(vpa(l));
    clear C;
    clear CMM;
    clear CN;
    clear denom;
    clear EAD;
    clear f;
    clear H;
    clear ind;
    clear l;
    clear LGC;
    clear Loss;
    clear LossMat;
    clear LR;
    clear LRE;
    clear LRZ;
    clear psi;
    clear pTheta;
    clear sampleZ;
    clear theta;
    clear ZDen;
    clear pncz;
    disp(strcat('FINISH COMPUTING LOSS...',num2str(cputime - t),'s'))
    disp(strcat('TOTAL RUNTIME...',num2str(cputime - totalT),'s'))
end

%[vpa(a); vpa(v)]'
disp('mean')
vpa(a)
vpa(mean(a))
disp('var')
vpa(v)
vpa(mean(v))

