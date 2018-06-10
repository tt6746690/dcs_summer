clear all;
N = 2500; % Number of creditors
NZ = 500; % Number of Z samples 
nE = 500; % Number of epsilion samples to take PER z sample
NRuns = 5; % Number of times to recompute integral before averaging results
S = 20; % Dimension of Z

a = zeros(1,NRuns);
v = zeros(1,NRuns);

for r=1:NRuns 
    totalT = cputime;
    disp(strcat('RUN NUMBER',num2str(r)))
    %Initialize data
    [H, BETA, tail, EAD, CN, LGC, CMM, C] = ProblemParams(N, S, true);

    disp('BEGIN SAMPLING')
    t = cputime;
    sampleZ = randn(S,NZ);
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
    Loss = reshape(Loss,1,nE*NZ);
    LR = LRE;
    l = double(Loss > tail).*LR;
    a(r) = vpa(mean(l));
    v(r) = vpa(var(l));
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
    clear MoGDen;
    clear psi;
    clear pTheta;
    clear sampleZ;
    clear theta;
    clear weights;
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

