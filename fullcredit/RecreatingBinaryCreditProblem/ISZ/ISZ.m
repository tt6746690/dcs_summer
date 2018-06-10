clear all;
N = 2500; % Number of creditors
NZ = 400; % Number of samples from MoG (pi*) 
nE = 400; % Number of epsilion samples to take PER z sample
NPi = 1800; % Number of samples from MCMC of pi
NRuns = 5; % Number of times to recompute integral before averaging results
S = 15; % Dimension of Z
k = 2; % Number of Gaussians in MoG
burninRatio = 0.1;

a = zeros(1,NRuns);
v = zeros(1,NRuns);

for r=1:NRuns
    totalT = cputime;
    disp(strcat('RUN NUMBER',num2str(r)))
    %Initialize data
    [H, BETA, tail, EAD, CN, LGC, CMM, C] = ProblemParams(N, S, true);
    %Sample from pi
    disp('BEGIN MCMC SAMPLING FROM PI')
    t = cputime;
    B = floor(NPi * burninRatio);
    f = @(z) DensityAtZ(z,H,BETA,tail,EAD,LGC);
    sampleZ = slicesample(rand(1,S), NPi, 'pdf', f, 'thin', 3, 'burnin', B);
    disp(strcat('FINISH MCMC SAMPLING FROM PI...',num2str(cputime - t),'s'))

    disp('BEGIN TRAINING MOG')
    t = cputime;
    [~, model, ~] = Emgm(sampleZ', k);  
    MoGWeights = model.weight;
    MoGMu = model.mu;
    MoGSigma = model.Sigma;
    disp(strcat('FINISH TRAINING MOG...',num2str(cputime - t),'s'))

    disp('BEGIN SAMPLING')
    t = cputime;
    sampleZ = SampleMoG(MoGWeights,MoGMu,MoGSigma,NZ)';
    %sampleZ = mixGaussRnd(S,k,NZ);
    MoGDen = arrayfun(@(i) EvalMoG(MoGWeights,MoGMu,MoGSigma,sampleZ(:,i)),1:NZ);
    %ZDen = arrayfun(@(i) f(sampleZ(:,i)'),1:NZ);
    %MoGIntegrand = ZDen./MoGDen;
    %disp('MoG Estimate')
    %vpa(mean(MoGIntegrand))
    %vpa(var(MoGIntegrand))
    clear MoGIntegrand;
    clear MoGMu;
    clear MoGSigma;
    clear MoGWeights; 
    clear model;
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

    disp('BEGIN SAMPLING PNCZ')
    t = cputime;
    cdf = cumsum(pncz,2);
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
    weights = EAD.*LGC;
    LossMat = repelem(weights,1,1,NZ*nE).*ind;
    Loss = sum(sum(LossMat,2),1);
    Loss = reshape(Loss,1,nE*NZ);
    LRZ = repelem(arrayfun(@(i) mvnpdf(sampleZ(:,i))/MoGDen(i),1:NZ),1,nE);
    LR = LRZ;
    l = double(Loss > tail).*LR;
    a(r) = vpa(mean(l));
    v(r) = vpa(var(l));
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
    clear LRZ;
    clear MoGDen;
    clear sampleZ;
    clear weights;
    clear ZDen;
    clear pncz;
    disp(strcat('FINISH COMPUTING LOSS...',num2str(cputime - t),'s'))%clear all;
    disp(strcat('TOTAL RUNTIME...',num2str(cputime - totalT),'s'))
end
%[vpa(a); vpa(v)]'
disp('mean')
vpa(a)
vpa(mean(a))
disp('var')
vpa(v)
vpa(mean(v))
