clear all;
N = 2500; % Number of creditors
NZ = 300; % Number of samples from MC
nE = 300; % Number of epsilion samples to take PER z sample
NPi = 600; % Number of samples from MCMC of pi
NRuns = 1; % Number of times to recompute integral before averaging results
S = 5; % Dimension of Z
k = 2; % Number of Gaussians in MoG
burninRatio = 0.1;
C = 4;
  
a = zeros(1,NRuns);
v = zeros(1,NRuns);
  
az = zeros(1,NZ);
vz = zeros(1,NZ);

%Initialize data
[H, BETA, tail, EAD, CN, LGC, CMM, C] = ProblemParams(N, S, true);
for r=1:NRuns
    totalT = cputime;
    disp(strcat('RUN NUMBER',num2str(r)))
    zIndex = 1;
    %Sample from pi
    %disp('BEGIN MCMC SAMPLING FROM PI')
    %t = cputime;
    B = floor(NPi * burninRatio);
    f = @(z) DensityAtZ(z,H,BETA,tail,EAD,LGC);
    sampleZ = slicesample(rand(1,S), NPi, 'pdf', f, 'thin', 3, 'burnin', B);
    %disp(strcat('FINISH MCMC SAMPLING FROM PI...',num2str(cputime - t),'s'))

    %disp('BEGIN TRAINING MOG')
    %t = cputime;
    [~, model, ~] = Emgm(sampleZ', k);  
    MoGWeights = model.weight;
    MoGMu = model.mu;
    MoGSigma = model.Sigma;
    %disp(strcat('FINISH TRAINING MOG...',num2str(cputime - t),'s'))
    for zIndex=1:NZ
    %while true
        %disp('BEGIN SAMPLING')
        t = cputime;
        sampleZ = SampleMoG(MoGWeights,MoGMu,MoGSigma,1)';
        %sampleZ = mixGaussRnd(S,k,1);
        MoGDen = arrayfun(@(i) EvalMoG(MoGWeights,MoGMu,MoGSigma,sampleZ(:,i)),1:1);
        %ZDen = arrayfun(@(i) f(sampleZ(:,i)'),1:1);
        %MoGIntegrand = ZDen./MoGDen;
        %disp('MoG Estimate')
        %vpa(mean(MoGIntegrand))
        %vpa(var(MoGIntegrand))
        %clear MoGIntegrand;
        %clear MoGMu;
        %clear MoGSigma;
        %clear MoGWeights; 
        %clear model;
        %disp(strcat('FINISH SAMPLING...',num2str(cputime - t),'s'))

        %disp('BEGIN COMPUTING PNCZ')
        %t = cputime;
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
        pncz = diff(PHI,1,2); %column wise diff
%         clear sampleZ;
%         clear BETA;
%         clear BZ;
%         clear CH;
%         clear CHZ;
%         clear CBZ;
%         clear PHI;
%         clear PINV;
%         clear denom;
        %disp(strcat('FINISH COMPUTING PNCZ...',num2str(cputime - t),'s'))

        %disp('BEGIN SAMPLING PNCZ')
        %t = cputime;
        cdf = cumsum(pncz,2);
        cdf = repelem(cdf,1,1,nE);
        u = rand([N,1,nE*1]);
        isOne = (cdf >= u) == 1;
        ind = (cumsum(isOne,2) == 1);
        %clear isOne;
        %clear u;
        %clear cdf;
        %disp(strcat('FINISH SAMPLING PNCZ...',num2str(cputime - t),'s'))

        %disp('BEGIN COMPUTING LOSS')
        %t = cputime;
        weights = EAD.*LGC;
        weights = repelem(weights,1,1,nE*1);
        LossMat = weights.*ind;
        Loss = sum(sum(LossMat,2),1);
        Loss = reshape(Loss,1,nE*1);
        LRZ = repelem(arrayfun(@(i) mvnpdf(sampleZ(:,i))/MoGDen(i),1:1),1,nE);
        LR = LRZ;
        l = double(Loss > tail).*LR;
        az(zIndex) = vpa(mean(l));
        vz(zIndex) = vpa(var(l));
        if (mod(zIndex,10000) == 0)
           vpa(mean(az))
        end
        %zIndex = zIndex + 1;

        %disp(strcat('FINISH COMPUTING LOSS...',num2str(cputime - t),'s'))
        %clear C;
        %clear CMM;
        %clear CN;
        %clear EAD;
        %clear H;
        %clear l;
        %clear LGC;
        %clear Loss;
        %clear LossMat;
        %clear pncz;
        %clear weights;
        
    end
    a(r) = vpa(mean(az));
    v(r) = vpa(mean(vz));
    disp(strcat('TOTAL RUNTIME...',num2str(cputime - totalT),'s'))
end

%[vpa(a); vpa(v)]'
disp('mean')
vpa(a)
vpa(mean(a))
disp('var')
vpa(v)
vpa(mean(v))


