clear all;
N = 2500; % Number of creditors
NZ = 300; % Number of samples from MC
nE = 300; % Number of epsilion samples to take PER z sample
NRuns = 5;
S = 20; % Dimension of Z

a = zeros(1,NRuns);
v = zeros(1,NRuns);

for r=1:NRuns 
    totalT = cputime;
    disp(strcat('RUN NUMBER',num2str(r)))
    %Initialize data
    [H, BETA, tail, EAD, CN, LGC, CMM, C] = ProblemParams(N, S, true);

%     disp('BEGIN COMPUTING CREDIT STATE SELECTION MATRIX')
%     t = cputime;
%     CSMat = sparse(1:N,CN,ones(N,1));
%     disp(strcat('FINISH COMPUTING CREDIT STATE SELECTION MATRIX...',num2str(cputime - t),'s'))

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
    clear CSMat;
    disp(strcat('FINISH COMPUTING PNCZ...',num2str(cputime - t),'s'))

    disp('BEGIN SAMPLING PNCZ')
    t = cputime;
    cdf = cumsum(pncz,2);
    cdf = repelem(cdf,1,1,nE);
    u = rand([N,1,nE*NZ]);
    isOne = (cdf >= u) == 1;
    ind = isOne & (cumsum(isOne,2) == 1);
    clear isOne;
    clear u;
    clear cdf;
    disp(strcat('FINISH SAMPLING PNCZ...',num2str(cputime - t),'s'))

    disp('BEGIN COMPUTING LOSS')
    t = cputime;
    weights = EAD.*LGC;
    weights = repelem(weights,1,1,nE*NZ);
    LossMat = weights.*ind;
    Loss = sum(sum(LossMat,2),1);
    Loss = reshape(Loss,1,nE*NZ);
    l = double(Loss > tail);
    a(r) = vpa(mean(l));
    v(r) = vpa(var(l));
    clear C;
    clear CSMat;
    clear CMM;
    clear CN;
    clear denom;
    clear EAD;
    clear H;
    clear ind;
    clear l;
    clear LGC;
    clear Loss;
    clear LossMat;
    clear sampleZ;
    clear weights;
    clear pncz;
    clear PINV;
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

