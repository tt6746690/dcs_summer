% naive MC integration with bernoulli variable

clear all;

N       = 2500;     % Number of creditors
NZ      = 40000;    % Number of samples from MC
nE      = 2;        % Number of epsilion samples to take PER z sample
S       = 5;       % Dimension of Z
NRuns   = 5;        % Number of times to recompute integral before averaging results
  
a = zeros(1,NRuns);
v = zeros(1,NRuns);
muci = zeros(2,NRuns);     % CI for MC estimate 
sigmaci = zeros(2, NRuns); % CI for MC estimate variance

for r=1:NRuns
    totalT = cputime;
    disp(strcat('RUN NUMBER',num2str(r)))
    %Initialize data
    [H, BETA, tail, EAD, CN, LGC, CMM, C] = ProblemParams(N, S, true);

    disp('BEGIN SAMPLING')
    t = cputime;
    sampleZ = randn(S,NZ);
    disp(strcat('FINISH SAMPLING...',num2str(cputime - t),'s'))

    disp('BEGIN COMPUTING PNCZ')    % PNCZ: p_n^c(z) probability of 1 for the bernoulli indicator
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
    clear sampleZ;
    clear BETA;
    clear BZ;
    clear CH;
    clear CHZ;
    clear CBZ;
    clear PHI;
    clear PINV;
    clear denom;
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
    l = double(Loss > tail);
    alpha_ = 0.05;
    [m, s, mci, sci] = normfit(l, alpha_);
    a(r) = vpa(m);
    v(r) = vpa(s);
    muci(:,r) = vpa(mci);
    sigmaci(:,r) = vpa(sci);
    disp(strcat('FINISH COMPUTING LOSS...',num2str(cputime - t),'s'))
    clear C;
    clear CMM;
    clear CN;
    clear EAD;
    clear H;
    clear l;
    clear LGC;
    clear Loss;
    clear LossMat;
    clear pncz;
    clear weights;
    disp(strcat('TOTAL RUNTIME...',num2str(cputime - totalT),'s'))
end


disp('mean')
vpa(a)
vpa(muci)
disp('sqrt(variance)')
vpa(v)
vpa(sigmaci)

