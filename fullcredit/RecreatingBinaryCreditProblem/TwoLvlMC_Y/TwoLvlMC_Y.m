% https://github.com/AdamSturge/FullCreditProblem/blob/master/RecreatingBinaryCreditProblem/TwoLvlMC_Y/TwoLvlMC_Y.m
% naive MC integration
clear all;

N       = 2500;     % Number of creditors
NZ      = 100;    % Number of samples from MC
nE      = 100;        % Number of epsilion samples to take PER z sample
S       = 5;       % Dimension of Z
% NRuns   = 5;        % Number of times to recompute integral before averaging results
NRuns   = 1;

a = zeros(1,NRuns);
v = zeros(1,NRuns);
muci = zeros(2,NRuns);     % CI for MC estimate 
sigmaci = zeros(2, NRuns); % CI for MC estimate variance

for r=1:NRuns
    totalT = cputime;
    disp(strcat('RUN NUMBER',num2str(r)))
    [H, BETA, tail, EAD, CN, LGC, CMM, C] = ProblemParams(N, S, true);  %Initialize data

    disp('BEGIN SAMPLING')
    t = cputime;
    % note sampling N(0, I_S) is equivalent to sampling N(0,1) S times, and concatenate the results
    % each column of sample{Z,E} is a valid sample from N(0, I_{S/N})
    sampleZ = randn(S,NZ);          % S x NZ
    sampleE = randn(N,nE*NZ);       % N x NZ*NE  (inner for loop, so needs to *NZ)
    disp(strcat('FINISH SAMPLING...',num2str(cputime - t),'s'))

    disp('BEGIN COMPUTING Y')
    t = cputime;
    denom = (1-sum(BETA.^2,2)).^(1/2);  % (N x S) -> (N x 1)
    BZ = BETA*sampleZ;                  % (N x S)*(S x NZ) -> (N x NZ)
    BZ = reshape(BZ,N,1,NZ);            % (N x 1 x NZ),  (:,:,i) returns i-th column with shape (N x 1)
    BZ = repelem(BZ,1,1,nE);            % (N x 1 x NZ*nE), repeat, because Z stays constant for nE iterations in the same loop
    sampleE = reshape(sampleE,N,1,nE*NZ);       % (N x 1 x NZ*nE), nE*NZ samples, both for loop ordered linearly in the last dimension
    Y = BZ + bsxfun(@times,sampleE,denom);      % (N x 1 x NZ*nE), y_n is Nx1 matrix, different at nE*NZ inner loop body
    clear sampleE;
    clear BETA;
    clear denom;
    clear BZ;
    disp(strcat('FINISH COMPUTING Y...',num2str(cputime - t),'s'))
    
    disp('BEGIN COMPUTING INDICATORS')
    t = cputime;
    CH = H;                                 %  (N x C)
    CHZE = repmat(CH,1,1,nE*NZ);            %  (N x C x NZ*nE), threshold H is same at inner loop body
    isOne = ((Y <= CHZE) == 1);             %  (N x C x NZ*nE), logical matrix, 1 if Y <= CHZE, note Y broadcast to (_ x C x _) 
    % conceptually, at inner loop body, i.e. one of NZ*nE iterations, 
    % we have (N x C) isOne matrix that looks like the following.
    % 1   1   0   0
    % 1   0   0   0
    % 0   0   0   0
    % 1   1   1   0
    % ...
    % cumsum(isOne,2)
    % 1   2   2   2
    % 1   1   1   1
    % 0   0   0   0
    % 1   2   3   3
    % ...
    % ind
    % 1   0   0   0
    % 1   0   0   0
    % 0   0   0   0
    % 1   0   0   0
    % ...
    % basically LGC is non-negtative for credit state 1, so care about indicator_{c(n) -> c=1} only
    ind = isOne & (cumsum(isOne,2) == 1);   %  (N x C x NZ*nE)
    disp(strcat('FINISH COMPUTING INDICATORS...',num2str(cputime - t),'s'))

    disp('BEGIN COMPUTING LOSS')
    t = cputime;
    weights = EAD.*LGC;                         % (N x C)
    LossMat = repelem(weights,1,1,NZ*nE).*ind;  % (N x C x NZ*nE)
    Loss = sum(sum(LossMat,2),1);               % (1 x 1 x NZ*nE)
    Loss = reshape(Loss,1,nE*NZ);               % (1 x NZ*nE),  loss for each inner loop
    l = double(Loss > tail);                    % (1 x NZ*nE),  a indicator, i.e. the integrand for 2 level MC integraion
    % TODO: store normfit(l) instead
    % a(r) = vpa(mean(l));                        % the MC estimate
    % v(r) = vpa(var(l));   
    alpha_ = 0.05;
    [m, s, mci, sci] = normfit(l, alpha_);
    a(r) = vpa(m);
    v(r) = vpa(s);
    muci(:,r) = vpa(mci);
    sigmaci(:,r) = vpa(sci);
    disp(strcat('TOTAL RUNTIME...',num2str(cputime - totalT),'s'))
end


disp('mean')
vpa(a)
vpa(muci)

disp('sqrt(variance)')
vpa(v)
vpa(sigmaci)
