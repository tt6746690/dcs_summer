clear all;

clear_intermediates = false;

N = 2500; % Number of creditors
NZ = 100; % Number of Z samples 
nE = 100; % Number of epsilion samples to take PER z sample
NRuns = 1; % Number of times to recompute integral before averaging results
S = 5; % Dimension of Z


a = zeros(1,NRuns);
v = zeros(1,NRuns);
muci = zeros(2,NRuns);     % CI for MC estimate 
sigmaci = zeros(2, NRuns); % CI for MC estimate variance


for r=1:NRuns 
    totalT = cputime;
    disp(strcat('RUN NUMBER',num2str(r)))
    %Initialize data
    [H, BETA, tail, EAD, CN, LGC, CMM, C] = ProblemParams(N, S, false);
    
    disp('BEGIN FINDING SHIFTED MEAN')
    t = cputime;
    mu = GlassermanMu(zeros(S,1), H, BETA, tail, EAD, LGC);
    disp('>>>>> fminsearch (un-constrained): '); mu
    [mu,~] = GlassermanMuCon(zeros(S,1),0, H, BETA, tail, EAD, LGC);
    disp('>>>>> fmincon (constrained): '); mu
    [mu,~] = GlassermanMuConHess(zeros(S,1),0, H, BETA, tail, EAD, LGC);
    disp('>>>>> fmincon with hessian (constrained): '); mu
    mu = GlassermanMuGradAscent(mu, H, BETA, tail, EAD, LGC, 0, 0.5, 1000, 200, 0.02)
    disp('>>>>> gradient ascent (constrained): '); mu

    % verify gradient of mu by computing the gradient
    % yield array of length 5 containing -0.2 as gradient for Z
    % pnc = ComputePNC(H,BETA,mu);
    % gradZ = ComputeGradZ(pnc, H, BETA, tail, EAD, LGC, mu);
    % gradZ

    % mu is optimal shifted parameter for sampling Z from N(mu, I_S)
    %   where mu has dimension (Sx1)
    disp(strcat('FINISH FINDING SHIFTED MEAN...',num2str(cputime - t),'s'))
    
    disp('BEGIN SAMPLING')
    t = cputime;
    % sample Z ~ N(mu, I_S) NZ times  ->  100 x 5  
    % --transpose--> 5 x 100, i.e. each column is one sampled Z
    % this is the outer level IS
    sampleZ = mvnrnd(mu,eye(S),NZ)';
    disp(strcat('FINISH SAMPLING...',num2str(cputime - t),'s'))

    disp('BEGIN COMPUTING PNCZ')
    t = cputime;
    % denom(n) = sqrt{1 - \beta_n^T \beta_n}
    %   -- N x 1
    denom = (1-sum(BETA.^2,2)).^(1/2);
    % BZ(n) = \beta_n^T * Z
    %   -- N x NZ  for each obligor n and each sampleZ
    BZ = BETA*sampleZ;
    % BZ reshaped
    %   -- N x 1 x NZ
    BZ = reshape(BZ,N,1,NZ);
    % CBZ second dimension element repeated 1 -> C
    %   -- N x C x NZ
    CBZ = repelem(BZ,1,C);
    CH = H;
    % CHZ is simply NZ copies of transition probability matrix H
    %   -- N x C x NZ
    CHZ = repmat(CH,1,1,NZ);
    % quantiles under corpula model for credit state transition
    %   -- N x C x NZ
    PINV = (CHZ - CBZ) ./ denom;
    % computes cdf of PINV under standard normal
    %   -- N x C x NZ
    PHI = normcdf(PINV);
    %   -- N x C+1 x NZ
    %  0    0.0032    1.0000    1.0000    1.0000
    %  0    0.0006    1.0000    1.0000    1.0000
    PHI = [zeros(N,1,NZ) PHI];
    % probabilities for bernoulli indicator variable
    %   -- N x C x NZ  (diff does column wise difference)
    pncz = diff(PHI,1,2); 

    if clear_intermediates
        clear BETA;
        clear BZ;
        clear CH;
        clear CHZ;
        clear CBZ;
        clear PHI;
        clear PInv;
    end
    disp(strcat('FINISH COMPUTING PNCZ...',num2str(cputime - t),'s'))
    
    disp('BEGIN COMPUTING THETA')
    t = cputime;
    weights = EAD.*LGC;
    % theta is the twisting parameter and pTheta is the twisted
    % distribution for the inner level IS
    [pTheta,theta] = GlassermanPTheta(pncz,weights,tail);
    disp(strcat('FINISH COMPUTING THETA...',num2str(cputime - t),'s'))

    disp('BEGIN SAMPLING PNCZ')
    t = cputime;
    % computes cumulative sum for each row, i.e. each obligor, cross
    % different credit states 
    %     -- N x C x NZ
    cdf = cumsum(pTheta,2);
    % use the same outer Z variable for `nE` number of inner level runs
    cdf = repelem(cdf,1,1,nE);
    % To sample from discrete distribution, Inner level IS
    %     -- sample from `u ~ Unif(0,1)`
    %     -- returns the first weight s.t. the `cdf > u`
    %     -- Do this for both all inner/outer level samples
    u = rand([N,1,nE*NZ]);
    %     -- N x C x NE*NZ
    isOne = (cdf >= u) == 1;
    % cumsum(isOne,2) will yield 1 for the first 1, and 2 for the second..
    ind = (cumsum(isOne,2) == 1);
    % indicator representing
    if clear_intermediates
        clear isOne;
        clear u;
        clear cdf;
    end
    disp(strcat('FINISH SAMPLING PNCZ...',num2str(cputime - t),'s'))

    disp('BEGIN COMPUTING LOSS')
    t = cputime;
    % weights shared in each iteration, mask by value of the indicator,
    % which represents the next credit state
    %     -- N x C x NE*NZ
    LossMat = repelem(weights,1,1,NZ*nE).*ind;
    % sum up the loss for all oligor and credit state, i.e. L_{ij} for 
    % i-th Z sample and j-th E sample
    %     -- 1 x NE*NZ
    Loss = sum(sum(LossMat,2),1);
    % theta needed to compute psi
    theta = reshape(theta,[1,1,NZ]);
    B = zeros([N C NZ]);
    for j=1:NZ
        B(:,:,j) = theta(:,:,j)*weights;
    end
    % psi used to compute the likelihood
    %     -- 1 x 1 x NZ
    psi = sum(log(sum(pncz.*exp(B),2)),1);
    % likelihood for E sampling, inner level
    LRE = reshape(exp(-repelem(theta,1,1,nE).*Loss + repelem(psi,1,1,nE)),1,nE*NZ,1);
    % likelihood for Z sampling, outer level
    %LRZ = repelem(arrayfun(@(i) exp(-mu'*sampleZ(:,i) + 0.5*(mu'*mu)),1:NZ),1,nE);
    LRZ = repelem(arrayfun(@(i) mvnpdf(sampleZ(:,i))/mvnpdf(sampleZ(:,i),mu,eye(S)),1:NZ),1,nE);
    LR = LRE.*LRZ;
    Loss = reshape(Loss,1,nE*NZ);
    l = double(Loss > tail).*LR;
    % l = double(Loss >= tail).*LR;

    % put to store
    alpha_ = 0.05;
    [m, s, mci, sci] = normfit(l, alpha_);
    a(r) = vpa(m);
    v(r) = vpa(s);
    muci(:,r) = vpa(mci);
    sigmaci(:,r) = vpa(sci);
    % a(r) = mean(vpa(l));
    % v(r) = var(vpa(l));

    if clear_intermediates
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
    end
    disp(strcat('FINISH COMPUTING LOSS...',num2str(cputime - t),'s'))
    disp(strcat('TOTAL RUNTIME...',num2str(cputime - totalT),'s'))
end

%[vpa(a); vpa(v)]'
% disp('mean')
% vpa(a)
% vpa(mean(a))
% disp('var')
% vpa(v)
% vpa(mean(v))

disp('mean')
vpa(a)
vpa(muci)

disp('sqrt(variance)')
vpa(v)
vpa(sigmaci)

