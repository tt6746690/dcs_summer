clear all;
N = 2500; % Number of creditors
NZ = 60000; % Number of samples from MoG (pi*) 
NPi = 600; % Number of samples from MCMC of pi
NRuns = 5; % Number of times to recompute integral before averaging results
S = 10; % Dimension of Z
k = 2; % Number of Gaussians in MoG
burninRatio = 0.1;
C = 4;

a = zeros(1,NRuns);
v = zeros(1,NRuns);

for r=1:NRuns
    totalT = cputime;
    disp(strcat('RUN NUMBER',num2str(r)))
    %Initialize data
    [H, BETA, tail, EAD, CN, LGC, CMM] = ProblemParams(N, S, true);
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
    MoGDen = arrayfun(@(i) EvalMoG(MoGWeights,MoGMu,MoGSigma,sampleZ(:,i)),1:NZ);
    clear MoGIntegrand;
    clear MoGMu;
    clear MoGSigma;
    clear MoGWeights; 
    clear model;
    disp(strcat('FINISH SAMPLING...',num2str(cputime - t),'s'))

    disp('BEGIN COMPUTING LOSS')
    t = cputime;
    l = arrayfun(@(i) f(sampleZ(:,i)')/MoGDen(i),1:NZ);
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
