clear all;
N = 2500; % Number of creditors
NZ = 5000; % Number of samples from MoG (pi*) 
NRuns = 1; % Number of times to recompute integral before averaging results
S = 30; % Dimension of Z

a = zeros(1,NRuns);
v = zeros(1,NRuns);

for r=1:NRuns
    disp(strcat('RUN NUMBER',num2str(r)))
    %Initialize data
    [H, BETA, tail, EAD, CN, LGC, CMM] = ProblemParams(N, S, true);
    %Sample from pi
    disp('BEGIN SAMPLING Z')
    t = cputime;
    sampleZ = randn(S,NZ);
    disp(strcat('FINISH SAMPLING Z...',num2str(cputime - t),'s'))

    disp('BEGIN COMPUTING LOSS')
    t = cputime;
    f = @(z) DensityAtZ(z,H,BETA,tail,EAD,LGC);
    l = arrayfun(@(i) f(sampleZ(:,i)'),1:NZ);
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
    disp(strcat('FINISH COMPUTING LOSS...',num2str(cputime - t),'s'))%clear all;
end
%[vpa(a); vpa(v)]'
vpa(a)
vpa(mean(a))
vpa(mean(v))
