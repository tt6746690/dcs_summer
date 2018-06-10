clear all;
N = 2500; % Number of creditors
NMC = 20000; % Number of samples from MC
nE = 2; % Number of epsilion samples to take PER z sample
NRuns = 1;
S = 5; % Dimension of Z
C = 4; % Number of credit states

%Initialize data
[H, BETA, tail, EAD, CN, LGC, CMM] = ProblemParams(N, S, true);

disp('BEGIN COMPUTING CREDIT STATE SELECTION MATRIX')
t = cputime;
CSMat = sparse(1:N,CN,ones(N,1));
disp(strcat('FINISH COMPUTING CREDIT STATE SELECTION MATRIX...',num2str(cputime - t),'s'))

disp('BEGIN SAMPLING')
t = cputime;
sampleZ = randn(S,NMC);
sampleE = randn(N,nE*NMC);
disp(strcat('FINISH SAMPLING...',num2str(cputime - t),'s'))

disp('BEGIN COMPUTING Y')
t = cputime;
denom = (1-sum(BETA.^2,2)).^(1/2); %Not used as denom but keeping notation consistant
BZ = BETA*sampleZ;
Y = repelem(BZ,1,2) + bsxfun(@times,sampleE,denom);
clear sampleE;
clear sampleZ;
clear BETA;
clear denom;
clear BZ;
disp(strcat('FINISH COMPUTING Y...',num2str(cputime - t),'s'))

disp('BEGIN COMPUTING INDICATORS')
t = cputime;
CH = CSMat*H;
CHZE = repmat(CH,1,1,nE*NMC);
Y = reshape(Y,N,1,nE*NMC);
isOne = ((Y <= CHZE) == 1);
ind = isOne & (cumsum(isOne,2) == 1);
disp(strcat('FINISH COMPUTING INDICATORS...',num2str(cputime - t),'s'))

disp('BEGIN COMPUTING LOSS')
t = cputime;
weights = EAD.*LGC;
weights = repelem(weights,1,1,nE*NMC);
LossMat = weights.*ind;
Loss = sum(sum(LossMat,2),1);
Loss = reshape(Loss,1,nE*NMC);
l = double(Loss > tail);
mean(l)
var(l)
disp(strcat('FINISH COMPUTING LOSS...',num2str(cputime - t),'s'))
%clear all;

