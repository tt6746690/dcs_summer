clear all;
N = 2500; % Number of creditors
NZ = 5000; % Number of samples from MC
nE = 2; % Number of epsilion samples to take PER z sample
NRuns = 1;
S = 20; % Dimension of Z
C = 4; % Number of credit states

%Initialize data
[H, BETA, tail, EAD, CN, LGC, CMM] = ProblemParams(N, S, true);

disp('BEGIN COMPUTING CREDIT STATE SELECTION MATRIX')
t = cputime;
CSMat = sparse(1:N,CN,ones(N,1));
disp(strcat('FINISH COMPUTING CREDIT STATE SELECTION MATRIX...',num2str(cputime - t),'s'))

disp('BEGIN SAMPLING')
t = cputime;
sampleZ = randn(S,NZ);
disp(strcat('FINISH SAMPLING...',num2str(cputime - t),'s'))

disp('BEGIN COMPUTING PNCZ')
t = cputime;
denom = (1-sum(BETA.^2,2)).^(1/2); %Not used as denom but keeping notation consistant
BZ = BETA*sampleZ;
CH = CSMat*H;
CHZ = repmat(CH,1,1,NZ);
BZ = reshape(BZ,N,1,NZ);
CBZ = repelem(BZ,1,C);
PHI = normcdf((CHZ - CBZ) ./ (denom));
PHI = [zeros(N,1,NZ) PHI];
pncz = diff(PHI,1,2); %column wise diff
clear sampleE;
clear sampleZ;
clear BETA;
clear BZ;
clear CH;
clear CHZ;
clear CBZ;
clear PHI;
disp(strcat('FINISH COMPUTING PNCZ...',num2str(cputime - t),'s'))

disp('BEGIN COMPUTING THETA')
t = cputime;
weights = EAD.*LGC;
[pTheta,theta] = GlassermanPTheta(pncz,weights,tail);
disp(strcat('FINISH COMPUTING THETA...',num2str(cputime - t),'s'))

disp('BEGIN SAMPLING PTHETA')
t = cputime;
cdf = cumsum(pTheta,2);
cdf = repelem(cdf,1,1,nE);
u = rand([N,1,nE*NZ]);
isOne = (cdf >= u) == 1;
ind = isOne & (cumsum(isOne,2) == 1);
disp(strcat('FINISH SAMPLING PTHETA...',num2str(cputime - t),'s'))

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
estimator = exp(-repelem(theta,1,1,nE).*Loss + repelem(psi,1,1,nE)).*(Loss > tail);
l = reshape(estimator,[nE*NZ,1,1]);
vpa(mean(l))
vpa(var(l))
disp(strcat('FINISH COMPUTING LOSS...',num2str(cputime - t),'s'))
%clear all;

