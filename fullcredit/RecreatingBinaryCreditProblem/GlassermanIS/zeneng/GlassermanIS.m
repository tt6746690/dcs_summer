clear all;
N = 2500; % Number of creditors
NZ = 10; % Number of Z samples 
nE = 10; % Number of epsilion samples to take PER z sample
Nrun=40; %Number of run times
S = 2; % Dimension of Z

%probability of loss bigger than tail using GlassermanLi
a1 = zeros(1,Nrun);
v1 = zeros(1,Nrun);
L1=[];
%probability of loss bigger than tail using Naive Bilevel
a2 = zeros(1,Nrun);
v2 = zeros(1,Nrun);
L2=[];



%Initialize data

[H, BETA, tail, EAD, CN, LGC, CMM, C] = ProblemParams(N, S, false);

totalT = cputime;
disp(strcat('RUN NUMBER',num2str(1)))
disp('BEGIN FINDING SHIFTED MEAN')
t = cputime;
[Mu1,obj1,exit1] = GlassermanMuCon(zeros(S,1),0, H, BETA, tail, EAD, LGC, true, false);
[Mu2,obj2,exit2] = GlassermanMuCon(zeros(S,1),0, H, BETA, tail, EAD, LGC, true, true);

[~,idx1] = min(obj1);
[~,idx2] = min(obj2);
if min(obj1)<min(obj2)
    mu=Mu1(:,idx1);
else
    mu=Mu2(:,idx2);
end
    
disp(strcat('FINISH FINDING SHIFTED MEAN...',num2str(cputime - t),'s'))

for r=1:Nrun 
    
    
    disp('BEGIN SAMPLING')
    t = cputime;
    sampleZ = mvnrnd(mu,eye(S),NZ)';
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
    disp(strcat('FINISH COMPUTING PNCZ...',num2str(cputime - t),'s'))
    
    disp('BEGIN COMPUTING THETA')
    t = cputime;
    weights = EAD.*LGC;
    [pTheta,theta] = GlassermanPTheta(pncz,weights,tail);
    disp(strcat('FINISH COMPUTING THETA...',num2str(cputime - t),'s'))

    disp('BEGIN SAMPLING PNCZ')
    t = cputime;
    cdf = cumsum(pTheta,2);
    cdf = repelem(cdf,1,1,nE);
    u = rand([N,1,nE*NZ]);
    isOne = (cdf >= u) == 1;
    ind = (cumsum(isOne,2) == 1);
    disp(strcat('FINISH SAMPLING PNCZ...',num2str(cputime - t),'s'))

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
    LRE = reshape(exp(-repelem(theta,1,1,nE).*Loss + repelem(psi,1,1,nE)),1,nE*NZ,1);
    LRZ = repelem(arrayfun(@(i) exp(-mu'*sampleZ(:,i) + 0.5*(mu'*mu)),1:NZ),1,nE);
    LR = LRE.*LRZ;
    Loss = reshape(Loss,1,nE*NZ);
    l1 = double(Loss > tail).*LR;
    L1=[L1 l1];
    a1(r) = mean(vpa(L1));
    v1(r) = var(vpa(L1));
   
    disp(strcat('FINISH COMPUTING LOSS...',num2str(cputime - t),'s'))
    disp(strcat('TOTAL RUNTIME...',num2str(cputime - totalT),'s'))
    
    
    %%%%%%%%%%%%%%%%%naive way
    disp('BEGIN Naive SAMPLING')
    t = cputime;
    sampleZ = randn(S,NZ);
    sampleE = randn(N,nE*NZ);
    disp(strcat('FINISH SAMPLING...',num2str(cputime - t),'s'))

    disp('BEGIN COMPUTING Y')
    t = cputime;
    denom = (1-sum(BETA.^2,2)).^(1/2); %Not used as denom but keeping notation consistant
    BZ = BETA*sampleZ;
    BZ = reshape(BZ,N,1,NZ);
    BZ = repelem(BZ,1,1,nE);
    sampleE = reshape(sampleE,N,1,nE*NZ);
    Y = BZ + bsxfun(@times,sampleE,denom);
    disp(strcat('FINISH COMPUTING Y...',num2str(cputime - t),'s'))
    
    disp('BEGIN COMPUTING INDICATORS')
    t = cputime;
    CH = H;
    CHZE = repmat(CH,1,1,nE*NZ);
    isOne = ((Y <= CHZE) == 1);
    ind = isOne & (cumsum(isOne,2) == 1);
    disp(strcat('FINISH COMPUTING INDICATORS...',num2str(cputime - t),'s'))

    disp('BEGIN COMPUTING LOSS')
    t = cputime;
    weights = EAD.*LGC;
    LossMat = repelem(weights,1,1,NZ*nE).*ind;
    Loss = sum(sum(LossMat,2),1);
    Loss = reshape(Loss,1,nE*NZ);
    l2 = double(Loss > tail);
    L2=[L2 l2];
    a2(r) = vpa(mean(L2));
    v2(r) = vpa(var(L2));
    disp(strcat('TOTAL RUNTIME...',num2str(cputime - totalT),'s'))
    
    
    
end

figure(1)
plot(NZ*nE:NZ*nE:NZ*nE*Nrun,a1,NZ*nE:NZ*nE:NZ*nE*Nrun,a2)
legend('GlassermanLi','Naive')
title('Mean')
xlabel('Run_number')

figure(2)
plot(NZ*nE:NZ*nE:NZ*nE*Nrun,v1,NZ*nE:NZ*nE:NZ*nE*Nrun,v2)
legend('GlassermanLi','Naive')
title('Variance')
xlabel('Run_number')

disp('GlassermanLi mean')
vpa(mean(a1))
disp('GlassermanLi var')
vpa(mean(v1))
disp('naive mean')
vpa(mean(a2))
disp('naive var')
vpa(mean(v2))

