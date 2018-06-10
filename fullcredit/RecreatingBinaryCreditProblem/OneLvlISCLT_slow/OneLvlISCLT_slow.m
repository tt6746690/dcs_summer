clear all;
N = 2500; % Number of creditors
NZ = 1400000; % Number of samples from MoG (pi*) 
NPi = 600; % Number of samples from MCMC of pi
NRuns = 3; % Number of times to recompute integral before averaging results
S = 10; % Dimension of Z
k = 2; % Number of Gaussians in MoG
burninRatio = 0.1;
C = 4;

a = zeros(1,NRuns);
v = zeros(1,NRuns);
try
    for r=1:NRuns
        totalT = cputime;
        %disp(strcat('RUN NUMBER',num2str(r)))
        %Initialize data
        [H, BETA, tail, EAD, CN, LGC, CMM] = ProblemParams(N, S, true);
        %Sample from pi
        %disp('BEGIN MCMC SAMPLING FROM PI')
        t = cputime;
        B = floor(NPi * burninRatio);
        f = @(z) DensityAtZ(z,H,BETA,tail,EAD,LGC);

        sampleZ = slicesample(rand(1,S), NPi, 'pdf', f, 'thin', 3, 'burnin', B);
        %disp(strcat('FINISH MCMC SAMPLING FROM PI...',num2str(cputime - t),'s'))

        %disp('BEGIN TRAINING MOG')
        t = cputime;
        [~, model, ~] = Emgm(sampleZ', k);  
        MoGWeights = model.weight;
        MoGMu = model.mu;
        MoGSigma = model.Sigma;
        %disp(strcat('FINISH TRAINING MOG...',num2str(cputime - t),'s'))
        zIndex = 1;
        l = zeros(NZ,1);
        for zIndex=1:NZ
        %while true
            %disp('BEGIN SAMPLING')
            t = cputime;
            sampleZ = SampleMoG(MoGWeights,MoGMu,MoGSigma,1)';
            MoGDen = arrayfun(@(i) EvalMoG(MoGWeights,MoGMu,MoGSigma,sampleZ(:,i)),1:1);
            %disp(strcat('FINISH SAMPLING...',num2str(cputime - t),'s'))

            %disp('BEGIN COMPUTING LOSS')
            t = cputime;
            l(zIndex) = arrayfun(@(i) f(sampleZ(:,i)')/MoGDen(i),1:1);
            %a(zIndex) = vpa(mean(l));
            %v(zIndex) = vpa(var(l));

            %zIndex = zIndex + 1;
            %disp(strcat('FINISH COMPUTING LOSS...',num2str(cputime - t),'s'))%clear all;
            %disp(strcat('TOTAL RUNTIME...',num2str(cputime - totalT),'s'))
        end
        a(r) = vpa(mean(l));
        disp(strcat('TOTAL RUNTIME...',num2str(cputime - totalT),'s'))
    end
catch ex
    disp(ex)
end
%[vpa(a); vpa(v)]'
disp('mean')
vpa(a)
vpa(mean(a))
%disp('var')
%vpa(v)
%vpa(mean(v))
