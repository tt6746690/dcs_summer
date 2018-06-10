% naive MC integration with bernoulli variable

clear all;

N       = 2500;     % Number of creditors
NZ      = 100;    % Number of samples from MC
nE      = 200;        % Number of epsilion samples to take PER z sample
S       = 10;       % Dimension of Z
NRuns   = 10;        % Number of times to recompute integral before averaging results

% N = 2500; % Number of creditors
% NZ = 200; % Number of samples from MC
% nE = 200; % Number of epsilion samples to take PER z sample
% S = 1; % Dimension of Z
% NRuns = 1; % Number of times to recompute integral before averaging results
  
a = zeros(1,NRuns);
v = zeros(1,NRuns);

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
    a(r) = vpa(mean(l));
    v(r) = vpa(var(l));
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
vpa(mean(a))
disp('var')
vpa(v)
vpa(mean(v))


% N       = 2500;     % Number of creditors
% NZ      = 10000;    % Number of samples from MC
% nE      = 2;        % Number of epsilion samples to take PER z sample
% S       = 10;       % Dimension of Z
% NRuns   = 5;        % Number of times to recompute integral before averaging results


% mean
% [ 0.00135, 0.0014, 0.0011, 0.0014, 0.0011]
% 0.00127

% var
% [ 0.0013482449122454364873191501317251, 0.0013981099054950218724818755688943, 0.0010988449422476020706646027136344, 0.0013981099054943891320940130285067, 0.0010988449422473763397722912316112]
% 0.0012684309215459650503621258366138


% mean
% [ 0.00115, 0.0009, 0.0013, 0.0007, 0.0012]
% 0.00105

% var
% [ 0.0011487349367467758476379913190613, 0.00089923496174773912253258556503965, 0.0012983749187465253981804691463253, 0.00069954497724870206860209309596144, 0.001198619930996165303116463007882]
% 0.0010489019450971815913820073262741




% N       = 2500;     % Number of creditors
% NZ      = 100;    % Number of samples from MC
% nE      = 200;        % Number of epsilion samples to take PER z sample
% S       = 10;       % Dimension of Z
% NRuns   = 10;        % Number of times to recompute integral before averaging results




% mean
% [ 0.00035, 0.0016, 0, 0.00125, 0.0003, 0.0012, 0.0009, 0.0002, 0.00035, 0.0001]
% 0.000625

% var
% [ 0.00034989499474991813321098077516069, 0.0015975198759928511092753833366942, 0, 0.00124849992499622363360811672095, 0.00029992499624980359720916034405036, 0.0011986199309959215744680882664852, 0.00089923496174783561652593677493428, 0.00019996999849997173343517375965206, 0.00034989499474984153432749467427243, 0.000099994999749968661384783019663303]
% 0.00062435546777323360813360952903395

% mean
% [ 0.00215, 0.00325, 0.0004, 0.0004, 0.00025, 0.0002, 0.00005, 0.0002, 0, 0.00005]
% 0.00069499999999999997977312427011043

% var
% [ 0.0021454847742392599632688110489198, 0.0032395994799719191048481548023119, 0.00039985999299957354028572931170515, 0.00039985999299958855648581823594157, 0.00024994999749980792045842004078793, 0.00019996999849994679678520659305008, 0.000049999999999990949307946547230941, 0.00019996999849992413695980164600741, 0, 0.000049999999999992765346585460450868]
% 0.00069346942347100039641394353395754