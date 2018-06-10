% https://github.com/AdamSturge/FullCreditProblem/blob/master/RecreatingBinaryCreditProblem/TwoLvlMC_Y/TwoLvlMC_Y.m
% naive MC integration
clear all;

N       = 2500;     % Number of creditors
NZ      = 40000;    % Number of samples from MC
nE      = 2;        % Number of epsilion samples to take PER z sample
S       = 10;       % Dimension of Z
NRuns   = 5;        % Number of times to recompute integral before averaging results
% NRuns   = 1;

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




% N       = 2500;     % Number of creditors
% NZ      = 10000;    % Number of samples from MC
% nE      = 2;        % Number of epsilion samples to take PER z sample
% S       = 10;       % Dimension of Z
% NRuns   = 5;        % Number of times to recompute integral before averaging results


% mean
% [ 0.0013, 0.00075, 0.00105, 0.00095, 0.00045]

% [ 0.00080058746728925455661551868757897, 0.00037056451333258432368628065844973, 0.0006011132763899369952156903629259, 0.00052300218729800171130889241410955, 0.00015604641963715940207457766319976]
% [  0.0017994125327107451071950716681158,  0.0011294354866674155991185246605824, 0.0014988867236100627666067763854585,   0.001376997812701998176454498690191, 0.00074395358036284046487213172937913]

% sqrt(variance)
% [ 0.03603296988517949023256292662154, 0.027376540573070085732299006053836, 0.032387496777271547465648637853519, 0.030808196267353243208431123889568, 0.021208960158374146054427455965197]

% [ 0.035683297646309772299666462913592, 0.027110872317991547425508613855527, 0.032073201049065351575517723858866, 0.030509226435039362573808929823826, 0.021003143524155280102982956691449]
% [ 0.036389611571505824094696635029322, 0.027647503960958264707103992918746, 0.032708056850541712690105811134345, 0.031113124978599457359784707932704, 0.021418879000484954133742832027565]


% reduce NZ: 10000 -> 10. 
% expect to see an increase in MC estimate variance
% conclusion: actually observed this...
% N       = 2500;     % Number of creditors
% NZ      = 10;    % Number of samples from MC
% nE      = 2;        % Number of epsilion samples to take PER z sample
% S       = 10;       % Dimension of Z
% NRuns   = 5;        % Number of times to recompute integral before averaging results


% mean
% [ 0, 0, 0.05, 0, 0]
% [ 0, 0, -0.054651202720415417712906958058738, 0, 0]
% [ 0, 0,   0.15465120272041543714180988899898, 0, 0]

% sqrt(variance)
% [ 0, 0, 0.22360679774997896964091736687313, 0, 0]
% [ 0, 0, 0.17005082179572186817928525215393, 0, 0]
% [ 0, 0, 0.32659374643706590157776759042463, 0, 0]



% NZ 10000-> 20000 ; nE 2->1
% N       = 2500;     % Number of creditors
% NZ      = 20000;    % Number of samples from MC
% nE      = 1;        % Number of epsilion samples to take PER z sample
% S       = 10;       % Dimension of Z
% NRuns   = 5;        % Number of times to recompute integral before averaging results

% mean
% [ 0.00105, 0.0012, 0.00125, 0.001, 0.00075]
% [ 0.00060111327638996410027000250053675, 0.00072015593726996215539343237921344, 0.00076027346258614072853354581837948, 0.00056192045059451573366104959461609, 0.00037056451333259505728778826494363]
% [  0.0014988867236100356615524642478476,  0.0016798440627300375259378650838471,  0.0017397265374138591066677239638238,  0.0014380795494054840911318793317264,  0.0011294354866674049739372343026389]

% sqrt(variance)
% [ 0.032387496777269590697567736015117, 0.034621090840641632890850587500609, 0.035334118426754372588316499559369, 0.031607751462223519778316216388703, 0.027376540573069315515075672351486]
% [ 0.032073201049063415624118533742148,  0.03428511980120371033242676617192, 0.034991228003396490642273874982493, 0.031301022562145279504175476859018, 0.027110872317990784147179184060406]
% [ 0.032708056850539735105343197574257, 0.034963758243834006100847489051375, 0.035683843126681345780948362289564, 0.031920593890099754896816364180268, 0.027647503960957487550986755309168]



% N       = 2500;     % Number of creditors
% NZ      = 40000;    % Number of samples from MC
% nE      = 2;        % Number of epsilion samples to take PER z sample
% S       = 10;       % Dimension of Z
% NRuns   = 5;        % Number of times to recompute integral before averaging results

% mean
% [ 0.0010625, 0.0013625, 0.001075, 0.0010125, 0.000875]
% [ 0.00083674051252920515417410562974965,  0.001106886037179408274408465295835, 0.00084791781939579395285799812143068, 0.00079211099423860979393263725967245, 0.0006701073770237360970047912900327]
% [  0.0012882594874707949984815602562094, 0.0016181139628205917172648620194764,  0.0013020821806042060485297806593508,  0.0012328890057613900964328390585933, 0.0010798926229762639394244017054802]

% sqrt(variance)
% [ 0.032578894411156966715026328529348, 0.03688713328527002538459456104647, 0.03276976957841804677373787058059, 0.031803891072391947403286138751355, 0.029567639457155577009661584497735]
% [ 0.032420043798213678265529580357907, 0.036707276238614708341323478180129, 0.032609988282041105223107990696008, 0.031648819278761619600004451058339,  0.02942347134660585478149741334164]
% [ 0.032739320215629384713817984220441, 0.037068773827066435477828321154448, 0.032931135295150071318204254566808, 0.031960500586092760455514394379861, 0.029713237164873985973745362798581]