function [H, BETA, tail, EAD, CN, LGC, CMM, C] = ProblemParams(N, S, loadFixed)
    % Generate parameters for simulating events
    %
    % N: 2500
    %       number of creditors
    % S: 10
    %       dimension of Z, the systematic risk factor with pdf ~ N(0, I_S)
    % loadFixed: 
    %       use a fixed randomly generated set of parameters for consistent benchmarking
    %
    % H:    (N x C)  ... (2500 x 4)
    %       represents threshold for each of n creditors to migrate from c(n) -> c, 
    %       i.e. H_{c(n)}^1 represents threshold for migrating from c(n) to state 1
    % BETA: (N x S) ... (2500 x 10)
    %       each row represents n-th creditor's sensitivity to the systematic risk factor Z
    % tail: (1,) ... 0.20
    %       tail quantile, i.e. P(L > l)
    % EAD:  (N x 1)
    %       exposure at default for each creditor
    % CN:   (N x 1)
    %       initial credit state for each creditor
    %       in this case, initially all at state 2
    % LGC:  (N x C)
    %       percentage gain/loss when creditor n moves to state c
    %       in this case, idea is higher n, higher LGC, i.e. [1 4 9 16] each appear ~2500/4 times, and 5 appear once
    % CMM:  (N x C)
    %       credit state matrix, representing probability of n-th creditor move to credit state c for THIS timestep
    %       in this case, large probability of transition to state 2 and small chance of transition to state 1, i.e. default
    % C:    (1,) ... 4
    %       number of credit states
    %
    %
        tail = 0.20; % cursive l in paper
        C = 4; % number of credit states
    
        if(loadFixed)
                filename = strcat(pwd,'/Experiments/S=',num2str(S),'/params.mat');
                load(filename); 
        else
    
            % (1:N)/N gives N quantiles, equally spaced
            % (1:5)/5 gives
            %    0.2000    0.4000    0.6000    0.8000    1.0000
    
            % We will be putting all creditors in credit state C with prob p of moving to D
            % sum(p) = 0.1
            p = 0.01*(1 + sin(16*pi*(1:N)/N));
            
            % Credit state matrix N x C. Contains probs of each creditor to move to each credit state for THIS timestep
            CMM = zeros(N,C);
            CMM(:,1) = p;
            CMM(:,2) = 1-p;  % no values for column 3,4, why?
            % 0.0096    0.9904         0         0
            % 0.0098    0.9902         0         0
            % 0.0100    0.9900         0         0
            % ...
    
            % EAD: (N x 1)
            %   exposure at default for each creditor
            EAD = 0.5 + rand(N, 1);
            EAD = EAD / sum(EAD);    % why not some smaller value, for N=2500?
            % 0.5770
            % 0.3113
            % 0.5001
            % ...
    
    
            % CN:   (N x 1)
            %       initial credit state for each creditor, initially all at state=2
            CN = 2 * ones(N, 1);
            % 2
            % 2
            % 2 
            % ...
    
    
            % BETA: (N x S) ... (2500 x 10)
            %   each row represents n-th creditor's sensitivity to the systematic risk factor Z
            %   1/\sqrt{S} * Unif(-1,1)
            BETA = 1/sqrt(S)*(-1 + (2)*rand(N,S));
            % -0.1585   -0.0688    0.1718    0.1612    0.3089   -0.1276   -0.2441   -0.2641    0.0603   -0.1755
            % 0.0174   -0.1777    0.2231    0.1167   -0.0419    0.1418    0.0270    0.0691   -0.2446    0.0803
            % 0.0956   -0.2260    0.1298   -0.3138   -0.1727   -0.0103    0.0977    0.0602    0.2281   -0.2874
            % ...
    
    
            % (N x C)
            LGC = zeros(N, C);
            LGC(:,1) = floor(5*(1:N)/N).^2';
    
            % cumsum(A,2) returns the cumulative sum of each row.
            % (N x C)
            cumCMM = cumsum(CMM, 2);
            % i.e. (N=2500, C=4)
            % 0.0096    1.0000    1.0000    1.0000
            % 0.0098    1.0000    1.0000    1.0000
            % 0.0100    1.0000    1.0000    1.0000
            % ...
    
    
            % norminv(p): inverse cdf evaluated at probability p
            %   norminv(0.025) = -1.96
            % H represent threshold for N creditors to migrate from CN -> c
            %   so H_{c(n)}^1 represents threshold for migrating from c(n) to state 1
            % (N x C)
            H = norminv(cumCMM, 0, 1);
            % (N=2500 x C=4)
            % -2.3417       Inf       Inf       Inf
            % -2.3340       Inf       Inf       Inf
            % -2.3263       Inf       Inf       Inf
            % ...
        end
    
    end
    
    