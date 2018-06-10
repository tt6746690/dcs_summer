function [H, BETA, tail, EAD, CN, LGC, CMM] = ProblemParams(N, S, loadFixed)
tail = 0.1; % cursive l in paper
C = 4; % number of credit states

if(loadFixed)
           filename = strcat(pwd,'\Experiments\S=',num2str(S),'\params.mat');
           load(filename); 
else
    % credit states explaination C x C
    % probability of transition from "row" state to "col" state
    % 1 - "D" default
    % 2 - "C" grade
    % 3 - "B" grade
    % 4 - "A" grade
    CMM = [1.0000 0.0000 0.0000 0.0000
           0.2550 0.6801 0.0649 0.0000
           0.5000 0.0250 0.1250 0.1250
           0.0008 0.4790 0.0202 0.5000]; %0.9796

    % CMM = [1.0000 0.0000 0.0000 0.0000
    %        0.0000 0.0000 0.0000 0.0000
    %        0.0000 0.0000 0.0000 0.0000
    %        0.9900 0.0000 0.0000 0.0000];

    EAD = 0.5 + rand(N, 1);           % exposure of the nth obligor
    EAD = EAD / sum(EAD);
    CN = 2 * ones(N, 1);        % initial credit state of the nth obligor
    BETA = 1/sqrt(S)*(-1 + (2)*rand(N,S));    % sensitivity to each individual risk
    %BETA = repelem(0.01,N,S);

    LGC = zeros(N, C);
    LGC(:, 1) = 0.1;
    LGC(:, 2) = 0.2;
    LGC(:, 3) = 0.05;
    % LGC(:, 1) = 0.2;
    % LGC(:, 2) = 0.2;
    % LGC(:, 3) = 0.2;

    cumCMM = cumsum(CMM, 2);
    %H_indc = zeros(C, C);
    H = norminv(cumCMM, 0, 1);
    %H_indc(cumCMM<=0) = -1e100;
    %H_indc(cumCMM>=1) = 1e100;
end

end

