function [H, BETA, tail, EAD, CN, LGC, CMM, C] = ProblemParams(N, S, loadFixed)
tail = 0.3; % cursive l in paper
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
           0.0270 0.0125 0.9397 0.0208
           0.0002 0.0000 0.0202 0.9796];

    % CMM = [1.0000 0.0000 0.0000 0.0000
    %        0.0000 0.0000 0.0000 0.0000
    %        0.0000 0.0000 0.0000 0.0000
    %        0.9900 0.0000 0.0000 0.0000];

    % Homogenerous
    EAD = ones(N, 1);           % exposure of the nth obligor
    EAD = EAD / sum(EAD);
    %CN = 4 * ones(N, 1);        % initial credit state of the nth obligor
    CN = randi([1 4],N,1);  
    BETA = 1/sqrt(S)*(-1 + (2)*rand(N,S));    % sensitivity to each individual risk
    %BETA = repelem(0.01,N,S);

    LGC = zeros(N, C);
    LGC(:, 1) = 0.8;
    LGC(:, 2) = 0.5;
    LGC(:, 3) = 0.3;
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

