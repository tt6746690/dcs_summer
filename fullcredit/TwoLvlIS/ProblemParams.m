function [H, BETA, tail, EAD, CN, LGC, CMM, C] = ProblemParams(N, S, loadFixed)
tail = 0.3; % cursive l in paper
C = 4; % number of credit states

if(loadFixed)
           filename = strcat(pwd,'\Experiments\S=',num2str(S),'\params.mat');
           load(filename); 
else
    p = 0.01*(1 + sin(16*pi*(1:N)/N));
    
    % Credit state matrix N x C. Contains probs of each creditor to move to each
    % credit state for THIS timestep
    CMM = zeros(N,C);
    CMM(:,1) = 1/3*p;
    CMM(:,2) = 1-p;
    CMM(:,3) = 1/3*p;
    CMM(:,4) = 1/3*p;

    EAD = 0.5 + rand(N, 1);           % exposure of the nth obligor
    EAD = EAD / sum(EAD);
    CN = 2 * ones(N, 1);        % initial credit state of the nth obligor
    BETA = 1/sqrt(S)*(-1 + (2)*rand(N,S));    % sensitivity to each individual risk
    %BETA = repelem(0.01,N,S);

    LGC = zeros(N, C);
    LGC(:,1) = floor(5*(1:N)/N).^2';
    %LGC(:,2) = LGC(:,1);
    LGC(:,3) = LGC(:,1);
    LGC(:,4) = LGC(:,1);

    cumCMM = cumsum(CMM, 2);
    %H_indc = zeros(C, C);
    H = norminv(cumCMM, 0, 1);
    %H_indc(cumCMM<=0) = -1e100;
    %H_indc(cumCMM>=1) = 1e100;
end

end

