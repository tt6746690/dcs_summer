function [H, BETA, tail, EAD, CN, LGC, p, CMM] = BinaryCreditProblemParams(N, S, loadFixed)
%ASSUMES C=2
%Credit states explaination on page 12 of Zhe Wang thesis

    tail = 0.20; 
    
    if(loadFixed)
       filename = strcat(pwd,'\Experiments\FixedParams\S=',num2str(S),'\params.mat');
       load(filename); 
    else
        p = 0.01*(1 + sin(16*pi*(1:N)/N));

        CMM = diag(1-p);
        CMM = [CMM; zeros(1,N)];
        CMM = [CMM [p 1]'];

        H = norminv(p, 0, 1);

        LGC = floor(5*(1:N)/N).^2'; % Note: This is ceil() in the paper, but looks like he used floor in the code

        EAD = 0.5 + rand(N, 1);                % exposure of obligors [0.5,1.5]
        EAD = EAD / sum(EAD);                  % normalize to get weights for each obligor

        CN = ones(N, 1);                       % initial credit state of the obligors

        BETA = 1/sqrt(S)*(-1 + (2)*rand(N,S)); % sensitivity to each individual risk [-1/sqrt(S),1/sqrt(S)]
    end
end

