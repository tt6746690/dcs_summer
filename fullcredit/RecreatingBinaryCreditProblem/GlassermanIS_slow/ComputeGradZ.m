function [ gradz ] = ComputeGradZ(pnc, H, BETA, tail, EAD, LGC, z)
    % All gradients are computed with respect to z 
    gPsi = @(pnc,weights,theta,gp) sum(sum(gp.*exp(weights.*theta),2) ./ sum(pnc.*exp(weights.*theta),2),1);
    
    gdPsidThetaL = @(pnc,weights,theta,gp) sum(sum(gp.*weights.*exp(weights.*theta),2) ./ ...
        sum(pnc.*exp(weights.*theta),2),1);
    gdPsidThetaR =  @(pnc,weights,theta,gp) sum(sum(gp.*exp(weights.*theta),2) .* ...
         sum(pnc.*weights.* exp(weights.*theta),2) ./ (sum(pnc.* exp(weights.*theta),2)).^2,1);
    gdPsidTheta = @(pnc,weights,theta,gp) gdPsidThetaL(pnc,weights,theta,gp) - gdPsidThetaR(pnc,weights,theta,gp);
    
    dPsidTheta = @(pnc,weights,theta) sum(sum(pnc.*weights.*exp(weights.*theta),2) ./ sum(pnc.*exp(weights.*theta),2),1);
    
    ddPsiddThetaL = @(pnc,weights,theta) sum(sum(pnc.*(weights.^2).*exp(weights.*theta),2) ./ ...
                       sum(pnc.*exp(weights.*theta),2),1);
    ddPsiddThetaR = @(pnc,weights,theta) sum(((sum(pnc.*weights.*exp(weights.*theta),2)).^2) ./ ...
                        (sum(pnc.*exp(weights.*theta),2).^2),1);
                    
    ddPsiddTheta = @(pnc,weights,theta) ddPsiddThetaL(pnc,weights,theta) - ddPsiddThetaR(pnc,weights,theta);
    
    gTheta = @(pnc,weights,theta,gp) (-1/ddPsiddTheta(pnc,weights,theta)).*gdPsidTheta(pnc,weights,theta,gp);
    
    %grad = @(z,pnc,weights,theta,gp) (-tail + gPsi(pnc,weights,theta,gp) - z + -tail + dPsidTheta(pnc,weights,theta)).*gTheta(pnc,weights,theta,gp);
    grad = @(z,pnc,weights,theta,gp) -tail + gPsi(pnc,weights,theta,gp) - z + (-tail + dPsidTheta(pnc,weights,theta)).*gTheta(pnc,weights,theta,gp);
    % DONE DECLARING UTILITY FUNCTIONS, ACTUALLY TIME TO COMPUTE GRAD
    
    weights = EAD.*LGC;
    
    [N,C] = size(H);
    S = size(BETA,2);
    denom = (1-sum(BETA.^2,2)).^(1/2);
    BZ = BETA*z;
    CH = H;
    CHZ = CH;
    BZ = reshape(BZ,N,1,1);
    CBZ = repelem(BZ,1,C);
    PINV = (CHZ - CBZ) ./ denom;
    PINV = [repelem(-inf,N,1,1) PINV];
    ePINV = exp(-0.5.*PINV.^(2));
    h = (1/sqrt(2*pi)).*diff(ePINV,1,2);
    gp = zeros(N,C,S);
    for n=1:N
       for k=1:C
           gp(n,k,:) = -h(n,k)*BETA(n,:) ./ denom(n);
       end 
    end
    
    [~,theta] = GlassermanPTheta(pnc,weights,tail);
    
    gradz = grad(reshape(z,1,1,S),pnc,weights,theta,gp);
    
    gradz = reshape(gradz,S,1,1);
end

