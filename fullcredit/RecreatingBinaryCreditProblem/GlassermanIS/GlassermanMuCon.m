function [ mu, thetaOpt ] = GlassermanMuCon(z0, theta0, H, BETA, tail, EAD, LGC )

    weights = EAD.*LGC;
    v0 = [z0;theta0];
    function [p] = psi(theta,pnc,weights)
        p = sum(log(sum(pnc.*exp(weights.*theta),2)),1);
    end

    function pncz = ComputePNC(H,BETA,z)
            [N,C] = size(H);
            NZ = 1;
            denom = (1-sum(BETA.^2,2)).^(1/2);
            BZ = BETA*z;
            CH = H;
            CHZ = repmat(CH,1,1,NZ);
            BZ = reshape(BZ,N,1,NZ);
            CBZ = repelem(BZ,1,C);
            PINV = (CHZ - CBZ) ./ denom;
            PHI = normcdf(PINV);
            PHI = [zeros(N,1,NZ) PHI];
            pncz = diff(PHI,1,2); %column wise diff
    end
    
    function [gp] = ComputeGradP(H,BETA,z)
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
    end

    function [energy, grad] = fu(v, H, BETA, tail, weights)
        theta = v(end);
        z = v(1:end-1);
        pnc = ComputePNC(H,BETA,z);
        energy = theta*tail - psi(theta,pnc,weights) + 0.5*(z'*z);
        
        gp = ComputeGradP(H,BETA,z);
        gPsi = sum(sum(gp.*exp(weights.*theta),2) ./ sum(pnc.*exp(weights.*theta),2),1);
        gPsi = reshape(gPsi,size(z));
        
        dpsidtheta = sum(sum(pnc.*weights.*exp(weights.*theta),2) ./ sum(pnc.*exp(weights.*theta),2),1);
        
        grad = [-gPsi + z; ...
            tail - dpsidtheta
            ];
    end

    function [c,ceq,gradc,gradceq] = fl(v, H, BETA, tail, weights)
        c = 0;
        
        theta = v(end);
        z = v(1:end-1);
        pnc = ComputePNC(H,BETA,z);
        dpsidtheta = sum(sum(pnc.*weights.*exp(weights.*theta),2) ./ sum(pnc.*exp(weights.*theta),2),1);
        ceq = -tail + dpsidtheta;
       
        gp = ComputeGradP(H,BETA,z);
        
        gradc = zeros(size(v));
        
        gdPsidThetaL = sum(sum(gp.*weights.*exp(weights.*theta),2) ./ ...
        sum(pnc.*exp(weights.*theta),2),1);
        
        gdPsidThetaR =   sum(sum(gp.*exp(weights.*theta),2) .* ...
         sum(pnc.*weights.* exp(weights.*theta),2) ./ (sum(pnc.* exp(weights.*theta),2)).^2,1);
        
        gdPsidTheta =  gdPsidThetaL - gdPsidThetaR;
        gdPsidTheta = reshape(gdPsidTheta,size(z));
        
        ddPsiddThetaL = sum(sum(pnc.*(weights.^2).*exp(weights.*theta),2) ./ ...
                       sum(pnc.*exp(weights.*theta),2),1);
        ddPsiddThetaR = sum(((sum(pnc.*weights.*exp(weights.*theta),2)).^2) ./ ...
                            (sum(pnc.*exp(weights.*theta),2).^2),1);

        ddPsiddTheta = ddPsiddThetaL - ddPsiddThetaR;

        
        gradceq = [gdPsidTheta;ddPsiddTheta];
    end

    function [hessian] = hessianPZ(z, H, BETA)
        [N,C] = size(H);
        S = size(BETA,2);
        denom = (1-sum(BETA.^2,2)).^(1/2);
        BZ = BETA*z;
        CH = H;
        CHZ = CH;
        BZ = reshape(BZ,N,1,1);
        CBZ = repelem(BZ,1,C);
        PINV = (CHZ - CBZ) ./ denom;
        PINV = [repelem(0,N,1,1) PINV];
        PINV(PINV == inf) = 0;
        ePINV = -1*PINV.*exp(-0.5.*PINV.^(2));
        h = (1/sqrt(2*pi)).*diff(ePINV,1,2);
        hessian = zeros(N,C,S,S);
        for n=1:N
           for k=1:C
               hessian(n,k,1:S,1:S) = h(n,k)*(BETA(n,:)'*BETA(n,:)) ./ (denom(n).^2);
           end 
        end
    end

    function [hessian] = hessianPsiZ(v, H, BETA, tail, weights, pnc)
        theta = v(end);
        z = v(1:end-1);
        S = size(z,1);
        gpz = ComputeGradP(H,BETA,z);
        hpz = hessianPZ(z, H, BETA);
        hLeft = sum(sum(hpz.*exp(weights.*theta),2),1) ./ sum(sum(pnc.*exp(weights.*theta),2),1);
        hLeft = reshape(hLeft,S,S,1,1);
        a = sum(sum(gpz.*exp(weights.*theta),2),1);
        a = reshape(a,S,1,1);
        b = sum(sum(pnc.*exp(weights.*theta),2),1)^.2;
        hRight = (a*a') ./ b;
        hessian = hLeft + hRight;
    end

    function [hessian] = hessianE(v, lambda, H, BETA, tail, weights)
        theta = v(end);
        z = v(1:end-1);
        S = size(z,1);
        pnc = ComputePNC(H,BETA,z);
        gp = ComputeGradP(H,BETA,z);
        hPsiZ = hessianPsiZ(v, H, BETA, tail, weights, pnc);
        hPsiZ = reshape(hPsiZ,S,S,1,1);
        
        gdPsidThetaL = sum(sum(gp.*weights.*exp(weights.*theta),2) ./ ...
        sum(pnc.*exp(weights.*theta),2),1);
        
        gdPsidThetaR =   sum(sum(gp.*exp(weights.*theta),2) .* ...
         sum(pnc.*weights.* exp(weights.*theta),2) ./ (sum(pnc.* exp(weights.*theta),2)).^2,1);
        
        gdPsidTheta =  gdPsidThetaL - gdPsidThetaR;
        gdPsidTheta = reshape(gdPsidTheta,size(z));
        
        ddPsiddThetaL = sum(sum(pnc.*(weights.^2).*exp(weights.*theta),2) ./ ...
                       sum(pnc.*exp(weights.*theta),2),1);
        ddPsiddThetaR = sum(((sum(pnc.*weights.*exp(weights.*theta),2)).^2) ./ ...
                            (sum(pnc.*exp(weights.*theta),2).^2),1);

        ddPsiddTheta = ddPsiddThetaL - ddPsiddThetaR;
        
        hessian = zeros(S+1,S+1);
        hessian(1:S,1:S) = -hPsiZ;
        hessian(S+1,1:S) = -gdPsidTheta;
        hessian(1:S,S+1) = -gdPsidTheta;
        hessian(S+1,S+1) = -ddPsiddTheta;
        
        hessian = hessian - lambda.eqnonlin*hessian; % hessian of constraint is negative of current hessian
        
        hessian(1:S,1:S) = hessian(1:S,1:S) + ones(S,S);
  
    end

    Energy = @(v) fu(v, H, BETA, tail, weights);
    Condition = @(v) fl(v, H, BETA, tail, weights);
    HessianFcn = @(v,lambda) hessianE(v, lambda, H, BETA, tail, weights);
    options = optimoptions ('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',13000,...
        'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'CheckGradients',false,...
        'maxIter',3000); %,'HessianFcn',HessianFcn
    [v,fval,exitflag,output,lambda,grad,hessian] = fmincon(Energy,v0,[],[],[],[],[],[],Condition,options)
    
    thetaOpt = v(end);
    mu = v(1:end-1);

end

