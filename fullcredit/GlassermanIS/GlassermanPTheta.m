function [pTheta,thetaVec] = GlassermanPTheta(pncz,weights,tail)

    [~,~,NMC] = size(pncz);
    pTheta = pncz;
    thetaVec = zeros(1,NMC);
    psi = @(theta,pnc) sum(log(sum(pnc.*exp(weights.*theta),2)),1);
    for i=1:NMC
        pnc = pncz(:,:,i);
        threshold = sum(sum(weights.*pnc,2),1);
        if tail > threshold
%             cvx_begin
%                 variable theta(1)
%                 minimize(psi(theta,pnc) - tail*theta)
%             cvx_end
            energy = @(theta) psi(theta,pnc) - tail*theta;
            options = optimset('LargeScale','off', 'display', 'off');
            [theta, fval, exitflag, output] = fminunc(energy, 0, options);
            twist = pnc.*exp(weights.*theta);
            s = sum(twist,2);
            pTheta(:,:,i) = bsxfun(@rdivide,twist,s);
            thetaVec(i) = theta;
        end
    end

end