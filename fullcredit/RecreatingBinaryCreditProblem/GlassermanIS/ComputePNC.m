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
