function objFuncVal = computeObjectiveFunc( Z, B, errVec )
    
    [dim, numVec] = size(Z);
    dist = zeros(numVec,1);
    for k = 1:numVec
        dist(k) = Z(:,k)'*B*Z(:,k);
    end
    objectiveFunc = norm (sum(dist-errVec), 2); 
end

