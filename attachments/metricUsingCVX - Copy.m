function [cvx_optval, M_CVX] = metricUsingCVX(normalizedError, Z)
% Optimization using CVX package %
cvx_begin sdp
    [dim,numVec] = size(Z);
    variable B(dim,dim);  
    expression dist(numVec);
    for k = 1:numVec
        dist(k) = Z(:,k)'*B*Z(:,k);
    end
    minimize (norm(sum(dist - normalizedError),2))
    subject to 
        B >= 0
cvx_end
M_CVX = B;
end 