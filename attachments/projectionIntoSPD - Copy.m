function A = projectionIntoSPD(B)

% projection into semi positive definite space %
    
    [m,n] = size(B);
    A = zeros(m,m);  
    [V D] = eig(B);
    for i=1:m
        D(i,i) = max(0,D(i,i));
    end;
    A = V*D*inv(V);
end
