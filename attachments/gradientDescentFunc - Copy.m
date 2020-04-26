function [ optValGradientDescent, A] = gradientDescentFunc( Z, error, maxIter, maxLoss, step_size)
%% gradient descent step %%
optValGradientDescent = 10;
[dim,numVec] = size(Z);
numIteration = 0;
A = eye(dim);
%A = generateSPDmatrix(dim);

while (numIteration <= maxIter && optValGradientDescent >= maxLoss) % 
    numIteration = numIteration + 1;   
    grad = zeros(dim, dim);
    for k = 1:numVec        
        %alpha = 2*(diag(Z'*A*Z) - error);   
        %grad = grad + 2*(Z(:,k)'*A*Z(:,k) - error(k))*(Z(:,k)*Z(:,k)');
        grad = grad + 2*(tanh(Z(:,k)'*A*Z(:,k)) - error(k))*(1- (tanh(Z(:,k)'*A*Z(:,k)))^2)*(Z(:,k)*Z(:,k)');
    end    
    grad_F = grad;    
    A = A - step_size*(grad_F);
        
    A = projectionIntoSPD(A); % projection into positive definite space %
    
    dist = zeros(numVec,1);
    %dist = diag(Z'*A*Z) - error;
    dist = (tanh(diag(Z'*A*Z)) - error); 
 
    optValGradientDescent = (1/numVec)*sum(dist.^2);
        
    A
    optValGradientDescent
    numIteration
end
%% Best pair of signal %%
[maxDist, index] = max(dist);
mostDistinctSignal = dist(index);

end
