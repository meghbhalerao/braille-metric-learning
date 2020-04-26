clear all;clc;
%% Initialization %%%
dim = 3;
numVec = 10;
numIteration = 0;
step_size = 0.001;
%A = eye(dim);
A = generateSPDmatrix(dim);
%Ae = generateSPDmatrix(dim);
Z = zeros(dim,numVec);
error = zeros(numVec,1);
%C = zeros(numVec,1);

%% Input %%
% for k = 1:numVec
%     X=rand(dim,1);
%     Y=rand(dim,1);
%     Z(:,k)=X-Y;
%     error(k) = rand(1);
% end
matfileZ = matfile('Input.mat');
Z = matfileZ.Z;
matfileError = matfile('error.mat');
error = matfileError.error;

% matfileAe = matfile('Ae.mat');
% Ae = matfileAe.Ae;
% A = Ae;
%matfileA = matfile('A.mat');
%A = matfileA.A;

%% Error Computation %%
% for k = 1:numVec
%    C(k) = Z(:,k)'*Ae*Z(:,k);
% end 
% error = C;

%%%%%%%  different error function %%%%%%

error = -error.*log(error);
% for k = 1:numVec
%     error(k) = 1/(1+exp(-error(k)));
%  end 
%% gradient descent step %%
optValGradientDescent = 10;
while (numIteration <= 5000 && optValGradientDescent >= 10^(-10)) % 
    numIteration = numIteration + 1;   
    grad = zeros(dim, dim);
    for k = 1:numVec
        grad = grad + 2*(Z(:,k)'*A*Z(:,k) - error(k))*(Z(:,k)*Z(:,k)');
    end    
    grad_F = grad;    
    A = A - step_size*(grad_F);
        
    A = projectionIntoSPD(A); % projection into positive definite space %
    
    dist = zeros(numVec,1);
    for k = 1:numVec
        dist(k) = (Z(:,k)'*A*Z(:,k)) - error(k);        
    end
    optValGradientDescent = sum(dist.^2);
        
    A
    optValGradientDescent
    numIteration
end

%% Best pair of signal %%
[maxDist, index] = max(dist);
mostDistinctSignal = dist(index);

%% Optimal Value %%
[optValClosedForm, M_closedForm] = closedForm( Z, error);
M_closedForm
optValClosedForm
%[cvx_optval, M_CVX] = metricUsingCVX(error, Z); % Using CVX







