function [error] = computedError(X, Z)
% Compute error for different funtions %

%Z = Y - X;
[dim,numVec] = size(Z);
%error = zeros(numVec,1);
% N =[ 1.4478    0.3629        0.6896;
%     0       1.9469          0.7065;
%     0       0           1.2857];
%  N =[ 0     0           1.3994;
%     0       1.9469       0.4944;
%     1.6896    0.7065    0.6896];
 N =[ 1.4478    0.7617    0.3994;
    0.3629    1.9469    0.4944;
    0.6896    0.7065    1.2857];

% N =  [4     2     3
%      2     8     6
%      3     6    11];
% N =  [2     0     0
%      2     4     0
%      3     6    5];
%N = [0,0,2;0,3,0;5,0,0];
%N = [ 1.4478,0,0;0,1.9469,0;0,0,1.2857];
%N = rand(3,3);
magZ = zeros(numVec, 1);
sigma = 1.15;
C1 = zeros(numVec,1);
errorType = 'Mahalanobis';

magZ = norm2Signal( Z );

switch(errorType)
    case('Mahalanobis')
        %% C(x,y) = Z'*Ae*Z and d(x,y) = Z'*M*Z and check if Ae = M %%        
        matfileAe = matfile('Ae.mat');
        Ae = matfileAe.Ae;
        C1 = diag(Z'*N*Z);
        error = C1;
    case('Weber')
        for k=1:numVec
            if(magZ(k) >= 0.6*norm(X(:,k),2))
                C1(k) = 1;
            else C1(k) = 0;
            end
        end
        error = C1;
    case('Gaussian')
        %% E(x,y) using Gaussian, hence C(x,y)=1-E(x,y), E(x,y) = exp(-norm(|x-y|,2)/2*sigma^2) %%
        %C1 = magZ.^2;
        %C1 =  exp(-(2*(sigma^2))./C1 + 1);
        C1 = diag(Z'*N*Z);
        C1 = exp(-C1/(2*(sigma^2)));
        error = C1;
    case('Entropy') % -plog(p)
        C1 = rand(numVec,1);
        C1 = -C1.*log(C1);
        error = C1;
    case('Sigmoid')
        %% C(x,y) = sigmoid function 1/(1+exp(-p))%%
        C1 = rand(numVec,1);
        C1 = 1./(1+exp(-C1));
        error = C1;
    %case('tanh&Weber')
        %% distance func by tanh(ZT*M*Z) and C by Weber's law %%
        
    otherwise
        fprintf('Invalid Input\n');
end

end