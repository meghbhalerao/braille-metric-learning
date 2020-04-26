function [ distance ] = computeDistance( Z, M )
% compute distance %
[dim, numVec] = size(Z);
D = zeros(numVec,1);
distanceType = 'tanh';
mahalonabisDistance = diag(Z'*M*Z);

switch(distanceType)
    case('Mahalanobis')
        distance = mahalonabisDistance;
    case('Exponential')
        D = mahalonabisDistance;
        distance = exp(-0.5*D);
    case('tanh')
        distance = tanh(mahalonabisDistance);
    otherwise
        fprintf('Invalid Input\n');
end
        
end


