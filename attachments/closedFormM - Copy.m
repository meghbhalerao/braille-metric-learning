function [ optValClosedForm, M_closedForm] = closedFormM( Z, error)

% Closed Form Solution %
[dim,numVec] = size(Z);
optValClosedForm = 0;
M_closedForm = zeros(dim, dim);

% M = generateSPDmatrix(dim);
% %M = eye(dim);
% distM = zeros(numVec,1);
% for k = 1:numVec
%     distM(k) = Z(:,k)'*M*Z(:,k);
% end
% sigmaC = sum(error);
% sigmaDist = sum(distM);
% scalingConst = sigmaC/sigmaDist;
% 
% M_closedForm = scalingConst*M;

alpha = zeros(numVec,1);
Am = zeros(dim,dim,numVec);
Af = zeros(dim*dim, numVec);
M_temp = zeros(dim,dim);

for k = 1:numVec
     Am(:,:,k) = Z(:,k)*Z(:,k)';
     M_temp = Am(:,:,k);
     yourvector = M_temp(:);
     Af(:,k) = yourvector;
end
alpha = null(Af);

%% 2nd step  to get M %%
b = zeros(numVec,1);
Y = zeros(dim*dim,1);
for k = 1:numVec
    b(k) = error(k) + alpha(k);
end

Y = pinv(Af')*b;
M_closedForm = vec2mat(Y,3);

distClosedForm = zeros(numVec,1);
for k = 1:numVec
    distClosedForm(k) = (Z(:,k)'*M_closedForm*Z(:,k)) - error(k);
end
optValClosedForm = sum(distClosedForm.^2);
end