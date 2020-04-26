function [ OptTrainingValue, M_closedForm] = closedForm( Z, error)
% Closed Form Solution from direct objective function %
[dim,numVec] = size(Z);
%OptTrainingValue = 0;
M_closedForm = zeros(dim, dim);
Z_temp = zeros(dim, dim);
H = zeros(numVec, (dim*dim));

for k=1:numVec
    Z_temp = Z(:,k)*Z(:,k)';
    yourvector = (Z_temp(:))';
    H(k,:) = yourvector;
end

%m = zeros(dim*dim,1);
m = pinv(H)*error;
M_closedForm = vec2mat(m,3);
%distClosedForm = zeros(numVec,1);

distClosedForm = diag(Z'*M_closedForm*Z) - error;
% for k = 1:numVec
%     distClosedForm(k) = (Z(:,k)'*M_closedForm*Z(:,k)) - error(k);
% end
%OptTrainingValue = diag(Z'*M_closedForm*Z);
OptTrainingValue = (1/numVec)*sum(distClosedForm.^2);
end