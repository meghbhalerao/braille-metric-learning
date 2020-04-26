function [ mag ] = norm2Signal( Z )
%% mag is a column vector %%
[dim,numVec] = size(Z);
for k = 1:numVec
    magZ(k) = norm(Z(:,k),2);
end
mag = magZ';
end

