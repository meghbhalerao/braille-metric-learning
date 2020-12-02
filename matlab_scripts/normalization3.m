function [ A_norm ] = normalization3( A )
%UNTITLED2 Summary of this function goes here
% %  A is a matrix of signals with row = no of observation, col = dim.
% Every signal is normalized to norm 1
[m,n] = size(A);
% A_norm = A/norm(A);
for i=1:m
  A_norm(i,:) = A(i,:)/norm(A(i,:));
end 
end

