function B = generateSPDmatrix(n)
% Generate a dense n x n symmetric, positive definite matrix

B = rand(n,n); % generate a random n x n matrix

% construct a symmetric matrix using either
B = 0.5*(B+B'); 
%B = B*B';
% The first is significantly faster: O(n^2) compared to O(n^3)

% since B(i,j) < 1 by construction and a symmetric diagonally dominant matrix
%   is symmetric positive definite, which can be ensured by adding nI
B = B + n*eye(n);

end 