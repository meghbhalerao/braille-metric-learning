%%%%% Error %%%%%%
clear all;clc;
N=10;
error = randi([0,N]);
%normalizedError = 1-(error/N);
normalizedError = 0.8;

%% Input %% 
X = rand(3,100);
Y = rand(3,100);
Z = Y-X;
error = ones(100,1);
%%%%%% Input %%%%%%

%distance = (X)'*A*(X);
cvx_begin sdp
    variable A(3,3);    
    %maximize(distance)
    minimize (norm((diag((X-Y)'*A*(X-Y)) - error),2))
    subject to 
        A >= 0
       %A == semidefinite(3);
       %(X'*A*X)>=1
      % (X-Y)'*A*(X-Y)>=1
cvx_end