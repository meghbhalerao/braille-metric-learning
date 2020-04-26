clear;clc;
load('brailletouch.mat');
A = A/1000;
Atemp = A;

load('FVrow.mat');
rearr_A = zeros(26,26);

rearr_A(1,:) = A(1,:);
rearr_A([2:6],:) = [A(2,:);A(3,:);A(5,:);A(9,:);A(11,:)];
rearr_A([7:15],:) = [A(4,:);A(6,:);A(8,:);A(10,:);A(12,:);A(13,:);A(15,:);A(19,:);A(21,:)];
rearr_A([16:24],:) =[A(7,:);A(14,:);A(16,:);A(18,:);A(20,:);A(22,:);A(23,:);A(24,:);A(26,:)];
rearr_A([25:26],:) = [A(17,:);A(25,:)];

rearr_A2 =  zeros(26,26);
A = rearr_A; 
rearr_A2(:,1) = A(:,1);
rearr_A2(:,[2:6]) = [A(:,2) A(:,3) A(:,5) A(:,9) A(:,11)];
rearr_A2(:,[7:15]) = [A(:,4) A(:,6) A(:,8) A(:,10) A(:,12) A(:,13) A(:,15) A(:,19) A(:,21)];
rearr_A2(:,[16:24]) =[A(:,7) A(:,14) A(:,16) A(:,18) A(:,20) A(:,22) A(:,23) A(:,24) A(:,26)];
rearr_A2(:,[25:26]) = [A(:,17) A(:,25)];

A = rearr_A2;

Cf5 = zeros(5,5);

Cf5(1,1) = sum(sum(A([1:1],[1:1])));
Cf5(1,2) = sum(sum(A([1:1],[2:6])));
Cf5(1,3) = sum(sum(A([1:1],[7:15])));
Cf5(1,4) = sum(sum(A([1:1],[16:24])));
Cf5(1,5) = sum(sum(A([1:1],[25:26])));
Cf5(2,1) = sum(sum(A([2:6],[1:1])));
Cf5(2,2) = sum(sum(A([2:6],[2:6])));
Cf5(2,3) = sum(sum(A([2:6],[7:15])));
Cf5(2,4) = sum(sum(A([2:6],[16:24])));
Cf5(2,5) = sum(sum(A([2:6],[25:26])));
Cf5(3,1) = sum(sum(A([7:15],[1:1])));
Cf5(3,2) = sum(sum(A([7:15],[2:6])));
Cf5(3,3) = sum(sum(A([7:15],[7:15])));
Cf5(3,4) = sum(sum(A([7:15],[16:24])));
Cf5(3,5) = sum(sum(A([7:15],[25:26])));
Cf5(4,1) = sum(sum(A([16:24],[1:1])));
Cf5(4,2) = sum(sum(A([16:24],[2:6])));
Cf5(4,3) = sum(sum(A([16:24],[7:15])));
Cf5(4,4) = sum(sum(A([16:24],[16:24])));
Cf5(4,5) = sum(sum(A([16:24],[25:26])));
Cf5(5,1) = sum(sum(A([25:26],[1:1])));
Cf5(5,2) = sum(sum(A([25:26],[2:6])));
Cf5(5,3) = sum(sum(A([25:26],[7:15])));
Cf5(5,4) = sum(sum(A([25:26],[16:24])));
Cf5(5,5) = sum(sum(A([25:26],[25:26])));

s = sum(Cf5(1,:));
Cf5(1,:) = Cf5(1,:)/s;

s = sum(Cf5(2,:));
Cf5(2,:) = Cf5(2,:)/s;

s = sum(Cf5(3,:));
Cf5(3,:) = Cf5(3,:)/s;

s = sum(Cf5(4,:));
Cf5(4,:) = Cf5(4,:)/s;

s = sum(Cf5(5,:));
Cf5(5,:) = Cf5(5,:)/s;

X = Xrow;
X = X';
Xg = zeros(5,6);

D = 1 - Cf5;

Xg(1,:) = X(1,:);
Xg(2,:) = (X(2,:)+X(3,:)+X(5,:)+X(9,:)+X(11,:))/5;
Xg(3,:) =  (X(4,:)+X(6,:)+X(8,:)+X(10,:)+X(12,:)+X(13,:)+X(15,:)+X(19,:)+X(21,:))/9;
Xg(4,:) = (X(7,:)+X(14,:)+X(16,:)+X(18,:)+X(20,:)+X(22,:)+X(23,:)+X(24,:)+X(26,:))/9 ;
Xg(5,:) =(X(17,:)+X(25,:))/2 ;

X = Xg;

[r , ~] = size(X);

pdiff = zeros((r)*(r+1)*0.5,6);

c = 1;
for i = 1:r
    for j =i:r
        pdiff(c,:) = X(i,:) - X(j,:);
        c = c + 1;
    end
end

H = zeros(r*(r+1)*0.5,36);
for k = 1: c-1
    H(k,:) = unroll((pdiff(k,:)' * pdiff(k,:))')';
end

i = 1;
j = 1;
E = [];

for i = 1:5
    E = horzcat(E,D(i:end,i)');
end

cvx_begin
variable M(6,6) symmetric semidefinite
minimize( norm(H*unroll(M)-E') )
cvx_end

Xp = (M^(0.5)*X')';
pdiffp = zeros((r)*(r+1)*0.5,6);
c = 1;
for i = 1:r
    for j =i:r
        pdiffp(c,:) = Xp(i,:) - Xp(j,:);
        c = c + 1;
    end
    
end

pnorm = zeros(15,1);
p = 1;
for n = 1:15
    pnorm(n,:) = norm(pdiffp(n,:),2);
end

c = 1;
Dest = zeros(5,5);
for j = 1:5
    Dest(j:end,j) = pnorm(c:c+5-j,:) ;
    c = c + 6 - j;
end

Dest = (Dest + Dest');

Cest = 1 - Dest;
Cf5 = 0.5*(Cf5 + Cf5');
