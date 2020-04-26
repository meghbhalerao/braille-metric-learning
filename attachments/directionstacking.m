clear;clc;
load('brailletouch.mat');
A = A/1000;
Atemp = A;

load('FVrow.mat');
rearr_A = zeros(26,26);

rearr_A([1:5],:) = [A(3,:);A(13,:);A(21,:);A(24,:);A(26,:)];
rearr_A([6:17],:) = [A(2,:);A(11,:);A(12,:);A(14,:);A(16,:);A(17,:);A(18,:);A(19,:);A(20,:);A(22,:);A(23,:);A(25,:)];
rearr_A([18:26],:) = [A(1,:);A(4,:);A(5,:);A(6,:);A(7,:);A(8,:);A(9,:);A(10,:);A(15,:)];

rearr_A2 =  zeros(26,26);
A = rearr_A;

rearr_A2(:,[1:5]) = [A(:,3) A(:,13) A(:,21) A(:,24) A(:,26)];
rearr_A2(:,[6:17]) = [A(:,2) A(:,11) A(:,12) A(:,14) A(:,16) A(:,17) A(:,18) A(:,19) A(:,20) A(:,22) A(:,23) A(:,25)];
rearr_A2(:,[18:26]) = [A(:,1) A(:,4) A(:,5) A(:,6) A(:,7) A(:,8) A(:,9) A(:,10) A(:,15)];

A = rearr_A2;

Cf3 = zeros(3,3);

Cf3(1,1) = sum(sum(A([1:5],[1:5])));
Cf3(1,2) = sum(sum(A([1:5],[6:17])));
Cf3(1,3) = sum(sum(A([1:5],[18:26])));
Cf3(2,1) = sum(sum(A([6:17],[1:5])));
Cf3(2,2) = sum(sum(A([6:17],[6:17])));
Cf3(2,3) = sum(sum(A([6:17],[18:26])));
Cf3(3,1) = sum(sum(A([18:26],[1:5])));
Cf3(3,2) = sum(sum(A([18:26],[6:17])));
Cf3(3,3) = sum(sum(A([18:26],[18:26])));


s = sum(Cf3(1,:));
Cf3(1,:) = Cf3(1,:)/s;

s = sum(Cf3(2,:));
Cf3(2,:) = Cf3(2,:)/s;

s = sum(Cf3(3,:));
Cf3(3,:) = Cf3(3,:)/s;


C = Cf3;
D = 1-C;

X = Xrow;
X = X';

Xg1 = [X(3,:);X(13,:);X(21,:);X(24,:);X(26,:)];
Xg2 = [X(2,:);X(11,:);X(12,:);X(14,:);X(16,:);X(17,:);X(18,:);X(19,:);X(20,:);X(22,:);X(23,:);X(25,:)];
Xg3 = [X(1,:);X(4,:);X(5,:);X(6,:);X(7,:);X(8,:);X(9,:);X(10,:);X(15,:)];

X_sum_unroll = zeros(6,36);
X_sum_temp = zeros(6,6);
% 1 with 1
for i = 1:5
    for j = i:5
        X_sum_temp = Xg1(i,:)' * Xg1(i,:) - Xg1(j,:)' * Xg1(j,:) + X_sum_temp;
    end
end


X_sum_unroll(1,:) = unroll(X_sum_temp)';
X_sum_temp = 0;
%1 with 2

for i = 1:5
    for j = 1:12
        X_sum_temp = (Xg1(i,:)' - Xg2(j,:)')* (Xg1(i,:) - Xg2(j,:)) + X_sum_temp;
    end
end

X_sum_unroll(2,:) = unroll(X_sum_temp)';
X_sum_temp = 0;

%1 with 3
for i = 1:5
    for j = 1:9
X_sum_temp = (Xg1(i,:)' - Xg3(j,:)')* (Xg1(i,:) - Xg3(j,:)) + X_sum_temp;    
    end
end
X_sum_unroll(3,:) = unroll(X_sum_temp)';
X_sum_temp = 0;

%2 with 2


for i = 1:12
    for j = i:12
        X_sum_temp = (Xg2(i,:)' - Xg2(j,:)')* (Xg2(i,:) - Xg2(j,:)) + X_sum_temp;
    end
end

X_sum_unroll(4,:) = unroll(X_sum_temp)';
X_sum_temp = 0;

%2 with 3


for i = 1:12
    for j = 1:9
        X_sum_temp = (Xg2(i,:)' - Xg3(j,:)')* (Xg2(i,:) - Xg3(j,:)) + X_sum_temp;
    end
end

X_sum_unroll(5,:) = unroll(X_sum_temp)';
X_sum_temp = 0;

for i = 1:9
    for j = i:9
        X_sum_temp = (Xg3(i,:)' - Xg3(j,:)')* (Xg3(i,:) - Xg3(j,:)) + X_sum_temp;
    end
end

X_sum_unroll(6,:) = unroll(X_sum_temp)';
X_sum_temp = 0;

E = [];
Cf3 = (Cf3 + Cf3')/2;
%Cf3 = Cf3/max(max(Cf3));
Df3 = 1 - Cf3;
for i = 1:3
    E = horzcat(E,Df3(i:end,i)');
end
X_sum_unroll = abs(X_sum_unroll);

cvx_begin
variable M(6,6) symmetric semidefinite
minimize( norm(X_sum_unroll*unroll(M)- E') )
cvx_end
Dest_mat = X_sum_unroll*unroll(M);
c = 1;
Dest = zeros(3,3);
for j = 1:3
    Dest(j:end,j) = Dest_mat(c:c+3-j,:) ;
    c = c + 4 - j;
end
Dests = Dest + Dest';
Dests(1,1) = Dest(1,1);
Dests(2,2) = Dest(2,2);
Dests(3,3) = Dest(3,3);