clc;clear;close all; clf;
load('brailletouch.mat');
load('FVrow.mat');
X = Xrow;
digits(3);
A = A*0.001;
C = (A + A')*0.5;
D = 1 - C;
[r , ~] = size(X');
pdiff = zeros((r)*(r+1)*0.5,6);
X = X';
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
for i = 1:26
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
pnorm = zeros(351,1);
p = 1;
for n = 1:351
    pnorm(n,:) = norm(pdiffp(n,:),2);
end
c = 1;
Dest = zeros(26,26);
for j = 1:26
    Dest(j:end,j) = pnorm(c:c+26-j,:) ;
    c = c + 27 - j;
end
Dest = (Dest + Dest');
maxx = max(max(Dest));
ascend_matrix = zeros(26,26);
ascend_matrix_index = zeros(26,26);
ascend_matrixD = zeros(26,26);
ascend_matrix_indexD = zeros(26,26);
for i = 1:26
    [ascend_matrix(i,:),ascend_matrix_index(i,:)] = sort(Dest(i,:)) ;
end
for i = 1:26
    [ascend_matrixD(i,:),ascend_matrix_indexD(i,:)] = sort(D(i,:)) ;
end
X_axis = 1:26;
for i = 1:1
    m = ascend_matrix_index(i,:);
    c = cellstr(char(m'+64));
    plot(X_axis,ascend_matrix(i,:)/maxx,'b.');
    hold on;
    text(X_axis,(ascend_matrix(i,:)+0.08)/maxx,c);
end

for i = 1:1
    m = ascend_matrix_indexD(i,:);
    c = cellstr(char(m'+64));
    plot(X_axis,ascend_matrixD(i,:),'r.');
    hold on;
    text(X_axis,ascend_matrixD(i,:)+0.08,c);
end
xlabel('Alphabets in ascending order of their distance from A','FontSize',17,'FontWeight','bold');
ylabel('The normalized numerical values of distances','FontSize',17,'FontWeight','bold');
legend('Estimated Distances','Human Response Distances');
C1 = X(1,:);
C2 = zeros(5,6);
C2(1,:) = X(2,:);
C2(2,:) = X(3,:);
C2(3,:) = X(5,:);
C2(4,:) = X(9,:);
C2(5,:) = X(11,:);
C3 = zeros(9,6);
C3(1,:) = X(4,:);
C3(2,:) =X(6,:) ;
C3(3,:) =X(8,:) ;
C3(4,:) =X(10,:) ;
C3(5,:) =X(12,:) ;
C3(6,:) = X(13,:);
C3(7,:) = X(15,:);
C3(8,:) = X(19,:);
C3(9,:) =X(21,:) ;
C4 = zeros(9,6);
C4(1,:) = X(7,:);
C4(2,:) = X(14,:);
C4(3,:) = X(16,:);
C4(4,:) = X(18,:);
C4(5,:) = X(20,:);
C4(6,:) = X(22,:);
C4(7,:) = X(23,:);
C4(8,:) = X(24,:);
C4(9,:) = X(26,:);
C5 = zeros(2,6);
C5(1,:) = X(17,:);
C5(2,:) = X(25,:);
Cf = zeros(5,5);
Cf(1,1) = sum(sum(pdist2(C1,C1,'mahalanobis',inv(M))));
Cf(1,2) = sum(sum(pdist2(C1,C2,'mahalanobis',inv(M))));
Cf(1,3) = sum(sum(pdist2(C1,C3,'mahalanobis',inv(M))));
Cf(1,4) = sum(sum(pdist2(C1,C4,'mahalanobis',inv(M))));
Cf(1,5) = sum(sum(pdist2(C1,C5,'mahalanobis',inv(M))));
Cf(2,1) = sum(sum(pdist2(C2,C1,'mahalanobis',inv(M))));
Cf(2,2) = sum(sum(pdist2(C2,C2,'mahalanobis',inv(M))));
Cf(2,3) = sum(sum(pdist2(C2,C3,'mahalanobis',inv(M))));
Cf(2,4) = sum(sum(pdist2(C2,C4,'mahalanobis',inv(M))));
Cf(2,5) = sum(sum(pdist2(C2,C5,'mahalanobis',inv(M))));
Cf(3,1) = sum(sum(pdist2(C3,C1,'mahalanobis',inv(M))));
Cf(3,2) = sum(sum(pdist2(C3,C2,'mahalanobis',inv(M))));
Cf(3,3) = sum(sum(pdist2(C3,C3,'mahalanobis',inv(M))));
Cf(3,4) = sum(sum(pdist2(C3,C4,'mahalanobis',inv(M))));
Cf(3,5) = sum(sum(pdist2(C3,C5,'mahalanobis',inv(M))));
Cf(4,1) = sum(sum(pdist2(C4,C1,'mahalanobis',inv(M))));
Cf(4,2) = sum(sum(pdist2(C4,C2,'mahalanobis',inv(M))));
Cf(4,3) = sum(sum(pdist2(C4,C3,'mahalanobis',inv(M))));
Cf(4,4) = sum(sum(pdist2(C4,C4,'mahalanobis',inv(M))));
Cf(4,5) = sum(sum(pdist2(C4,C5,'mahalanobis',inv(M))));
Cf(5,1) = sum(sum(pdist2(C5,C1,'mahalanobis',inv(M))));
Cf(5,2) = sum(sum(pdist2(C5,C2,'mahalanobis',inv(M))));
Cf(5,3) = sum(sum(pdist2(C5,C3,'mahalanobis',inv(M))));
Cf(5,4) = sum(sum(pdist2(C5,C4,'mahalanobis',inv(M))));
Cf(5,5) = sum(sum(pdist2(C5,C5,'mahalanobis',inv(M))));

Cfn = zeros(5,5);
Cfn(1,1) = sum(sum(pdist2(C1,C1,'mahalanobis',inv(M))))/1;
Cfn(1,2) = sum(sum(pdist2(C1,C2,'mahalanobis',inv(M))))/5;
Cfn(1,3) = sum(sum(pdist2(C1,C3,'mahalanobis',inv(M))))/9;
Cfn(1,4) = sum(sum(pdist2(C1,C4,'mahalanobis',inv(M))))/9;
Cfn(1,5) = sum(sum(pdist2(C1,C5,'mahalanobis',inv(M))))/2;
Cfn(2,1) = sum(sum(pdist2(C2,C1,'mahalanobis',inv(M))))/5;
Cfn(2,2) = sum(sum(pdist2(C2,C2,'mahalanobis',inv(M))))/25;
Cfn(2,3) = sum(sum(pdist2(C2,C3,'mahalanobis',inv(M))))/45;
Cfn(2,4) = sum(sum(pdist2(C2,C4,'mahalanobis',inv(M))))/45;
Cfn(2,5) = sum(sum(pdist2(C2,C5,'mahalanobis',inv(M))))/10;
Cfn(3,1) = sum(sum(pdist2(C3,C1,'mahalanobis',inv(M))))/9;
Cfn(3,2) = sum(sum(pdist2(C3,C2,'mahalanobis',inv(M))))/45;
Cfn(3,3) = sum(sum(pdist2(C3,C3,'mahalanobis',inv(M))))/81;
Cfn(3,4) = sum(sum(pdist2(C3,C4,'mahalanobis',inv(M))))/81;
Cfn(3,5) = sum(sum(pdist2(C3,C5,'mahalanobis',inv(M))))/18;
Cfn(4,1) = sum(sum(pdist2(C4,C1,'mahalanobis',inv(M))))/9;
Cfn(4,2) = sum(sum(pdist2(C4,C2,'mahalanobis',inv(M))))/45;
Cfn(4,3) = sum(sum(pdist2(C4,C3,'mahalanobis',inv(M))))/81;
Cfn(4,4) = sum(sum(pdist2(C4,C4,'mahalanobis',inv(M))))/81;
Cfn(4,5) = sum(sum(pdist2(C4,C5,'mahalanobis',inv(M))))/18;
Cfn(5,1) = sum(sum(pdist2(C5,C1,'mahalanobis',inv(M))))/2;
Cfn(5,2) = sum(sum(pdist2(C5,C2,'mahalanobis',inv(M))))/10;
Cfn(5,3) = sum(sum(pdist2(C5,C3,'mahalanobis',inv(M))))/18;
Cfn(5,4) = sum(sum(pdist2(C5,C4,'mahalanobis',inv(M))))/18;
Cfn(5,5) = sum(sum(pdist2(C5,C5,'mahalanobis',inv(M))))/4;