clc;clear;
load('brailletouch.mat');
load('FVcolumn.mat');
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

ascend_matrix = zeros(26,26);
ascend_matrix_index = zeros(26,26);
ascend_matrixD = zeros(26,26);
ascend_matrix_indexD = zeros(26,26);
maxx = max(max(Dest));
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
    plot(X_axis,ascend_matrix(i,:)/maxx,'.');
    hold on;
    text(X_axis,(ascend_matrix(i,:)+0.08)/maxx,c);
end

for i = 1:1
    m = ascend_matrix_indexD(i,:);
    c = cellstr(char(m'+64));
    plot(X_axis,ascend_matrixD(i,:),'.');
    hold on;
    text(X_axis,(ascend_matrixD(i,:)+0.08),c);
end
xlabel('Alphabets in ascending order of their distance from A','FontSize',17,'FontWeight','bold');
ylabel('The normalized numerical values of distances','FontSize',17,'FontWeight','bold');
legend('Estimated Distances','Human Response Distances');







