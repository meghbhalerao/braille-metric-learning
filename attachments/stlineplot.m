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
plot(unroll(D),unroll(Dest)/(max(max(Dest))),'r.');
hold on;
plot([0,1],[0,1]);
legend('Estimated Scatter Plot','Ideal Plot')

xlabel('Cell by Cell entries of human response distance matrix','FontSize',17,'FontWeight','bold');
ylabel('Cell by Cell entries of estimated distance matrix','FontSize',17,'FontWeight','bold');

