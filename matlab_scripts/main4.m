load('X.mat');
load('textureMatrix69.mat');
load('w.mat');
ro = 1;
Ipc=0;
x =0;
xd =0;
y  = 0;
flag = 0;
a = [];
dist = [];
epsilonijks = [];
epsilontill = 0;
Xn = normalization3(Xn);
disMatrix = 30 - textureMatrix69;
disMatrixNorm = disMatrix/30;
n = 1000;
gg =0;
Ips = [];
Ipprev = 0;
Ars = zeros(6*n,6);
triplet =  gkt(disMatrix, n);
normar = [];
sumArtill = zeros(6,6);
triplet_sorted = [];
%initial guess for M matrix
M = zeros(6,6);
M = (M+M')/2;
[V,L] = eig(M);
L = max(L, 0);
M = V*L*V';

I = eye(6,6);
alpha = 0.100;
% for i = 1:n
%     Ar = (Xn(triplet(i,1),:) - Xn(triplet(i,3),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,3),:)) - (Xn(triplet(i,1),:) - Xn(triplet(i,2),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,2),:));
%     Ars([ro:ro+5],[1:6]) = Ar;
%     ro = 6*i + 1;
%     normar(i,:) = sum(sum(Ar/norm(Ar)));
% end
% [normar,order] = sort(normar,'ascend');
% for p = 1:n
%  triplet_sorted(p,:) = triplet(order(p,:),:);
%  end
%  triplet = triplet_sorted;


ro =1;
for i = 1:n
    Ar = (Xn(triplet(i,1),:) - Xn(triplet(i,3),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,3),:)) - (Xn(triplet(i,1),:) - Xn(triplet(i,2),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,2),:));
    Ars([ro:ro+5],[1:6]) = Ar;
    ro = 6*i + 1;
    epsilonijk = disMatrixNorm(triplet(i,1),triplet(i,3)) - disMatrixNorm(triplet(i,1),triplet(i,2));
    epsilonijks(i,:) = epsilonijk ;
    Ipcheck = sum(sum(Ar.*M));
    if Ipcheck<1-epsilonijk
        sumArtill = sumArtill + Ar;
        epsilontill = epsilontill + epsilonijk;
        Ip = sum(sum(sumArtill.*M));
        y = y+1;
        while(Ip<i-epsilontill)
            M = M + (i-Ip)*(sumArtill)/(norm(sumArtill,'fro'));
            Ip = sum(sum(sumArtill.*M));
        end
        if sum(sum(M.*Ar))>1-epsilonijk
            x = x+1;
        end
    end
    M = M - alpha*(I)/norm(I,'fro');
    a(i,:)=trace(M);
    M = (M + M')/2;
    [V,L] = eig(M);
    L = max(L, 0);
    M = V*L*V';
    if sum(sum(M.*Ar))>1-epsilonijk
        xd = xd+1;
    end
    M = M/norm(M,'fro');
end
if sum(sum(sumArtill.*M)>y - epsilontill)
    flag = 1;
end

X_projected = ((((M)^(0.5))'*Xn')');
count = 0;
count1 = 0;
count2 = 0;
sar = [];
[r,c] = size(triplet);
for i = 1:r
    
    if (X_projected(triplet(i,1),:) - X_projected(triplet(i,3),:)) * (X_projected(triplet(i,1),:) - X_projected(triplet(i,3),:))' - (X_projected(triplet(i,1),:) - X_projected(triplet(i,2),:)) * (X_projected(triplet(i,1),:) - X_projected(triplet(i,2),:))'>1 - epsilonijks(i,:)
        count = count + 1;
    end
    if(norm(Xn(triplet(i,3),:)-Xn(triplet(i,1),:)) > norm(Xn(triplet(i,2),:)-Xn(triplet(i,1),:)) )
        count1 = count1 + 1;
    end
    
end
M,y;x;xd;

