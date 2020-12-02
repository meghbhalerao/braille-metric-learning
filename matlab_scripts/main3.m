load('X.mat');
load('textureMatrix69.mat');
load('w.mat');
ro = 1;
Ipc=0;
x =0;
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
%code to sort the triplets in order of difference of distance between similar and dissimilar signals
%  for h = 1:n
%  dist(h,:) = (Xn(triplet(h,1),:) - Xn(triplet(h,3),:)) * (Xn(triplet(h,1),:) - Xn(triplet(h,3),:))' - (Xn(triplet(h,1),:) - Xn(triplet(h,2),:)) * (Xn(triplet(h,1),:) - Xn(triplet(h,2),:))'
%  end
%  [dist,order] = sort(dist,'ascend');
triplet_sorted = [];
%  for p = 1:n
%  triplet_sorted(p,:) = triplet(order(p,:),:);
%  end
%  triplet = triplet_sorted;

%initial guess for M matrix
M = 10*rand(6,6);
M = (M+M')/2;
[V,L] = eig(M);
    L = max(L, 0);
    M = V*L*V';

I = eye(6,6);
alpha = 0.100;
alpha2 = 0.05;% step size some appropriate fixed number
%loop for running through all the triplets
epsiloniter = 0.01;
for i = 1:n
    Ar = (Xn(triplet(i,1),:) - Xn(triplet(i,3),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,3),:)) - (Xn(triplet(i,1),:) - Xn(triplet(i,2),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,2),:));
    Ars([ro:ro+5],[1:6]) = Ar;
    ro = 6*i + 1;
    normar(i,:) = norm(Ar,'fro');
end
[normar,order] = sort(normar,'descend');
% for p = 1:n
%  triplet_sorted(p,:) = triplet(order(p,:),:);
%  end
%  triplet = triplet_sorted;
ro =1;
% while(flag~=1)
%  triplet =  gkt(disMatrix, n);
% gg = gg+1;
for i = 1:n
    M_prev = M;
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
    M = M/norm(M,'fro');

end
if sum(sum(sumArtill.*M)>y - epsilontill)
    flag = 1;
end
% M  = M + 0.01*eye(6,6);
%  end

X_projected = ((((M)^(0.5))'*Xn')');
count = 0;
count1 = 0;
count2 = 0;
sar = [];
[r,c] = size(triplet);
for i = 1:r
    
    if(sum(sum(((X_projected(triplet(i,1),:) - X_projected(triplet(i,3),:))' * (X_projected(triplet(i,1),:) - X_projected(triplet(i,3),:)) - (X_projected(triplet(i,1),:) - X_projected(triplet(i,2),:))' * (X_projected(triplet(i,1),:) - X_projected(triplet(i,2),:))).*M ))>1-epsilonijks(i,:))
        count = count + 1;
    end
    if(norm(Xn(triplet(i,3),:)-Xn(triplet(i,1),:)) > norm(Xn(triplet(i,2),:)-Xn(triplet(i,1),:)) )
        count1 = count1 + 1;
    end
    if(norm(X_projected(triplet(i,3),:)-X_projected(triplet(i,1),:)) > norm(X_projected(triplet(i,2),:)-X_projected(triplet(i,1),:))+1-epsilonijks(i,:))
        count2 = count2 + 1;
    end
    
end

M
eigenval =eig(M)
y,x

% cc = 0;
% for l = 1:6:6*n
% cc = cc + 1;
% sar(cc,:) = sum(sum(Ars([l:l+5],[1:6]).*M));
% if(sar(cc,:)>1 - epsilonijks(cc,:))
%  s= s+1
% endAr
%
% end
%
