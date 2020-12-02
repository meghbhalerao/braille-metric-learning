 load('X.mat');
 load('textureMatrix69.mat');
 load('w.mat');
 ro = 1;
 dist = [];
 epsilonijks = [];
 Xn = normalization3(Xn);
 disMatrix = 30 - textureMatrix69;
 disMatrixNorm = disMatrix/30;
 n = 100;
 Ips = [];
 Ipprev = 0;
 Ars = zeros(6*n,6);
 triplet =  gkt(disMatrix, n);
 %code to sort the triplets in order of difference of distance between similar and dissimilar signals
 for h = 1:n
 dist(h,:) = (Xn(triplet(h,1),:) - Xn(triplet(h,3),:)) * (Xn(triplet(h,1),:) - Xn(triplet(h,3),:))' - (Xn(triplet(h,1),:) - Xn(triplet(h,2),:)) * (Xn(triplet(h,1),:) - Xn(triplet(h,2),:))'; 
 end
 [dist,order] = sort(dist,'descend');
 triplet_sorted = [];
 for p = 1:n
 triplet_sorted(p,:) = triplet(order(p,:),:);
 end
 triplet = triplet_sorted;
 sumAr = zeros(6,6);
 %initial guess for M matrix
 M = eye(6,6);
 I = eye(6,6);
 alpha = 0.25;
 alpha2 = 0.05;% step size some appropriate fixed number 
 %loop for running through all the triplets
 epsiloniter = 0.01;
 for i = 1:n
     % projection onto the second constrain space
     %um = unroll(M)
     Mi = M;
     Ar = (Xn(triplet(i,1),:) - Xn(triplet(i,3),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,3),:)) - (Xn(triplet(i,1),:) - Xn(triplet(i,2),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,2),:));
     Ars([ro:ro+5],[1:6]) = Ar;
     ro = 6*i + 1;
     Ip = sum(sum(Ar.*M));
     Ips(i,:) = Ip;
     epsilonijk = disMatrixNorm(triplet(i,1),triplet(i,3)) - disMatrixNorm(triplet(i,1),triplet(i,2));
     epsilonijks(i,:) = epsilonijk ;
     sumAr= sumAr + Ar;
    end
   sumepsilon = sum(epsilonijks);
   Ip = sum(sum(Ar.*M));
for maxiter = 1:100
   while(Ip<n-sumepsilon)
M = M + (n - Ip )*sumAr/norm(sumAr)
 M = (M + M')/2;
    [V,L] = eig(M);
    L = max(L, 0);
    M = V*L*V';
end
M = M - alpha*I;
end  
   
X_projected = ((M^(0.5))'*Xn')';
count = 0;
sar = [];
[r,c] = size(triplet);
for i = 1:r
    if(norm(X_projected(triplet(i,3),:)-X_projected(triplet(i,1),:)) > norm(X_projected(triplet(i,2),:)-X_projected(triplet(i,1),:)) )
        count = count + 1;
    end
    
end
s = 0;
cc = 0;
for l = 1:6:6*n
cc = cc + 1;
sar(cc,:) = sum(sum(Ars([l:l+5],[1:6]).*M));
if(sar(cc,:)>1 - epsilonijks(cc,:))
 s= s+1
end

end
      

 
 
 