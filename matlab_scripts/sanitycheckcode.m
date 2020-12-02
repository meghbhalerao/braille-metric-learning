load('X.mat');
Xn = normalization3(Xn);
Mi = generateSPDmatrix(6);
C = inv(Mi);
discheck = pdist2(Xn,Xn,'mah',C);
triplet = gkt(discheck,1000); 
ro = 1;
 Ipc=0;
 x =0;
 y  = 0;
 flag = 0;
 a = [];
 dist = [];
 epsilonijks = [];
 epsilontill = 0;
 n = 1000;
 Ips = [];
 Ipprev = 0;
 Ars = zeros(6*n,6);
 normar = [];
 sumArtill = zeros(6,6);
 M = eye(6,6);
 I = eye(6,6);
 alpha = 0.100;
 epsiloniter = 0.01;
 for i = 1:n
 Ar = (Xn(triplet(i,1),:) - Xn(triplet(i,3),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,3),:)) - (Xn(triplet(i,1),:) - Xn(triplet(i,2),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,2),:));
     Ars([ro:ro+5],[1:6]) = Ar;
     ro = 6*i + 1;
 normar(i,:) = norm(Ar,'fro');
 end
ro =1;
     discheck = discheck/max(max(discheck));

for i = 1:n
     Ar = (Xn(triplet(i,1),:) - Xn(triplet(i,3),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,3),:)) - (Xn(triplet(i,1),:) - Xn(triplet(i,2),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,2),:));
     Ars([ro:ro+5],[1:6]) = Ar;
     ro = 6*i + 1;
     epsilonijk = discheck(triplet(i,1),triplet(i,3)) - discheck(triplet(i,1),triplet(i,2));
     epsilonijks(i,:) = epsilonijk ;
     Ipcheck = sum(sum(Ar.*M));
    if Ipcheck<i-epsilonijk
     sumArtill = sumArtill + Ar;
     epsilontill = epsilontill + epsilonijk;
     Ip = sum(sum(sumArtill.*M));
        y = y+1;
        while(Ip<i-epsilontill)
         
         M = M + (i-Ip)*(sumArtill)/(norm(sumArtill,'fro'));
         Ip = sum(sum(sumArtill.*M));
        end
       if sum(sum(M.*Ar))>Ipcheck
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
sumArtill =0;


X_projected = ((((M)^(0.5))'*Xn')');
count1 = 0;
count2 = 0;
count3 = 0;
sar = [];
[r,c] = size(triplet);
for i = 1:r
    
    if(norm(Xn(triplet(i,3),:)-Xn(triplet(i,1),:)) > norm(Xn(triplet(i,2),:)-Xn(triplet(i,1),:)) )
    count1 = count1 + 1;
    end
    if(norm(X_projected(triplet(i,3),:)-X_projected(triplet(i,1),:)) > norm(X_projected(triplet(i,2),:)-X_projected(triplet(i,1),:)) + 1 - epsilonijks(i,:) )
    count2 = count2 + 1;
    end
    if(norm(X_projected(triplet(i,3),:)-X_projected(triplet(i,1),:)) > norm(X_projected(triplet(i,2),:)-X_projected(triplet(i,1),:)))
    count3 = count3 + 1;
    end
    
end     
M

x
 
 
 