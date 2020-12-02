Xn = rand(100,6);
Xn = normalization3(Xn);
Mi = generateSPDmatrix(6);
C = inv(Mi);
dis = pdist2(Xn,Xn,'mah',C);
n = 1000;
y= 0;x = 0;xd = 0;
sumArtill = zeros(6,6);
epsilontill = 0;
I = eye(6,6);
epsilonijks = [];
triplet = gkt(dis,1000);
dis = dis/max(max(dis));
M = eye(6,6);
ro = 1;
alpha = 0.05;
no = 0;
for i = 1:n
    Ar = (Xn(triplet(i,1),:) - Xn(triplet(i,3),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,3),:)) - (Xn(triplet(i,1),:) - Xn(triplet(i,2),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,2),:));
    epsilonijk = dis(triplet(i,1),triplet(i,3)) - dis(triplet(i,1),triplet(i,2));
    epsilonijks(i,:) = epsilonijk ; 
    Ars([ro:ro+5],[1:6]) = Ar;
    ro = 6*i + 1;
end

for i = 1:n
    Ar = (Xn(triplet(i,1),:) - Xn(triplet(i,3),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,3),:)) - (Xn(triplet(i,1),:) - Xn(triplet(i,2),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,2),:));
    Ars([ro:ro+5],[1:6]) = Ar;
    ro = 6*i + 1;
    epsilonijk = dis(triplet(i,1),triplet(i,3)) - dis(triplet(i,1),triplet(i,2));
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
%     M = (M + M')/2;
%     [V,L] = eig(M);
%     L = max(L, 0);
%     M = V*L*V';
    if sum(sum(M.*Ar))>1-epsilonijk
        xd = xd+1;
    end
    M = M/norm(M,'fro');
end
X_projected = ((((M)^(0.5))*Xn')');
acce = 0;
acc0 = 0;
[r,c] = size(triplet);
for i = 1:r
    if (X_projected(triplet(i,1),:) - X_projected(triplet(i,3),:)) * (X_projected(triplet(i,1),:) - X_projected(triplet(i,3),:))' - (X_projected(triplet(i,1),:) - X_projected(triplet(i,2),:)) * (X_projected(triplet(i,1),:) - X_projected(triplet(i,2),:))'> 1 - epsilonijks(i,:)
        acce = acce + 1;
    end
    if sum(sum(M.*((Xn(triplet(i,1),:) - Xn(triplet(i,3),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,3),:)) - (Xn(triplet(i,1),:) - Xn(triplet(i,2),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,2),:)))))> 0
        acc0 = acc0 + 1;
    end
end
M
eig(M)
