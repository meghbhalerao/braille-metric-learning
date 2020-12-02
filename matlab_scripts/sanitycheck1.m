Xn = rand(100,6);
Mi = generateSPDmatrix(6);
C = inv(Mi);
dis = pdist2(Xn,Xn,'mah',C);
n = 1000;
I = eye(6,6);
epsilonijks = [];
triplet = gkt(dis,1000);
Xn = normalization3(Xn);
dis = dis/max(max(dis));
M = eye(6,6);
ro = 1;
alpha = 0.05;
no = 0;
for i = 1:n
    Ar = (Xn(triplet(i,1),:) - Xn(triplet(i,3),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,3),:)) - (Xn(triplet(i,1),:) - Xn(triplet(i,2),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,2),:));
    newAr(i,:) = reshape(Ar,[1 36]);
    epsilonijk = dis(triplet(i,1),triplet(i,3)) - dis(triplet(i,1),triplet(i,2));
    epsilonijks(i,:) = epsilonijk ;
    Ars([ro:ro+5],[1:6]) = Ar;
    ro = 6*i + 1;
end
minAr = min(newAr);
minAr = reshape(minAr,[6 6]);
epmin = min(epsilonijks);
for i = 1:n
    M = M - alpha*I/(norm(I));
    Ip = sum(sum(minAr.*M));
    if(Ip<1-epmin)
        no = no + 1;i
    end
    while(Ip<1 - epmin)
        M = M + (1-Ip)*(minAr)/(norm(minAr,'fro'));
        Ip = sum(sum(minAr.*M));
    end
%     
%     [V,L] = eig(M);
%     L = max(L, 0);
%     I = eye(6,6);
%     M = V*L*V';
%     M = M/norm(M,'fro');
end
X_projected = ((((M)^(0.5))*Xn')'); 
count = 0;
count2 = 0;
[r,c] = size(triplet);
for i = 1:r
    if (X_projected(triplet(i,1),:) - X_projected(triplet(i,3),:)) * (X_projected(triplet(i,1),:) - X_projected(triplet(i,3),:))' - (X_projected(triplet(i,1),:) - X_projected(triplet(i,2),:)) * (X_projected(triplet(i,1),:) - X_projected(triplet(i,2),:))'> 1 - epsilonijks(i,:)
        count = count + 1;
    end
    if sum(sum(M.*((Xn(triplet(i,1),:) - Xn(triplet(i,3),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,3),:)) - (Xn(triplet(i,1),:) - Xn(triplet(i,2),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,2),:)))))> 0
        count2 = count2 + 1;
    end
end
M
eig(M)
