load('X.mat');
load('textureMatrix69.mat');
load('w.mat');
Xn = normalization3(Xn);
disMatrix = 30 - textureMatrix69;
disMatrixNorm = disMatrix/30;
n = 1000;
M = eye(6,6);
I = eye(6,6);
Mi = generateSPDmatrix(6);
Ars = zeros(6*n,6);
triplet =  gkt(disMatrix, n);
newAr = zeros(n,36);
alpha = 0.05;
ro =1;
no = 0;
for i = 1:n
    Ar = (Xn(triplet(i,1),:) - Xn(triplet(i,3),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,3),:)) - (Xn(triplet(i,1),:) - Xn(triplet(i,2),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,2),:));
    newAr(i,:) = reshape(Ar,[1 36]);
    epsilonijk = disMatrixNorm(triplet(i,1),triplet(i,3)) - disMatrixNorm(triplet(i,1),triplet(i,2));
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
        no = no + 1;i;
    end
    while(Ip<1 - epmin)
        M = M + (1-Ip)*(minAr)/(norm(minAr,'fro'));
        Ip = sum(sum(minAr.*M));
    end
    sum(sum(M.*minAr))
%     [V,L] = eig(M);
%     L = max(L, 0);
%     I = eye(6,6);
%     M = V*L*V';
%     M = M/norm(M,'fro');
end

X_projected = ((((M)^(0.5))*Xn')');
acc0 = 0;
acc = 0;
acci = 0;
[r,c] = size(triplet);
for i = 1:r
    if (X_projected(triplet(i,1),:) - X_projected(triplet(i,3),:)) * (X_projected(triplet(i,1),:) - X_projected(triplet(i,3),:))' - (X_projected(triplet(i,1),:) - X_projected(triplet(i,2),:)) * (X_projected(triplet(i,1),:) - X_projected(triplet(i,2),:))'> 0
        acc0 = acc0 + 1;
    end
    if sum(sum(M.*((Xn(triplet(i,1),:) - Xn(triplet(i,3),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,3),:)) - (Xn(triplet(i,1),:) - Xn(triplet(i,2),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,2),:)))))> 0
        acc = acc + 1;
    end
    if (Xn(triplet(i,1),:) - Xn(triplet(i,3),:)) * (Xn(triplet(i,1),:) - Xn(triplet(i,3),:))' - (Xn(triplet(i,1),:) - Xn(triplet(i,2),:)) * (Xn(triplet(i,1),:) - Xn(triplet(i,2),:))'> 0
        acci = acci + 1;
    end
end
M