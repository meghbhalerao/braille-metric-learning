Xn = rand(100,6);
Xn = normalization3(Xn);
Mi = generateSPDmatrix(6);
C = inv(Mi);
dis = pdist2(Xn,Xn,'mah',C);
n = 1000;
triplet = gkt(dis,1000);
dis = dis/max(max(dis));
M = eye(6,6);
I = eye(6,6);
Mo= M;
dis = dis/max(max(dis));
iter1 = 0;
iter2 = 0;
alpha = 0.001;
sumAr = 0;
esum = 0;
ro = 1;
err2 = 100000000000;
err2s = [];
c2 =0;
err1 = 100000000000;
c1 = 0;
err1s = [];
x = 0; y = 0;
for i = 1:n
    Ar = (Xn(triplet(i,1),:) - Xn(triplet(i,3),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,3),:)) - (Xn(triplet(i,1),:) - Xn(triplet(i,2),:))' * (Xn(triplet(i,1),:) - Xn(triplet(i,2),:));
    sumAr = sumAr + Ar;
    epsilonijk = dis(triplet(i,1),triplet(i,3)) - dis(triplet(i,1),triplet(i,2));
    epsilonijks(i,:) = epsilonijk ;
    Ars([ro:ro+5],[1:6]) = Ar;
    ro = 6*i + 1;
    esum = esum + epsilonijk;
end
while (iter2<500 && err2>=0.1)
    while(iter1<1000 && err1>=0.1)
        if(sum(sum(sumAr.*M)) > n - esum)
            M = M;
            x = x+1;
        else
            M = M + (sumAr/norm(sumAr,'fro'))*(n - esum -  sum(sum(sumAr.*M)))
            y= y+1;
        end
        M = (M + M')/2;
        [V,L] = eig(M);
        L = max(L, 0);
        M = V*L*V';
        iter1 = iter1 + 1;
        err1 = abs(sum(sum(sumAr.*M)) - n + esum);
        err1s(c1+1,:) = err1;
        c1 = c1 + 1;
        M= M /norm(M,'fro');
    end
    if err1<0.01
        M = M - alpha*I;
    else
        M = Mo - alpha*I;
        Mo = M;
    end
    iter2 = iter2 + 1;
    err2 = abs(sum(sum(sumAr.*M)) - n + esum);    
    iter1 = 0;
end













