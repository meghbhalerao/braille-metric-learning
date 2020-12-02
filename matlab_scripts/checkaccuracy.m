X_projected = (M*Xn')';
count = 0;
sar = [];
[r,c] = size(triplet);
for i = 1:r
    if(norm(X_projected(triplet(i,3),:)-X_projected(triplet(i,1),:)) > norm(X_projected(triplet(i,2),:)-X_projected(triplet(i,1),:)) )
        count = count + 1;
    end
    
end
cc = 0;
for l = 1:6:6*n
cc = cc + 1;
sar(cc,:) = sum(sum(Ars([l:l+5],[1:6]).*M));
end
