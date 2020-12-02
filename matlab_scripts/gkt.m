function [triplets] = gkt(dist_mat, num)

nData = length(dist_mat);
m = 1;
n = 0;
triplets = [];

while(n~=num)
    n = n + 1;
    temp = randperm(nData, 3);
    
    i = temp(1);
    j = temp(2);
    k = temp(3);
    dist = dist_mat(i,k) - dist_mat(i,j);                 
    if (dist_mat(i,k)>dist_mat(i,j))
        triplets(n, :) = [i, j, k];
        m = m + 1;
%     else
%         triplets(n, :) = [i, k, j];
%         m = m + 1;
    end
    
   triplets = unique(triplets, 'rows');
   n = size(triplets, 1);
end

end

