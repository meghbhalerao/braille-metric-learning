clear;clc;
load('brailletouch.mat');
A = A/1000;
Atemp = A;

load('FVrow.mat');
rearr_A = zeros(26,26);

rearr_A(1,:) = A(1,:);
index_group1 = [1];
rearr_A([2:6],:) = [A(2,:);A(3,:);A(5,:);A(9,:);A(11,:)];
index_group2 = [2 3 5 9 11];
rearr_A([7:15],:) = [A(4,:);A(6,:);A(8,:);A(10,:);A(12,:);A(13,:);A(15,:);A(19,:);A(21,:)];
index_group3 = [4 6 8 10 12 13 15 19 21];
rearr_A([16:24],:) =[A(7,:);A(14,:);A(16,:);A(18,:);A(20,:);A(22,:);A(23,:);A(24,:);A(26,:)];
index_group4 = [7 14 16 18 20 22 23 24 26];
rearr_A([25:26],:) = [A(17,:);A(25,:)];
index_group5 = [17 25];

rearr_A2 =  zeros(26,26);
A = rearr_A;
rearr_A2(:,1) = A(:,1);
rearr_A2(:,[2:6]) = [A(:,2) A(:,3) A(:,5) A(:,9) A(:,11)];
rearr_A2(:,[7:15]) = [A(:,4) A(:,6) A(:,8) A(:,10) A(:,12) A(:,13) A(:,15) A(:,19) A(:,21)];
rearr_A2(:,[16:24]) =[A(:,7) A(:,14) A(:,16) A(:,18) A(:,20) A(:,22) A(:,23) A(:,24) A(:,26)];
rearr_A2(:,[25:26]) = [A(:,17) A(:,25)];

A = rearr_A2;

Cf5 = zeros(5,5);

Cf5(1,1) = sum(sum(A([1:1],[1:1])));
Cf5(1,2) = sum(sum(A([1:1],[2:6]))) ;
Cf5(1,3) = sum(sum(A([1:1],[7:15])));
Cf5(1,4) = sum(sum(A([1:1],[16:24])));
Cf5(1,5) = sum(sum(A([1:1],[25:26])));
Cf5(2,1) = sum(sum(A([2:6],[1:1])));
Cf5(2,2) = sum(sum(A([2:6],[2:6])));
Cf5(2,3) = sum(sum(A([2:6],[7:15])));
Cf5(2,4) = sum(sum(A([2:6],[16:24])));
Cf5(2,5) = sum(sum(A([2:6],[25:26])));
Cf5(3,1) = sum(sum(A([7:15],[1:1])));
Cf5(3,2) = sum(sum(A([7:15],[2:6])));
Cf5(3,3) = sum(sum(A([7:15],[7:15])));
Cf5(3,4) = sum(sum(A([7:15],[16:24])));
Cf5(3,5) = sum(sum(A([7:15],[25:26])));
Cf5(4,1) = sum(sum(A([16:24],[1:1])));
Cf5(4,2) = sum(sum(A([16:24],[2:6])));
Cf5(4,3) = sum(sum(A([16:24],[7:15])));
Cf5(4,4) = sum(sum(A([16:24],[16:24])));
Cf5(4,5) = sum(sum(A([16:24],[25:26])));
Cf5(5,1) = sum(sum(A([25:26],[1:1])));
Cf5(5,2) = sum(sum(A([25:26],[2:6])));
Cf5(5,3) = sum(sum(A([25:26],[7:15])));
Cf5(5,4) = sum(sum(A([25:26],[16:24])));

s = sum(Cf5(1,:));
Cf5(1,:) = Cf5(1,:)/s;

s = sum(Cf5(2,:));
Cf5(2,:) = Cf5(2,:)/s;

s = sum(Cf5(3,:));
Cf5(3,:) = Cf5(3,:)/s;

s = sum(Cf5(4,:));
Cf5(4,:) = Cf5(4,:)/s;

s = sum(Cf5(5,:));
Cf5(5,:) = Cf5(5,:)/s;

C = Cf5;
D = 1-C;

X = Xrow;
X = X';


Xg1 = [X(1,:)];
Xg2 = [X(2,:);X(3,:);X(5,:);X(9,:);X(11,:)];
Xg3 = [X(4,:);X(6,:);X(8,:);X(10,:);X(12,:);X(13,:);X(15,:);X(19,:);X(21,:)];
Xg4 = [X(7,:);X(14,:);X(16,:);X(18,:);X(20,:);X(22,:);X(23,:);X(24,:);X(26,:)];
Xg5 = [X(17,:);X(25,:)];


X_sum_unroll = zeros(15,36);
X_sum_temp = zeros(6,6);

%1 with 1

for i = 1:1
    for j = 1:1
        X_sum_temp = Xg1(i,:)' * Xg1(i,:) - Xg1(j,:)' * Xg1(j,:) + X_sum_temp;
    end
end


X_sum_unroll(1,:) = unroll(X_sum_temp)';
X_sum_temp = 0;
%1 with 2

for i = 1:1
    for j = 1:5
        X_sum_temp = (Xg1(i,:)' - Xg2(j,:)')* (Xg1(i,:) - Xg2(j,:)) + X_sum_temp;
    end
end

X_sum_unroll(2,:) = unroll(X_sum_temp)'/15;
X_sum_temp = 0;

%1 with 3
for i = 1:1
    for j = 1:9
X_sum_temp = (Xg1(i,:)' - Xg3(j,:)')* (Xg1(i,:) - Xg3(j,:)) + X_sum_temp;    
    end
end

X_sum_unroll(3,:) = unroll(X_sum_temp)'/9;
X_sum_temp = 0;
%1 with 4


for i = 1:1
    for j = 1:9
    X_sum_temp = (Xg1(i,:)' - Xg4(j,:)')* (Xg1(i,:) - Xg4(j,:)) + X_sum_temp;
    end
end

X_sum_unroll(4,:) = unroll(X_sum_temp)'/9;
X_sum_temp = 0;

%1 with 5

for i = 1:1
    for j = 1:2
        X_sum_temp = (Xg1(i,:)' - Xg5(j,:)')* (Xg1(i,:) - Xg5(j,:)) + X_sum_temp;
    end
end

X_sum_unroll(5,:) = unroll(X_sum_temp)'/2;
X_sum_temp = 0;
%2 with 2


for i = 1:5
    for j = 1:5
        X_sum_temp = (Xg2(i,:)' - Xg2(j,:)')* (Xg2(i,:) - Xg2(j,:)) + X_sum_temp;
    end
end

X_sum_unroll(6,:) = unroll(X_sum_temp)'/25;
X_sum_temp = 0;

%2 with 3


for i = 1:5
    for j = 1:9
        X_sum_temp = (Xg2(i,:)' - Xg3(j,:)')* (Xg2(i,:) - Xg3(j,:)) + X_sum_temp;
    end
end

X_sum_unroll(7,:) = unroll(X_sum_temp)'/45;
X_sum_temp = 0;

%2 with 4

for i = 1:5
    for j = 1:9
        X_sum_temp = (Xg2(i,:)' - Xg4(j,:)')* (Xg2(i,:) - Xg4(j,:)) + X_sum_temp;
    end
end

X_sum_unroll(8,:) = unroll(X_sum_temp)'/45;
X_sum_temp = 0;

%2 with 5
for i = 1:5
    for j = 1:2
        X_sum_temp = (Xg2(i,:)' - Xg5(j,:)')* (Xg2(i,:) - Xg5(j,:)) + X_sum_temp;
    end
end


X_sum_unroll(9,:) = unroll(X_sum_temp)'/10;
X_sum_temp = 0;

%3 with 3

for i = 1:9
    for j = i:9
        X_sum_temp = (Xg3(i,:)' - Xg3(j,:)')* (Xg3(i,:) - Xg3(j,:)) + X_sum_temp;
    end
end

X_sum_unroll(10,:) = unroll(X_sum_temp)'/81;
X_sum_temp = 0;

%3 with 4

for i = 1:9
    for j = 1:9
        X_sum_temp = (Xg3(i,:)' - Xg4(j,:)')* (Xg3(i,:) - Xg4(j,:)) + X_sum_temp;
    end
end

X_sum_unroll(11,:) = unroll(X_sum_temp)'/81;
X_sum_temp = 0;

%3 with 5

for i = 1:9
    for j = 1:2
        X_sum_temp = (Xg3(i,:)' - Xg5(j,:)')* (Xg3(i,:) - Xg5(j,:)) + X_sum_temp;
    end
end

X_sum_unroll(12,:) = unroll(X_sum_temp)'/18;
X_sum_temp = 0;

%4 with 4

for i = 1:9
    for j = i:9
        X_sum_temp = (Xg4(i,:)' - Xg4(j,:)')* (Xg4(i,:) - Xg4(j,:)) + X_sum_temp;
    end
end

X_sum_unroll(13,:) = unroll(X_sum_temp)'/81;
X_sum_temp = 0;



%4 with 5

for i = 1:9
    for j = 1:2
        X_sum_temp = (Xg4(i,:)' - Xg5(j,:)')* (Xg4(i,:) - Xg5(j,:)) + X_sum_temp;
    end
end

X_sum_unroll(14,:) = unroll(X_sum_temp)'/18;

X_sum_temp = 0;

%5 with 5

for i = 1:2
    for j = i:2
        X_sum_temp = (Xg5(i,:)' - Xg5(j,:)')* (Xg5(i,:) - Xg5(j,:)) + X_sum_temp;
    end
end

X_sum_unroll(15,:) = unroll(X_sum_temp)'/4;
% X_sum_temp = 0;
E = [];
Co = Cf5;
Cf5 = (Cf5 + Cf5')/2;
%Cf5 = Cf5/max(max(Cf5));
Df5 = 1 - Cf5;
for i = 1:5
    E = horzcat(E,Df5(i:end,i)');
end
%X_sum_unroll = abs(X_sum_unroll);
%E = unroll(Df5);

cvx_begin
variable M(6,6) symmetric semidefinite
minimize( norm(X_sum_unroll*unroll(M)- E') )
cvx_end

Dest_mat = X_sum_unroll*unroll(M);
c = 1;
Dest = zeros(5,5);
for j = 1:5
    Dest(j:end,j) = Dest_mat(c:c+5-j,:) ;
    c = c + 6 - j;
end
Dests = Dest + Dest';


Dests(2,2) = Dest(2,2);
Dests(3,3) = Dest(3,3);
Dests(4,4) = Dest(4,4);
Dests(5,5) = Dest(5,5);
M = M/norm(M);

% accuracy check

T = [];
t_count = 1;
%code to find all possible triplets from the grouped up alphabets
%for 2 raised dots
for i = 1:5
    for j = 1:5
        if i-j        
        %writing loops for taking the alphabets except for that particular
        %group
        %for group 1
        for k = 1:1
            T(t_count,1) = index_group2(1,i);
            T(t_count,2) = index_group2(1,j);
            T(t_count,3) = index_group1(1,k);
            t_count = t_count + 1;
        end
        %for group 3
         for k = 1:9
            T(t_count,1) = index_group2(1,i);
            T(t_count,2) = index_group2(1,j);
            T(t_count,3) = index_group3(1,k);
            t_count = t_count + 1;
         end
         %for group 4
         for k = 1:9
            T(t_count,1) = index_group2(1,i);
            T(t_count,2) = index_group2(1,j);
            T(t_count,3) = index_group4(1,k);
            t_count = t_count + 1;
         end
         %for group 5
         for k = 1:2
             T(t_count,1) = index_group2(1,i);
             T(t_count,2) = index_group2(1,j);
             T(t_count,3) = index_group5(1,k);
             t_count = t_count + 1;
        end
       end  
    end
end


for i = 1:9
    for j = 1:9
        if i-j        
        %writing loops for taking the alphabets except for that particular
        %group
        %for group 1
        for k = 1:1
            T(t_count,1) = index_group3(1,i);
            T(t_count,2) = index_group3(1,j);
            T(t_count,3) = index_group1(1,k);
            t_count = t_count + 1;
        end
        %for group 2
         for k = 1:5
            T(t_count,1) = index_group3(1,i);
            T(t_count,2) = index_group3(1,j);
            T(t_count,3) = index_group2(1,k);
            t_count = t_count + 1;
         end
         %for group 4
         for k = 1:9
            T(t_count,1) = index_group3(1,i);
            T(t_count,2) = index_group3(1,j);
            T(t_count,3) = index_group4(1,k);
            t_count = t_count + 1;
         end
         %for group 5
         for k = 1:2
             T(t_count,1) = index_group3(1,i);
             T(t_count,2) = index_group3(1,j);
             T(t_count,3) = index_group5(1,k);
             t_count = t_count + 1;
        end
       end  
    end
end




for i = 1:9
    for j = 1:9
        if i-j        
        %writing loops for taking the alphabets except for that particular
        %group
        %for group 1
        for k = 1:1
            T(t_count,1) = index_group3(1,i);
            T(t_count,2) = index_group3(1,j);
            T(t_count,3) = index_group1(1,k);
            t_count = t_count + 1;
        end
        %for group 2
         for k = 1:5
            T(t_count,1) = index_group3(1,i);
            T(t_count,2) = index_group3(1,j);
            T(t_count,3) = index_group2(1,k);
            t_count = t_count + 1;
         end
         %for group 4
         for k = 1:9
            T(t_count,1) = index_group3(1,i);
            T(t_count,2) = index_group3(1,j);
            T(t_count,3) = index_group4(1,k);
            t_count = t_count + 1;
         end
         %for group 5
         for k = 1:2
             T(t_count,1) = index_group3(1,i);
             T(t_count,2) = index_group3(1,j);
             T(t_count,3) = index_group5(1,k);
            t_count = t_count + 1;
        end
       end  
    end
end



for i = 1:9
    for j = 1:9
        if i-j        
        %writing loops for taking the alphabets except for that particular
        %group
        %for group 1
        for k = 1:1
            T(t_count,1) = index_group4(1,i);
            T(t_count,2) = index_group4(1,j);
            T(t_count,3) = index_group1(1,k);
            t_count = t_count + 1;
        end
        %for group 2
         for k = 1:5
            T(t_count,1) = index_group4(1,i);
            T(t_count,2) = index_group4(1,j);
            T(t_count,3) = index_group2(1,k);
            t_count = t_count + 1;
         end
         %for group 3
         for k = 1:9
            T(t_count,1) = index_group4(1,i);
            T(t_count,2) = index_group4(1,j);
            T(t_count,3) = index_group3(1,k);
            t_count = t_count + 1;
         end
         %for group 5
         for k = 1:2
             T(t_count,1) = index_group4(1,i);
             T(t_count,2) = index_group4(1,j);
             T(t_count,3) = index_group5(1,k);
            t_count = t_count + 1;
        end
       end  
    end
end




for i = 1:2
    for j = 1:2
        if i-j        
        %writing loops for taking the alphabets except for that particular
        %group
        %for group 1
        for k = 1:1
            T(t_count,1) = index_group5(1,i);
            T(t_count,2) = index_group5(1,j);
            T(t_count,3) = index_group1(1,k);
            t_count = t_count + 1;
        end
        %for group 2
         for k = 1:5
            T(t_count,1) = index_group5(1,i);
            T(t_count,2) = index_group5(1,j);
            T(t_count,3) = index_group2(1,k);
            t_count = t_count + 1;
         end
         %for group 3
         for k = 1:9
            T(t_count,1) = index_group5(1,i);
            T(t_count,2) = index_group5(1,j);
            T(t_count,3) = index_group3(1,k);
            t_count = t_count + 1;
         end
         %for group 4
         for k = 1:9
             T(t_count,1) = index_group5(1,i);
             T(t_count,2) = index_group5(1,j);
             T(t_count,3) = index_group4(1,k);
            t_count = t_count + 1;
        end
       end  
    end
end




num_count = 0;
for i = 1: t_count-1
    d12 = sqrt((X(T(i,1),:)-X(T(i,2),:))*M*transpose((X(T(i,1),:)-X(T(i,2),:))));
    d13 = sqrt((X(T(i,1),:)-X(T(i,3),:))*M*transpose((X(T(i,1),:)-X(T(i,3),:))));
    if d12<d13
        num_count = num_count + 1;
    end
end
acc = num_count/t_count;


