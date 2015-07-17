function [ A_new ] = restriction( A )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
n = length(A);
k = log(sqrt(n)-1)/log(2);
n = sqrt(n);

A_new = zeros( ( 2^(k-1)+1 )^2, 1);

counter = 1;
for iter1 = 0:2:n-1
    for iter2 = 0:n-1
        if( mod(iter2,2) == 0 )
            A_new(counter) = A(iter1*n+iter2+1);
            counter = counter + 1;
        end
    end
end

end

