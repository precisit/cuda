function [ A_out ] = interpolation( A_vec )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

A = zeros(sqrt(length(A_vec)));
n = length(A);
for x=0:n-1
    for y=0:n-1
        A(x +1, y +1) = A_vec(y+x*n +1);
    end
end


A_new = zeros(2*n-1, 2*n-1);

for x=1:n
    for y=1:n
        A_new(2*x-1,2*y-1) = A(x,y); %this will need a rewrite to work in 0-indexing
    end
end


for x=1:n-1
    for y=1:n
        A_new(2*x,2*y-1) = (A(x,y)+A(x+1,y))/2; %this will need a rewrite to work in 0-indexing
    end
end

for x=1:2*n-1
    for y=2:2:2*n-1
        A_new(x,y) = (A_new(x,y-1)+A_new(x,y+1))/2;
    end
end


A_out = zeros(length(A_new)^2,1);
for x=0:n
    for y=0:n
        A_out(y+x*n +1) = A_new(y +1,x +1);
    end
end

end

