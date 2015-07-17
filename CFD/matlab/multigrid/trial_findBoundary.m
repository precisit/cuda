close all
clc
clear all

n = 5;

allNodes = zeros(n*n, 2);

counter = 0;

for x=1:n
    for y=1:n
        counter = counter + 1;
        allNodes(counter, 1) = x;
        allNodes(counter, 2) = y;
    end
end

 
