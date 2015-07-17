close all
clear all
%clc

k = 2;
n = 2^k+1;

A = discMat(n);

b = zeros(n*n,1);
x = ones(size(b));

%for iter = 1:1
%    x = mgcyc(k,2,x,A,b,100,100);
%end
%
%sum(abs(x).^2)/length(x).^2

x = gaussSeidel(A,x,b)
