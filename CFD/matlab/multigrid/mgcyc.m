function [ u ] = mgcyc(k, gamma, u, L, f, v1, v2)
%MGCYC A multigrid cycle.
%   Detailed explanation goes here
k

%pre-smoothing
for iter = 1: v1
    u = gaussSeidel(L, u, f);
end

%Calculate error
d = f-L*u;
%Go to a smaller grid
d = restriction(d);

L_small = discMat(2^(k-1)+1);

if(k==1)
    v = L_small\d;
else
    for iter =1:gamma
        v = mgcyc(k-1, gamma, zeros(size(d)), L_small, d, v1, v2);
    end
end
k
v = interpolation(v);
u = u + v;

%post-smoothing
for iter =1:v2
   u=gaussSeidel(L,u,f); 
end

end

