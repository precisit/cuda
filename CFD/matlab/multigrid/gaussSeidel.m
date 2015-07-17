% function [ x_new ] = gaussSeidel( A, x, b )
% %GaussSeidel Does a GS smoothing of the equation Ax=b
% %   
% 
% n = length(x);
% 
% x_new = zeros(size(x));
% 
% for i= 1:n
%    
%     for j = 1:(i-1)
%         x_new(i) = x_new(i) - A(i,j)*x_new(j);    
%     end
%     
%     for j=(i+1):n
%        x_new(i) = x_new(i) - A(i,j)*x(j); 
%     end
%     
%     x_new(i) = x_new(i) + b(i);
%     x_new(i) = x_new(i) / A(i,i);
%     
% end
% 
% 
% end

function [ x ] = gaussSeidel( A, x, b )
%GaussSeidel Does a GS smoothing of the equation Ax=b
%

n = length(x);

for i= 1:n
   
    for j = 1:(i-1)
        x(i) = x(i) - A(i,j)*x(j);    
    end
    
    for j=(i+1):n
       x(i) = x(i) - A(i,j)*x(j); 
    end
    
    x(i) = x(i) + b(i);
    x(i) = x(i) / A(i,i);
    
end


end

% function [x] = gauss_seidel(A, b, x0, iters)
%     n = length(A);
%     x = x0;
%     for k = 1:iters
%         for i = 1:n
%             x(i) = (1/A(i, i))*(b(i) - A(i, 1:n)*x + A(i, i)*x(i));
%         end
%     end
% end