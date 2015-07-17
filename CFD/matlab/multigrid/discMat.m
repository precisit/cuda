function [ L ] = discMat( n )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

L = zeros(n*n, n*n);

for i=0:n*n-1 %this shit HAS to be 0-indexed
    
   %convert to x and y 
    y = mod(i,n);
    x = (i-y) / n;
    
    
    %The one above
    if(y+1 >= n)
        %do something bc related here. I don't care.
    else
        L(i +1,(y+1)+x*n +1) = L(i +1,(y+1)+x*n +1) +1;
    end
    
    
    %The one below
    if(y-1 < 0)
        %do something bc related here. I don't care.
    else
        L(i +1,(y-1)+x*n +1) = L(i +1,(y-1)+x*n +1) +1;
    end
    
    
    %The one to the right
    if(x+1 >= n)
        %do something bc related here. I don't care.
    else
        L(i +1,y+(x+1)*n +1) = L(i +1,y+(x+1)*n +1) +1;
    end
    
    
    %The one to the left
    if(x-1 < 0)
        %do something bc related here. I don't care.
    else
        L(i +1,y+(x-1)*n +1) = L(i +1,y+(x-1)*n +1) +1;
    end
    
    %And itself.
    L(i +1,i +1) = L(i +1,i +1) - 4;
    
end

    

end

