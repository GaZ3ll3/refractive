function [mat_new] = refinement(mat)
%REFINEMENT refine mesh by decrease size to half and linear interpolation
%   Input:
%           mat : input matrix
%   Output:
%           mat_new : output matrix with higher dimension.
m = size(mat);

n = 2*m - 1;

mat_new = zeros(n,n);

for i = 1 : m
    for j = 1 : m
        mat_new (2*i-1,2*j-1) = mat(i,j);
    end
end

for i = 1: m
    for j = 1 : m
        if (j<m)
            mat_new (2*i-1 , 2*j) = (mat_new(2*i-1,2*j-1)+mat_new(2*i-1,2*j+1))/2;
        end
        if (i<m)
            mat_new (2*i , 2*j-1) = (mat_new(2*i+1,2*j-1)+mat_new(2*i-1,2*j-1))/2;
        end
        if i<m & j<m
        	mat_new (2*i ,   2*j) = (mat_new(2*i+1,2*j+1)+mat_new(2*i+1,2*j-1)+mat_new(2*i-1,2*j-1)+mat_new(2*i-1,2*j+1))/4;
        end
    end
end
end

