function [mat_new] = derefinement(mat)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[m,~] = size(mat);
% m should be odd number

n = round((m+1)/2);

mat_new = zeros(n,n);

for i = 1 : n
    for j = 1 : n
        mat_new (i,j) = mat(2*i-1,2*j-1);
    end
end

end

