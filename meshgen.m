function [point,boundary,inner,freenodes,mesh, inner_inner] = meshgen(m)
%MESHGEN Generate mesh of size m x m
%
%   Input:
%           m : size of mesh grid
%   Ouput:
%           p: coordinates of grid points
%           mesh: triangle elements

[x,y]=ndgrid((0:(m-1))/(m-1),(0:(m-1))/(m-1));
point = [x(:),y(:)];
mesh=[1,2,m+2; 1,m+2,m+1];
mesh=kron(mesh,ones(m-1,1))+kron(ones(size(mesh)),(0:m-2)');
mesh=kron(mesh,ones(m-1,1))+kron(ones(size(mesh)),(0:m-2)'*m);
boundary=[1:m-1,m+1:m:m*m,m:m:m*m-1,m*m-m+2:m*m];
inner= [m+1:2*m-1,m+2:m:m*m,m-1:m:m*m-m,m*m-2*m+2:m*m-m];
inner_inner = [2*m+1:3*m-1,m+3:m:m*m, m-2:m:m*m-m, m*m-3*m+2:m*m-2*m];
freenodes=setdiff(1:m^2,boundary);
end

