function [M, K] = buildmk(m,point,mesh)
%   FEM_HELMHOLTZ Solve Helmholtz equation with force source.

%   Input: 
%           m: size of mesh
%           k: wave number 
%           force: source vector evaluate at each point.
%           refc: refractive index vecetor of the media.
%           Dirichlet: Dirichlet boundary condition vector, defined on
%           whole domain. That means the size is m^2.
%   Output:
%           U: solution of Helmoltz equation


%   Setup mesh, generate grid of m x m
%[point,boundary,inner_boundary,freenodes,mesh] = meshgen(m);
%   number of nodes, number of triangles
N = m^2;T = size(mesh,1); 
%   generate stiff matrix and mass matrix
K = sparse(N,N); 
M = sparse(N,N);
for e = 1:T  
  nodes = mesh(e,:); 
  Pe = [ones(3,1),point(nodes,:)]; 
  Area = abs(det(Pe))/2;
  C = inv(Pe); 
  grad = C(2:3,:);Ke = Area*(grad'*grad); 
  Me = Area*(1/12)*[[2 1 1];[1 2 1 ];[1 1 2]];
  K(nodes,nodes) = K(nodes,nodes) + Ke; 
  M(nodes,nodes) = M(nodes,nodes) + Me;  
end 
end


