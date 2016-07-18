function [U, U_nu,refc_M, elapse] = FEM_Helmholtz(M, K, m,k,force,Dirichlet,refc,boundary,inner_boundary,freenodes,mesh)
%   FEM_HELMHOLTZ Solve Helmholtz equation with force source.
%   Rebuilt.
%   Input: 
%           M: mass matrix
%           K: stiff matrix
%           m: size of mesh
%           k: wave number 
%           force: source vector evaluate at each point.
%           refc: refractive index vecetor of the media.
%           Dirichlet: Dirichlet boundary condition vector, defined on
%           whole domain. That means the size is m^2.
%   Output:
%           U: solution of Helmoltz equation
%   time counter
elapse = cputime;
%   number of nodes, number of triangles
N = m^2;T = size(mesh,1); 
%   involving refractive index in mass matrix
refc_M = sparse(N,N);
%   generate solution vector
U = zeros(N,1);
%   generate matrix by element, sparse indexing may take time.
for e = 1:T  
  nodes = mesh(e,:); 
  Refc_i = (1/120/(m-1)^2)*[[6 2 2];[2 2 1];[2 1 2]];
  Refc_j = (1/120/(m-1)^2)*[[2 2 1];[2 6 2];[1 2 2]];
  Refc_k = (1/120/(m-1)^2)*[[2 1 2];[1 2 2];[2 2 6]];
  refc_M(nodes,nodes) = refc_M(nodes,nodes) + [Refc_i*refc(nodes), Refc_j*refc(nodes), Refc_k*refc(nodes)]; 
end 
%   if force is not point source, we can do this.
%   generate load vector by 

F = M*force;

%   Dirichlet boundary condition
U(boundary) = Dirichlet(boundary);
F = F + (K-k^2*refc_M)*U;
U(freenodes) = (-K(freenodes,freenodes) + k^2*refc_M(freenodes,freenodes))\F(freenodes);
U_nu = sparse(m^2,1);
U_nu(boundary) = (U(boundary) - U(inner_boundary))*(m-1);
elapse = cputime - elapse;
end