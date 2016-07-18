function [U, U_nu,refc_M, elapse] = point_solver(M,K,m,k,smooth_force,point_force,Dirichlet,refc,boundary,inner_boundary,freenodes,mesh)
%   Solve second order equation with Dirichlet boundary condition
%   D(Du) + k^2 refc u = force
%                    u = Dirichlet
%   in Helmholtz equation form. k and refc can be arbitrary.
%
%   Rebuilt.
%
%   Input: 
%           M: mass matrix
%           K: stiff matrix
%           m: size of mesh, default value as 41.
%           k: wave number 
%           smooth_force: smooth source vector evaluate at each point.
%           point_force : point sources, with location and intensity,
%           sparse.
%           refc: refractive index vecetor of the media.
%           Dirichlet: Dirichlet boundary condition vector, defined on
%           whole domain. That means the size is m^2.
%           flag : indicates the force is point source or not.
%           mesh : the coordinates of points in dimain
%           innerboundary: the boundary after removing the real boundary.
%
%   Output:
%           U: solution.
%           U_nu : normal derivative at boundary.
%           refc_M: mass matrix associated with refc
%           elapse: CPU time on generating.
%
%
%   time counter
%
if nargin == 0
m= 41; k =1; N= m^2 ;
[point,boundary,inner_boundary,freenodes,mesh] = meshgen(m);
[M,K] = buildmk(m,point,mesh);
%smooth
smooth_force = zeros(N,1);
%point
point_force = [-1;0.5;0.5];
%translation
X = point(:,1)-.5;
Y = point(:,2)-.5;
%Dirichlet
Dirichlet = zeros(N,1);
%refc
refc = ones(N,1);
%solve the equation
end
elapse = cputime;
%   number of nodes, number of triangles
N = m^2;T = size(mesh,1); 
%   involving refractive index in mass matrix
refc_M = sparse(N,N);
%   generate solution vector
U = zeros(N,1);
%   generate matrix by element, sparse indexing may take time.
%
for e = 1:T  
  nodes = mesh(e,:); 
  Refc_i = (1/120/(m-1)^2)*[[6 2 2];[2 2 1];[2 1 2]];
  Refc_j = (1/120/(m-1)^2)*[[2 2 1];[2 6 2];[1 2 2]];
  Refc_k = (1/120/(m-1)^2)*[[2 1 2];[1 2 2];[2 2 6]];
  refc_M(nodes,nodes) = refc_M(nodes,nodes) + [Refc_i*refc(nodes), Refc_j*refc(nodes), Refc_k*refc(nodes)]; 
end 
%   if force is not point source, we can do this.
%   generate load vector by 
smooth_F = M*smooth_force;
%
    %force take value as point sources.
    %force should have form as a 3xL matrix, has intensity, coordinates.
quant = size(point_force,2);
point_F = sparse(N,1);
for i = 1:quant        
    %locate the point
    xnum = floor(point_force(2,i)*(m-1));
    ynum = floor(point_force(3,i)*(m-1));
    if xnum > m || xnum<0 ||ynum>m ||ynum<0
        continue
    end
    xextra = point_force(2,i)*(m-1)-xnum;
    yextra = point_force(3,i)*(m-1)-ynum;%scaling to be accurate
    
    if xextra>=yextra
        % in lower triangle            
        point_F((ynum-1)*m+xnum) = point_F((ynum-1)*m+xnum)+(1-xextra)*point_force(1,i);
        point_F((ynum-1)*+xnum+1) = point_F((ynum-1)*+xnum+1) + (xextra-yextra)*point_force(1,i);
        point_F(ynum*m+xnum+1) = point_F(ynum*m+xnum+1)+yextra*point_force(1,i);
    else
        point_F((ynum-1)*m+xnum) = point_F((ynum-1)*m+xnum) + (1-yextra)*point_force(1,i);
        point_F(ynum*m+xnum) = point_F(ynum*m+xnum) + (yextra-xextra)*point_force(1,i);
        point_F(ynum*m+xnum+1) = point_F(ynum*m+xnum+1) + xextra*point_force(1,i);
    end
end
%   Dirichlet boundary condition
%
U(boundary) = Dirichlet(boundary);
F = smooth_F+point_F;
F = F + (K-k^2*refc_M)*U;
U(freenodes) = (-K(freenodes,freenodes) + k^2*refc_M(freenodes,freenodes))\F(freenodes);
U_nu = sparse(m^2,1);
U_nu(boundary) = (U(boundary) - U(inner_boundary))*(m-1);
elapse = cputime - elapse;
end

