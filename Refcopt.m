function [f, g] = Refcopt( refc )
%REFCOPT Provide user supply gradient and objective function
%   Input:
%           refc: refractive index vector
%   Output:
%           f: objective function
%           g: gradient vector at refc

global m mic k refc_true refc_trueic point pointic boundary boundaryic inner_boundary inner_boundaryic inner_inner_ic freenodes freenodesic mesh meshic M Mic K Kic R E G
%   make movie of refc
% figure(1);
% surf(reshape(refc,m,m));colorbar;

%   Defind maximal scattering directions.
max_dir = 12;
%   Define regularization level.
gamma = 1e-7;
%   Gradient and objective function
f = 0;
g = zeros(m^2,1);
Q = zeros(m^2);
% noise argument
noise = 0.05;
try 
    load('NOISE.mat','NOISE');
catch err
    if strcmp(err.identifier,'MATLAB:load:couldNotReadFile')
        NOISE = (1+noise*2*(rand(size(pointic,1),1)-0.5));
        save('NOISE.mat','NOISE');
    end
end

%   force pre-defined, not necessary in this case.
force = zeros(m^2,1);
forceic = zeros(mic^2,1);
%   Gradient not accurate, divide into real and imag part.
for ind = 1 : max_dir
    direction = [cos(ind*2*pi/max_dir) sin(ind*2*pi/max_dir)];
    boundary_condition = exp(1i*k*point*direction');
    boundary_conditionic = exp(1i*k*pointic*direction');%data
    [~,U_true_nuic] = FEM_Helmholtz(Mic,Kic,mic,k,forceic,boundary_conditionic,refc_trueic,boundaryic,inner_inner_ic,freenodesic,meshic);
    %noise od data
    U_true_nuic = U_true_nuic.*NOISE;
    U_true_nu = reshape(derefinement(reshape(U_true_nuic,mic,mic)),m^2,1);
    U_true_nu = U_true_nu/2;
    %[~,U_true_nu] = FEM_Helmholtz(M,K, m,k,force,boundary_condition,refc_true,boundary,inner_boundary,freenodes,mesh);
    [U,U_nu] = FEM_Helmholtz(M,K, m,k,force,boundary_condition,refc,boundary,inner_boundary,freenodes,mesh);
    % corner should be matched, or just ignored.
    U_nu(1) = U_true_nu(1);
    U_nu(m) = U_true_nu(m);
    U_nu(m*m-m+1) = U_true_nu(m*m-m+1);
    U_nu(m*m) = U_true_nu(m*m);
    % matching finished
    W_boundary = conj(U_nu - U_true_nu);
    [W] = FEM_Helmholtz(M,K, m,k,force,W_boundary,refc,boundary,inner_boundary,freenodes,mesh);
    f = f + real((U_nu-U_true_nu)'*(U_nu - U_true_nu));
    % use 4 threads for accelerate the evaluation.
    for index = 1: m^2
        Q(:,index) = U.*M(:,index);
    end
    g = g - 2*k^2*real(Q*W)*(m-1);
end
  % normalize
  f = f/(m^2) + 0.5*gamma*(R*refc)'*(R*refc)+0.5*gamma*(E*refc)'*(E*refc);
  g = g/(m^2) + gamma*G*refc;
end

