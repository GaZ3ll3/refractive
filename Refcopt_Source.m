function [f g] = Refcopt_Source(source)
%   REFCOPT_SOURCE Recover source position and amplitude.
%   Input:  
%           source: vector value function with arguments [x_0, y_0, h_0]
%   Output:
%           f: objective function
%           g: gradient of f(source)

%   For reconstruction of source, only one data is needed.
%   Thus choose a random direction and use plain wave.

%   Global variables according to the mesh and initial data as well as
%   regularization level settings
global m k refc_true point boundary inner_boundary freenodes mesh M K source_true refc source_num source_num_true epsilon Denom count
%   If there is no other requirements, suppose there is only one
%   source[point]. We also focus on bi-point source.
%   initialization
f = 0;
g = zeros(source_num*3,1);
%   deviation of source function.
epsilon = (1/m)/Denom;
%   force vector defined for each tuple source.
force_true = zeros(m^2,1);
force = zeros(m^2,source_num);
force_all = zeros(m^2,1);
%   END OF SETTING
for ord = 1 : source_num_true 
    force_true = force_true + source_true(3*ord)*exp(-((point(:,1)-source_true(3*ord-2)).^2+(point(:,2)-source_true(3*ord-1)).^2)/epsilon);
end
%   STORE EACH, TAKES MEMORY!
for ord = 1 : source_num
    force(:,ord) = source(3*ord)*exp(-((point(:,1)-source(3*ord-2)).^2+(point(:,2)-source(3*ord-1)).^2)/epsilon);
    force_all = force_all + force(:,ord);
end
%   END OF SETTINGS
%locate points

%   figure
%-------------------------------
count = count + 1;
new_fig = figure(1);
% surf(reshape(force_all,m,m));title('numerical result');colorbar;
% view(2);
plot(source(1:3:end),source(2:3:end),'bo','MarkerSize',10,'MarkerFaceColor',[.49 1 .63])
hold on
plot(source_true(1:3:end),source_true(2:3:end),'rs','MarkerSize',5,'MarkerFaceColor',[.49 1 .63]);
hold off
set(gca, 'color', [0.5 0.5 0.8]);
axis([0,1,0,1]);
%file_name=sprintf('%d.gif',count);
%saveas(new_fig,file_name)  % here you save the figure

%----------------------------
%   Adjoint function W.
W_force = zeros(m^2,1);
%   Gradient not accurate, divide into real and imag part.
%   Only one direction is needed for recovering tuple source.
%rnd = 0;
%direction = [cos(rnd*2*pi) sin(rnd*2*pi)];
%boundary_condition = exp(1i*k*point*direction');
boundary_condition = zeros(m^2,1);
[~,U_true_nu] = FEM_Helmholtz(M,K,m,k,force_true,boundary_condition,refc_true,boundary,inner_boundary,freenodes,mesh);
[~,U_nu] = FEM_Helmholtz(M,K,m,k,force_all,boundary_condition,refc,boundary,inner_boundary,freenodes,mesh);
W_boundary = conj(U_nu - U_true_nu);
[W] = FEM_Helmholtz(M,K, m,k,W_force,W_boundary,refc,boundary,inner_boundary,freenodes,mesh);
%   Gradient of source function

grad = zeros(m^2,3*source_num);
for ord = 1 : source_num
    grad(:,3*ord-2) = -2*force(:,ord).*(point(:,1)-source(3*ord-2))/epsilon;
    grad(:,3*ord-1) = -2*force(:,ord).*(point(:,2)-source(3*ord-1))/epsilon;
    grad(:,3*ord)   = -exp(-((point(:,1)-source(3*ord-2)).^2+(point(:,2)-source(3*ord-1)).^2)/epsilon);
end
%   use f to force the result as a double type.
f = f + real((U_nu-U_true_nu)'*(U_nu - U_true_nu));

for ord = 1 : 3*source_num
    g(ord) = -2*real(W'*M*grad(:,ord)*(m-1));
end

