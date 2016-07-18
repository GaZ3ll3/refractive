function [source] = Refc_source_Solver()

global point boundary  inner_boundary freenodes  mesh  m  k refc_true refc source_true source_num  source_num_true M  K  R E G Denom count
%   set up mesh grid
 m = 40;
% frequency, wave number
 k = 2;
[point,boundary,inner_boundary,freenodes,mesh] = meshgen(m);
[M,K] = buildmk(m,point,mesh);
% refc_true on fine mesh
%refc_true = 1 + (0.5+0.5*cos(sqrt((point(:,1)-0.5).^2+(point(:,2)-0.5).^2)*pi*5)).*((point(:,1)-0.5).^2+(point(:,2)-0.5).^2<0.04);
% load('refc_true_corner_single_mesh_40_2__1_dir_6_radius_0.25.mat');
% load('refc_corner_single_mesh_40_2__1_dir_6_radius_0.25.mat');
% load('refc_mesh_40_2__1_dir_4_radii_0.4.mat')
% load('refc_true_mesh_40_2__1_dir_4_radii_0.4.mat')
% load('refc_mesh_40_2__1_dir_12_radius_0.2.mat')
% load('refc_true_mesh_40_2__1_dir_12_radii_0.2.mat')
% load('refc_mesh_40_diag_2__1_dir_6_radii_0.3.mat')
% load('refc_true_mesh_40_diag_2__1_dir_6_radius_0.3.mat')
load('refc_true_mesh_40_quad_2__1_dir_6_radius_0.3.mat');
load('refc_mesh_40_quad_2__1_dir_6_radius_0.3.mat');
% load('refc_true_mesh_40_jump_dir_6.mat');
% load('refc_mesh_40_jump_dir_6.mat');
%refc =refc_true;
source_true = [0.2;0.2;10;0.1;0.55;15;0.45;0.9;12;0.8;0.6;10;0.9;0.3;10];
options = optimset('Diagnostics','on','DerivativeCheck','off','FinDiffType','central','LargeScale','off',...
    'GradObj','on','Display','iter-detailed','TolFun',1e-16,'TolX',1e-16,'MaxFunEvals',40000,'MaxIter',40000,'HessUpdate','bfgs','LineSearchType','cubicpoly');

source = [0.5;0.5;10;0.2;0.4;10;0.5;0.8;10;0.8;0.4;10;0.8;0.5;10];


Denom = 32;
count = 0;

source_num = length(source)/3;
source_num_true = length(source_true)/3;
 

%   cutoff the sides
%surf(reshape(force_true,m,m));view(2);title('real source');colorbar;
%   TODO.

[source,fval,exitflag,output] = fminunc(@Refcopt_Source,source,options);
fprintf('final result: %s\n\n', num2str(fval));
fprintf('exitflag: %s\n\n', num2str(exitflag));
display(output);

if source_num == source_num_true
    % REARRANGE THE ORDER ,THEN COMPARE.
    inten = sort(source(3:3:end));
    x_cor = sort(source(1:3:end));
    y_cor = sort(source(2:3:end));    
    tr_inten = sort(source_true(3:3:end));
    tr_x_cor = sort(source_true(1:3:end));
    tr_y_cor = sort(source_true(2:3:end));
    fprintf('error_x:%s\nerror_y:%s\nerror_h:%s\n',num2str(max(abs(tr_x_cor-x_cor))),num2str(max(abs(tr_y_cor-y_cor)))...
    ,num2str(max(abs(inten-tr_inten))));
end

%   UNEQUAL CASES FOR TESTING
if source_num == 2 && source_num_true == 1
    fprintf('error_x:%s\nerror_y:%s\nerror_h:%s\n',num2str(max(abs(source_true(1)-source(1)),abs(source_true(1)-source(4)))),...
        num2str(max(abs(source_true(2)-source(2)),abs(source_true(2)-source(5))))...
    ,num2str(abs(source_true(3)-source(3)-source(6))));
end

if source_num == 1 && source_num_true == 2
    fprintf('weighted error_x:%s\nweighted error_y: %s\nweighted error_h:%s',...
        num2str(abs((source_true(3)*source_true(1)+source_true(6)*source_true(4))/(source_true(3)+source_true(6))-source(1))),...
        num2str(abs((source_true(3)*source_true(2)+source_true(6)*source_true(5))/(source_true(3)+source_true(6))-source(2))),...
        num2str(abs(source_true(3)+source_true(6)-source(3))));
end

end

