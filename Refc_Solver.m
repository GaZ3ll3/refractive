function refc = Refc_Solver(refc_init)
%REFC_SOLVER Solve refractive index through Newton's method[bfgs].
%   Output:
%           refc: refractive index vector at each grid point.
global point pointic boundary boundaryic inner_boundary inner_boundaryic inner_inner_ic freenodes freenodesic mesh meshic m k  mic  refc_true refc_trueic M Mic K Kic R E G


%-------------------------------FILE CHECKING------------------------------

%   Avoid crash or overlap, program needs sleep delay on file checking
pause(1);
%default as 1 seconds, with overheads, the system will have enough time on
%checking files.

for refc_type =refc_init:7
    filename=sprintf('%d.lock',refc_type);
    if ~exist(filename,'file')
        % generate the lock file
        fid = fopen(filename,'w');
        fclose(fid);
        break;
    end
end %while

disp('refc type is:');
disp(refc_type);
disp('\n'); 
% --------------------------BEGINNING OF PROGRAM---------------------------





%matlabpool open local 4;
%   set up mesh grid
 m = 79;
%   for inverse crime
 mic = 2*m-1;
%   wave number
 k = 8;
[point,boundary,inner_boundary,freenodes,mesh] = meshgen(m);
[pointic,boundaryic,inner_boundaryic,freenodesic,meshic,inner_inner_ic] = meshgen(mic);

!
[M,K] = buildmk(m,point,mesh);
[Mic,Kic] = buildmk(mic,pointic,meshic);
[R,E,G] = regularization(m);
%   exact refractive index
% TYPE-D-0.25
if refc_type==1
    radius = 0.25;
    refc_trueic =1+...
         (0.5+0.5*cos(sqrt((pointic(:,1)-0.75).^2+(pointic(:,2)-0.75).^2)*pi*1/radius)).*((pointic(:,1)-0.75).^2+(pointic(:,2)-0.75).^2<radius^2)+...
         (0.5+0.5*cos(sqrt((pointic(:,1)-0.25).^2+(pointic(:,2)-0.25).^2)*pi*1/radius)).*((pointic(:,1)-0.25).^2+(pointic(:,2)-0.25).^2<radius^2)+...
         (0.0+0.0*cos(sqrt((pointic(:,1)-0.75).^2+(pointic(:,2)-0.25).^2)*pi*1/radius)).*((pointic(:,1)-0.75).^2+(pointic(:,2)-0.25).^2<radius^2)+...
         (0.0+0.0*cos(sqrt((pointic(:,1)-0.25).^2+(pointic(:,2)-0.75).^2)*pi*1/radius)).*((pointic(:,1)-0.25).^2+(pointic(:,2)-0.75).^2<radius^2);
end
%TYPE-D-0.2
if refc_type==2
    radius = 0.2;
    refc_trueic =1+...
         (0.5+0.5*cos(sqrt((pointic(:,1)-0.75).^2+(pointic(:,2)-0.75).^2)*pi*1/radius)).*((pointic(:,1)-0.75).^2+(pointic(:,2)-0.75).^2<radius^2)+...
         (0.5+0.5*cos(sqrt((pointic(:,1)-0.25).^2+(pointic(:,2)-0.25).^2)*pi*1/radius)).*((pointic(:,1)-0.25).^2+(pointic(:,2)-0.25).^2<radius^2)+...
         (0.0+0.0*cos(sqrt((pointic(:,1)-0.75).^2+(pointic(:,2)-0.25).^2)*pi*1/radius)).*((pointic(:,1)-0.75).^2+(pointic(:,2)-0.25).^2<radius^2)+...
         (0.0+0.0*cos(sqrt((pointic(:,1)-0.25).^2+(pointic(:,2)-0.75).^2)*pi*1/radius)).*((pointic(:,1)-0.25).^2+(pointic(:,2)-0.75).^2<radius^2);
end
%TYPE-S-0.25
% if refc_type==3
%     radius=0.25;
%     refc_trueic =1+...
%          (0.5+0.5*cos(sqrt((pointic(:,1)-0.65).^2+(pointic(:,2)-0.65).^2)*pi*1/radius)).*((pointic(:,1)-0.65).^2+(pointic(:,2)-0.65).^2<radius^2)+...
%          (0.0+0.0*cos(sqrt((pointic(:,1)-0.25).^2+(pointic(:,2)-0.25).^2)*pi*1/radius)).*((pointic(:,1)-0.25).^2+(pointic(:,2)-0.25).^2<radius^2)+...
%          (0.0+0.0*cos(sqrt((pointic(:,1)-0.75).^2+(pointic(:,2)-0.25).^2)*pi*1/radius)).*((pointic(:,1)-0.75).^2+(pointic(:,2)-0.25).^2<radius^2)+...
%          (0.0+0.0*cos(sqrt((pointic(:,1)-0.25).^2+(pointic(:,2)-0.75).^2)*pi*1/radius)).*((pointic(:,1)-0.25).^2+(pointic(:,2)-0.75).^2<radius^2);
% end
%TYPE-C-0.2
if refc_type==3
    radius=0.2;
    refc_trueic = 1 + (0.5+0.5*cos(sqrt((pointic(:,1)-0.5).^2+(pointic(:,2)-0.5).^2)*pi/radius)).*((pointic(:,1)-0.5).^2+(pointic(:,2)-0.5).^2<radius^2);
end
%TYPE-C-0.3
if refc_type==4
    radius=0.3;
    refc_trueic = 1 + (0.5+0.5*cos(sqrt((pointic(:,1)-0.5).^2+(pointic(:,2)-0.5).^2)*pi/radius)).*((pointic(:,1)-0.5).^2+(pointic(:,2)-0.5).^2<radius^2);
end
%TYPE-C-0.4
if refc_type==5
    radius=0.4;
    refc_trueic = 1 + (0.5+0.5*cos(sqrt((pointic(:,1)-0.5).^2+(pointic(:,2)-0.5).^2)*pi/radius)).*((pointic(:,1)-0.5).^2+(pointic(:,2)-0.5).^2<radius^2);
end
%
if refc_type==6
    refc_trueic = 1+0.2*((pointic(:,1)-0.25).^2+(pointic(:,2)-0.5).^2<0.04)+0.4*((pointic(:,1)<0.75).*(pointic(:,1)>0.5).*(pointic(:,2)<0.75).*(pointic(:,2)>0.25))+...
        0.0*(1-((pointic(:,1)-0.25).^2+(pointic(:,2)-0.5).^2<0.01)).*(1-((pointic(:,1)<0.75).*(pointic(:,1)>0.5).*(pointic(:,2)<0.75).*(pointic(:,2)>0.25)));
end

if refc_type==7
    radius = 0.25;
    refc_trueic =1+...
         (0.5+0.5*cos(sqrt((pointic(:,1)-0.75).^2+(pointic(:,2)-0.75).^2)*pi*1/radius)).*((pointic(:,1)-0.75).^2+(pointic(:,2)-0.75).^2<radius^2)+...
         (0.5+0.5*cos(sqrt((pointic(:,1)-0.25).^2+(pointic(:,2)-0.25).^2)*pi*1/radius)).*((pointic(:,1)-0.25).^2+(pointic(:,2)-0.25).^2<radius^2)+...
         (0.5+0.5*cos(sqrt((pointic(:,1)-0.75).^2+(pointic(:,2)-0.25).^2)*pi*1/radius)).*((pointic(:,1)-0.75).^2+(pointic(:,2)-0.25).^2<radius^2)+...
         (0.5+0.5*cos(sqrt((pointic(:,1)-0.25).^2+(pointic(:,2)-0.75).^2)*pi*1/radius)).*((pointic(:,1)-0.25).^2+(pointic(:,2)-0.75).^2<radius^2);
end



refc_trans = derefinement(reshape(refc_trueic,mic,mic));
refc_true = reshape(refc_trans,m^2,1);

try 
    outputname=sprintf('%d.mat',refc_type);
    load(outputname,'refc');
catch err
    if strcmp(err.identifier,'MATLAB:load:couldNotReadFile')
        refc = 1 + (0.0+0.0*cos(sqrt((point(:,1)-0.5).^2+(point(:,2)-0.5).^2)*pi*5)).*((point(:,1)-0.5).^2+(point(:,2)-0.5).^2<0.04);
        %refc = 1 + (0.5+0.5*cos(sqrt((point(:,1)-0.5).^2+(point(:,2)-0.5).^2)*pi*5/2)).*((point(:,1)-0.5).^2+(point(:,2)-0.5).^2<0.16);
    end
end
% show result.
Diff(refc, refc_true);
% end of showing result


% refc = refc_true;
%   Newton's iteration
options = optimset('Diagnostics','off','DerivativeCheck','off','FinDiffType','central','LargeScale','off',...
    'GradObj','on','Display','iter-detailed','TolFun',1e-12,'TolX',1e-12,'MaxIter',60,'HessUpdate','bfgs');
[refc_st,fval,~,output] = fminunc(@Refcopt,refc,options);
%   iterations to deform
refc_id = fopen('refc.dat','a+');
counter = 0;
%while (refc_st~=refc)
    refc = refc_st;
    save( outputname, 'refc' );
    counter = counter + 1;
%    [refc_st,fval,exitflag,output] = fminunc(@Refcopt,refc,options);  
    fprintf(refc_id,'%s : objective\t %s \t error\t %6.4f\n',num2str(counter),num2str(fval),norm(refc-refc_true)/norm(refc_true));
%end
%refc = refc_st;
display(output);
%   check relative error L2 and L-inf
Diff(refc,refc_true);
% figure(2)
% trisurf(mesh,point(:,1),point(:,2),real(refc),real(refc),'edgecolor','none','facecolor','interp');
% view(3),axis equal;colorbar,title('numerical refractive index');
fclose(refc_id);
%matlabpool close;
end

