function visualize(refc,refc_true)
%VISUALIZE 
clf;
m = sqrt(size(refc_true,1));
[point,~,~,~,mesh] = meshgen(m);
%subplot(3,2,1);
figure(1);
trisurf(mesh,point(:,1),point(:,2),0*refc,real(refc_true),'edgecolor','none','facecolor','interp');view(2);colorbar;title('true refractive index');axis('square');
%subplot(3,2,3);
figure(2);
trisurf(mesh,point(:,1),point(:,2),0*refc,real(refc),'edgecolor','none','facecolor','interp');view(2);colorbar;title('numerical refractive index');axis('square');
%subplot(3,2,5);
figure(3);
trisurf(mesh,point(:,1),point(:,2),0*refc,real(refc-refc_true)./real(refc_true),'edgecolor','none','facecolor','interp');view(2);colorbar;title('relative error');axis('square');
%subplot(3,2,2);
figure(4);
plot(refc(floor(m/2)*m+1:floor(m/2)*m+m),'-b');hold on;plot(refc_true(floor(m/2)*m+1:floor(m/2)*m+m),'-r');title('cross in middle x-axis');axis('square');
%subplot(3,2,4);
figure(5);
plot(refc(floor(m/2):m:m^2),'-b');hold on;plot(refc_true(floor(m/2):m:m^2),'-r');title('cross in middle y-axis');axis('square');
%subplot(3,2,6);
figure(6);
plot(refc(1:m+1:m^2),'-b');hold on;plot(refc_true(1:m+1:m^2),'-r');title('cross in diagonal');axis('square');
end

