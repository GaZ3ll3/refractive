function [ R ,E , G ] = regularization(m)
% REGULARIZATION generate regularization operator
%   Define TV regularization
Seed=diag(2*ones(m,1))+diag(-1*ones(m-1,1),1)+diag(-1*ones(m-1,1),-1);
Seed(1,1)=1;
Seed(m,m)=1;

GL=kron(speye(m),Seed);
GR=kron(Seed,speye(m));
G=GR+GL;

R=diag(1*ones(m^2,1))+diag(-1*ones(m^2-1,1),1);
E=diag(1*ones(m^2,1))+diag(-1*ones(m^2-m,1),m);
for l=1:m-1
    R(m*l,:)=0;
    E(m*l,:)=0;
end

for l=1:m
    R(m*(m-1)+l,:)=0;
    E(m*(m-1)+l,:)=0;
end

end

