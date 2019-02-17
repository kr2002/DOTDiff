%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DiffSolve: FEM solver for forward diffusion equation
%
% The diffusion model:
%
% -\nabla\cdot\gamma\nabla u + \sigma u = 0  in \Omega
% \bnu\cdot\gamma \nabla u+\kappa u = f, on \partial \Omega
%
% The measurement quantity:
% 
% g=u  on \partial\Omega
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u meas]=DiffSolve(P,E,T,gamma,sigma,kappa,F,...
    ks,Nd,srcinfo,detinfo)

% interpolation to triangle middle point
gammam=pdeintrp(P,T,gamma);
sigmam=pdeintrp(P,T,sigma);
Fm=pdeintrp(P,T,F);

% construct boundary conditions
pdebound =@(p,e,u,time)DiffBC(p,e,[],[],kappa,ks,srcinfo,detinfo);
[Q,G,H,R] = assemb(pdebound,P,E);

% construct mass matrices
[K,M,F]=assema(P,T,gammam,sigmam,Fm);

% solve the PDE
u = assempde(K,M,F,Q,G,H,R);
    
% solve PDEs (old way of solving PDE)
%u=assempde(geobc,P,E,T,gamma,sigmam,Fm);

% compute measured data
meas=zeros(1,Nd);
for j=1:Nd
    meas(j)=u(detinfo(1,j));
end