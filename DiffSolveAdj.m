%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DiffSolveAdj: FEM solver for the adjoint diffusion
%               equation
%
% The adjoint diffusion model:
%
% -\nabla\cdot\gamma\nabla w + \sigma w = 0  in \Omega
% \bnu\cdot\gamma \nabla w+\kappa w = g, on \partial\Omega
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u=DiffSolveAdj(P,E,T,gamma,sigma,kappa,g,ks,srcdetpair)

% interpolation to triangle middle point
gammam=pdeintrp(P,T,gamma);
sigmam=pdeintrp(P,T,sigma);

% construct boundary conditions
pdebound =@(p,e,u,time)DiffBCAdj(p,e,[],[],kappa,g,ks,srcdetpair);
[Q,G,H,R] = assemb(pdebound,P,E);

% construct mass matrices
[K,M,F]=assema(P,T,gammam,sigmam,0);

% solve the PDE
u = assempde(K,M,F,Q,G,H,R);
    
% solve PDEs (old way of solving the PDE)
%u=assempde(geobc,P,E,T,gammam,sigmam,0);