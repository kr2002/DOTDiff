%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function setup boundary conditions for the 
% adjoint diffusion problem.
%
% The adjoint diffusion model:
%
% -\nabla\cdot\gamma\nabla w + \sigma w = 0  in \Omega
% \bnu\cdot\gamma\nabla w +kappa w=g, on \partial\Omega
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qmatrix,gmatrix,hmatrix,rmatrix] = DiffBCAdj(p,e,u,time,kappa,g,ks,srcdetpair)

ne = size(e,2); % number of edges
hmatrix = zeros(1,2*ne);
rmatrix = zeros(1,2*ne);

qmatrix = kappa*ones(1,ne);
gmatrix = zeros(1,ne);

for k = 1:ne

	x1 = p(1,e(1,k)); % x at first point in segment
	y1 = p(2,e(1,k)); % y at first point in segment
	x2 = p(1,e(2,k)); % x at second point in segment
	y2 = p(2,e(2,k)); % y at second point in segment
	xm = (x1 + x2)/2; % x at segment midpoint
	ym = (y1 + y2)/2; % y at segment midpoint

	gmatrix(k) = g(k); 
    
end