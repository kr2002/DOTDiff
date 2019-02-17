%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function setup boundary conditions for the 
% diffusion problem.
%
% The diffusion model:
%
% -\nabla\cdot\gamma\nabla u + \sigma u = F  in \Omega
% \bnu\cdot\gamma\nabla u+\kappa u = f, on \partial\Omega
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qmatrix,gmatrix,hmatrix,rmatrix] = DiffBC(p,e,u,time,kappa,ks,srcinfo,detinfo)

ne = size(e,2); % number of edges on the domain boundary
hmatrix = zeros(1,2*ne);
rmatrix = zeros(1,2*ne);

qmatrix = kappa*ones(1,ne);
gmatrix = zeros(1,ne);

% The sources we use are Gaussian sources (to mimic point sources)
xs=srcinfo(1,ks);
ys=srcinfo(2,ks);
srcseg=srcinfo(3,ks);
for k = 1:ne

	x1 = p(1,e(1,k)); % x at first point in segment
	y1 = p(2,e(1,k)); % y at first point in segment
	x2 = p(1,e(2,k)); % x at second point in segment
	y2 = p(2,e(2,k)); % y at second point in segment
	xm = (x1 + x2)/2; % x at segment midpoint
	ym = (y1 + y2)/2; % y at segment midpoint
    
    if detinfo(2,k)==srcseg % if the edge lives on the same side with the source
        gmatrix(k) = 5*exp(-((xm-xs)^2+(ym-ys)^2)/0.01); % sources are Gaussians
    end

end