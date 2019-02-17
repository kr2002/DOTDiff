%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOTDiff: Diffuse optical tomography (DOT) based on 
%          the diffusion model for light propagation
%
%          In the case of steady-state measurement, we 
%          can only reconstruct one of the absorption
%          and diffusion coefficient. This code aims 
%          at reconstructing the absorption coefficient.
%
%          The reconstruction algorithm is based on 
%          the least-square formulation of the inverse
%          problem minimization.
%
% Author:    Kui Ren
% Address:   Math and ICES, UT Austin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The diffusion model:
%
% -\nabla\cdot\gamma \nabla u + \sigma u = 0  in \Omega
% \bnu\cdot\nabla u+ \kappa u = f,  on \partial \Omega
%
% The measurement quantity:
% 
% g= u  on \partial\Omega
%
% The data:
%
% (f_j, g_j), 1\le j\le N_s
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The algorithm: minimize the functional \Phi given by
%
% \Phi_1:= 1/2*\sum_{j=1}^{N_s} \int_{\partial\Omega} (u_j-g_j^*)^2 ds
%
% or
%
% \Phi_2:= 1/2*\sum_{j=1}^{N_s} \int_{\partial\Omega} ((u_j-g_j^*)/g_j^*)^2 ds
%
% Gradient of the objective functions is computed with the adjoint state
% method where the adjoint problems are:
%
%  -\nabla\cdot\gamma \nabla w_j + \sigma w_j =S in \Omega
%  \bnu\cdot\gamma \nabla w_j +\kappa w_j = -r_j
%
%  where r_j:=(u_j-g_j^*) for \Phi_1 ad r_j:=(u_j-g_j^*)/(g_j^*)^2 for \Phi_2  
% 
% The gradient in direction \delta\sigma is given as:
% \Phi' = \int_\Omega u_j w_j \delta\sigma dx
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: 
%
% Step 1: Generate geometry using PDETOOL, and save the data
%
% Step 2: Setup simulations parameters accordingly
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;  clear all; 

tic;
tb=toc;

% Load information on the domain
disp(' ');
disp(' ');
disp('Setting simulation geometry and parameters .......');
disp(' ');

load 'geo-2b2';

MaxIT=200; % max # of iterations allowed for the optimization algorithm
Ns=36; % number of illumination sources used


% Create a Cartesian grid for inversion
dx=0.05; x=0:dx:2;
dy=0.05; y=0:dy:2;
Nx=length(x);
Ny=length(y);
% [X,Y]=meshgrid(x,y);

% Generate regular finite element mesh on rectangular geometry
% The mesh is not necessarily uniform. One can use other mesh here
[P,E,T]=poimesh(geo,Nx-1,Ny-1);
M=Nx*Ny; % total number of nodes in the mesh

srcinfo=SetSources(Ns);
[detinfo srcdetpair]=SetDetectors(P,E,Ns,srcinfo);
Nd=length(E(1,:)); % # of possible detectors, i.e number of boundary nodes, 
%                    also # of boundary edges

% Set gamma, kappa and true sigma
rec1=[0.2 0.4; 0.4 0.9];
rec2=[1.4 1.7; 1.4 1.7];
rec3=[1.0 1.8; 0.3 0.6];
circ1=[0.5 0.5 0.25];
circ2=[1.5 1.5 0.25];
%circ3=[0.5 1.2 0.3];
%circ4=[1.4 1.6 0.2];
%circ5=[1.7 0.4 0.2];
%circ6=[0.5 0.4 0.2];

kappa=1.0;
gamma=zeros(M,1);
sigmat=zeros(M,1);

gamma=0.02*ones(M,1);
%gamma=0.02+0.02*ind_rec(P,rec1)+0.02*ind_rec(P,rec2)+0.02*ind_rec(P,rec3);
r1=(P(1,:)-0.5).^2+(P(2,:)-0.5).^2;
r1=sqrt(r1);
r2=(P(1,:)-1.5).^2+(P(2,:)-1.5).^2;
r2=sqrt(r2);
sigmat=0.1+0.2*cos(pi*r1/0.5)'.*ind_circ(P,circ1)+0.2*cos(pi*r2/0.5)'.*ind_circ(P,circ2);

% Interpolate to Cartesian grid
sigmatg=tri2grid(P,T,sigmat,x,y); % true sigma on Cartesian grid

figure;
pcolor(x,y,sigmatg); axis tight; colorbar('SouthOutside');
axis square; axis off; shading interp;
title('true \sigma');
drawnow;

disp('Finished setting simulation geometry and parameters .......');

% Generating synthetic data
disp(' ');
disp(' ');
disp('Generating synthetic data .......');
disp(' ');

zerosrc=zeros(M,1); % interior source term is 0
measn=zeros(Nd,Ns); % measured noisy data

noiselevel=0.0; % noise level in polluted data
for ks=1:Ns
    
    % Solve the diffusion equation
    [ut meas]=DiffSolve(P,E,T,gamma,sigmat,kappa,zerosrc,ks,Nd,srcinfo,detinfo);
       
    % Plot forward solution
    %utg=tri2grid(P,T,ut,x,y);
    %figure;
    %pcolor(x,y,utg); axis tight; colorbar('SouthOutside');
    %axis square; axis off; shading interp;
    %drawnow;
    %pause
    
    % Add multiplicative noise to data
    noise=noiselevel*2*(rand(1,Nd)-0.5);
	measn(:,ks)=meas.*(1+noise);
        
    disp(['Synthetic data generated for source #: ' num2str(ks)]);
    disp('  ');
    clear ut meas noise;
end

disp('Finished generating synthetic data .......');

% Setup initial guess
disp(' ');
disp(' ');
disp('Setting initial guess .......');
disp(' ');

sigma0=0.1*ones(M,1); % initial guess of sigma
sigma0g=tri2grid(P,T,sigma0,x,y); % interpolate onto Cartesian grid

figure;
pcolor(x,y,sigma0g); axis tight; colorbar('SouthOutside');
axis square; axis off; shading interp;
%caxis([0.0 0.4]);
title('initial guess of \sigma');
drawnow;

X0=sigma0;

disp('Finished setting initial guess .......');

% This short part is only for debugging: checking gradient calculation
%[f0 g0]=DOTDiffObj(X0,x,y,dx,dy,Nx,Ny,P,E,T,gamma,kappa,Ns,Nd,...
%                                BdaryNode,NormVecNode,ds,measn)                        
%g0g=tri2grid(P,T,g0,x,y);
%figure;
%pcolor(x,y,g0g); axis tight; colorbar('SouthOutside');
%axis square; shading interp;
%title('Gradient at initial guess');
%drawnow;

OptimMethod='UNCON'; % use uncontstrained (UNCON) or constrained (CON) minimization

% Setup the minimization algorithm
disp(' ');
disp(' ');
disp('Minimizing objective function .......');
disp(' ');

f=@(X) DOTDiffObj(X,x,y,dx,dy,Nx,Ny,P,E,T,gamma,kappa,Ns,Nd,srcinfo,detinfo,srcdetpair,measn);

if strcmp(OptimMethod,'UNCON')
    options=optimoptions(@fminunc,'Algorithm','quasi-newton', ...
    'Display','iter-detailed','GradObj','on','TolFun',1e-12,...
    'MaxIter',MaxIT);
    [X,fval,exitflag,output,grad]=fminunc(f,X0,options);
else
    % set inequality constraint
    Aieq=zeros(1,M);
    Bieq=0;
    % set equality constraint
    Aeq=zeros(1,M);
    Beq=0;
    % set upper and lower bounds
    LB=0.05*ones(1,M);
    UB=0.4*ones(1,M);

    options=optimoptions(@fmincon,'Algorithm','sqp', ...
    'Display','iter-detailed','GradObj','on','TolFun',1e-12,...
    'MaxIter',MaxIT);
    %options=optimset('Display','iter-detailed','GradObj','on','TolFun',1e-12,'MaxIter',MaxIT);
    %options = optimset('algorithm','sqp','maxfunevals',5000,'maxiter',100);
    %options = optimset(options,'tolx',1e-9,'tolcon',1e-9,'tolfun',1e-6);
    %options = optimset(options,'GradObj','on','GradConstr','off');
    
    [X,fval,exitflag,output,lambda]=fmincon(f,X0,Aieq,Bieq,Aeq,Beq,LB,UB,[],options);
    %[X,fval,exitflag,output]=fmincon(f,X0,zeros(M,M),zeros(M,1),[],[],LB,UB);
end
sigmar=X; % reconstructed sigma

disp(' ');
disp(' ');
disp('Finished minimizing objective function .......');

disp(' ');
disp(' ');
disp('Plotting final results .......');
disp(' ');

% Plot reconstruction results
sigmarg=tri2grid(P,T,sigmar,x,y); % interpolate onto Cartesian grid
figure;
pcolor(x,y,sigmarg); axis tight; colorbar('SouthOutside');
axis square; axis off; shading interp;
title('reconstructed \sigma');
drawnow;

disp('Finished plotting final results .......');

% Save simulation results
save Exp01-info geo P E T srcinfo detinfo MaxIT noiselevel dx dy
save Exp01-result sigmat sigma0 sigmar -ASCII

te=toc;
disp(' ');
disp(' ');
disp(['The code run for: ' num2str(te-tb) ' seconds']);
disp(' ');
disp(' ');

% This last line is used to close MATLAB after the computation. It is 
% only used when runing the code in background.

%exit; % to exit MATLAB 