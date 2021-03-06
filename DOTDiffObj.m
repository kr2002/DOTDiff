%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The objective functional for the minimization problem
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f g]=DOTDiffObj(X,x,y,dx,dy,Nx,Ny,P,E,T,gamma,kappa,...
                            Ns,Nd,srcinfo,detinfo,srcdetpair,meas)

M=Nx*Ny;
sigmac=X; % current value of sigma

f=0.0;
g=zeros(M,1);

zerosrc=zeros(M,1);
for ks=1:Ns
    
    pred=zeros(1,Nd); % predicted data on measurement locations
    rz=zeros(1,Nd); % residual on measurement locations
 
    [uc pred]=DiffSolve(P,E,T,gamma,sigmac,kappa,zerosrc,ks,Nd,srcinfo,detinfo);
        
    %ucg=tri2grid(P,T,uc,x,y);
    %figure;
    %pcolor(x,y,ucg); axis tight; colorbar('SouthOutside');
    %axis square; axis off; shading interp;
    %drawnow;
    
    measR=meas(:,ks)';
    
    %rz=(pred-measR).*srcdetpair(ks,:);  % residual for \Phi_2
    rz=(pred-measR)./measR.*srcdetpair(ks,:); % residual for \Phi_1

    % calculate the objective function (the part for source ks)
    f=f+0.5*sum(rz.^2.*detinfo(5,:)); % rectangular rule for boundary integral
    
    if nargout > 1 % calculate gradient
        
        % solve the adjoint diffusion problem
        
        %adjsrc=-rz; % source for adjoint problem
        adjsrc=-rz./measR; % source for adjoint problem if normalized obj is used
        
        wc=DiffSolveAdj(P,E,T,gamma,sigmac,kappa,adjsrc,ks,srcdetpair);
        % calculate the gradient w.r.t sigma
        g=g+uc.*wc*dx*dy;
 
    end
    
end

beta=1e-16; % the regularization parameter, change accoring to noise level

% Add regularization to the objective function
[Gx,Gy] = pdegrad(P,T,sigmac);
Gx1=pdeprtni(P,T,Gx); Gy1=pdeprtni(P,T,Gy);
f=f+0.5*beta*sum(Gx1.^2+Gy1.^2)*dx*dy;

% Add regularization term to the gradient
if nargout >1   
    % The Laplacian part
    [Gxx, Gxy]=pdegrad(P,T,Gx1); [Gyx, Gyy]=pdegrad(P,T,Gy1);
    Gx2=pdeprtni(P,T,Gxx); Gy2=pdeprtni(P,T,Gyy);
    DeltaGamma=Gx2+Gy2;
    g=g-beta*DeltaGamma*dx*dy;
    
    % The boundary integral part: this part should be removed if boundary 
    %      value of n is not to be reconstructed
    dndn=zeros(1,Nd);
    for j=1:Nd
        dndn(j)=Gx1(detinfo(1,j))*detinfo(3,j)+Gy1(detinfo(1,j))*detinfo(4,j);
    end
    g=g+beta*sum(dndn.*detinfo(5,:));

end