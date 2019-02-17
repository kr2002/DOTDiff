%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function set up boundary sources for the simulations
%
% srcloc: location of sources (i.e. center)
% srcseg: which part of boundary the source is located at: 
%                      bottom 1, right 2, top 3 and left 4
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function srcinfo=SetSources(Ns)

% The following setup is for the domain (0,1)x(0,1). This can
% be changed for other type of domains.

% Source locations
location=[0.2:0.2:1.8 2*ones(1,9) 1.8:-0.2:0.2 zeros(1,9);...
                zeros(1,9) 0.2:0.2:1.8 2*ones(1,9) 1.8:-0.2:0.2];
% The side of the boundary where the each source is located
segment=[ones(1, 9) 2*ones(1,9) 3*ones(1,9) 4*ones(1,9)]; 

if Ns>length(segment)
	disp('Error in setting # of sources!');
	exit;
end

srcinfo(1:2,:)=location(:,1:Ns);
srcinfo(3,:)=segment(1:Ns);