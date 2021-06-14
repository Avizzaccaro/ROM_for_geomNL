function [parameters]=QM_MD_import_results_from_code_aster()
%
%
% this function imports results from code aster such as eigenvectors of the
% selected modes, natural frequencies of the system, nonlinear stiffness
% tensors, etc. 
% and generates stiffness, mass, and damping matrices in modal coordinates
% as well as the force vector (200N middle node) in modal coordinates
%
%
%
%

% import the tip node from the Mesh folder
force_dof=247*3-2;
zetam=0;
zetak=3E-7;
%%
% Import eigenvectors matrix
load('../astk/mapping_and_red_dyn.mat');

n_modes=length(Mr);

%% Build modal force vector
F0= zeros(length(phi),1);% in real coordinates
F0(force_dof)=30; % force in Newton
Fr=phi*F0; % in modal coordinates


%% save parameters for Manlab analysis
parameters.n_modes = n_modes;

parameters.stiffness=Kr;
parameters.mass=Mr;
parameters.damping =zetam*Mr+zetak*Kr;

% quadratic conservative
parameters.GD=thKphi+1/2*phiKth+phiGphiphi; 	% on Ri Rj
parameters.GV=phiMth;				% on Si Sj
parameters.GA=phiMth+thMphi;			% on \dot{S}i Rj
% quadratic non-conservative
parameters.GVd=zetam*(thMphi+phiMth)+zetak*(thKphi+phiKth); % on Si Rj 

% cubic conservative
parameters.HD=AH+1/2*thKth+thGphiphi;		% on Ri Rj Rk
parameters.HV=thMth;				% on Si Sj Rk
parameters.HA=thMth;				% on \dos{S}i Rj Rk
% cubic non-conservative
parameters.HVd=zetam*thMth+zetak*thKth;		% on Si Rj Rk

parameters.force =Fr;
parameters.force_dof=force_dof;

parameters.Phi=phi;
parameters.theta=th;





end
