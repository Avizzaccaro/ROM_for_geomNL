function [parameters]=DNF_SO_import_results_from_code_aster()
%
%
% this function imports results from code aster such as eigenvectors of the
% selected modes, natural frequencies of the system, nonlinear stiffness
% tensors, etc. 
% and generates stiffness, mass, and damping matrices in modal coordinates
% as well as the force vector in modal coordinates
%
%
%
%

% import the tip node from the Mesh folder
force_dof=247*3-2; 	% force on node 247 in direction x
zetam=0;
zetak=3E-7;
%%
% Import eigenvectors matrix
load('../astk/mapping_and_red_dyn.mat');

n_modes=length(Mr);

%% Build modal force vector
F0= zeros(length(phi),1);% in real coordinates
F0(force_dof)=30; % force in newton N 
Fr=phi*F0; % in modal coordinates


%% save parameters for Manlab analysis
parameters.n_modes = n_modes;

parameters.stiffness=Kr;
w=sqrt(diag(Kr));
parameters.mass=Mr;
parameters.damping =zetam*Mr+zetak*Kr;


parameters.AH=AH;
parameters.B=B;
parameters.C=-2*zetak*(AH-H);
for j = 1:n_modes
    for k = 1:n_modes
        parameters.C(:,:,j,k)=parameters.C(:,:,j,k)+...
            (zetam+3*w(j)^2*zetak)*parameters.B(:,:,j,k)+...
            (-zetam+2*w(j)^2*zetak)*(Cs(:,:,j,k)+Cd(:,:,j,k))+...
            (-zetam+2*w(k)^2*zetak)*w(j)/w(k)*(Cs(:,:,j,k)-Cd(:,:,j,k));
    end
end
parameters.force =Fr;
parameters.force_dof=force_dof;

parameters.Phi=phi;
parameters.a=a;
parameters.b=b;
parameters.c=c;

parameters.ad=a*0;
parameters.bd=b*0;
parameters.cd=2*zetak*parameters.a;
for j = 1:n_modes
    for k = 1:n_modes
        parameters.cd(:,j,k)=parameters.cd(:,j,k)+...
            (zetam+3*w(j)^2*zetak)*parameters.b(:,j,k)+...
            (-zetam+2*w(j)^2*zetak)*(Zss(:,j,k)+Zdd(:,j,k))+...
            (-zetam+2*w(k)^2*zetak)*w(j)/w(k)*(Zss(:,j,k)-Zdd(:,j,k));
        parameters.ad(:,j,k)=-parameters.cd(:,j,k)*w(k)^2;
        parameters.bd(:,j,k)=parameters.cd(:,j,k)+(2*zetam+zetak*w(j)^2+zetak*w(k)^2)*parameters.b(:,j,k);
    end
end


end
