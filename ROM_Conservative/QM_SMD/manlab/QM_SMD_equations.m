function [Rtot,dRaux,Forcing] = QM_SMD_equations(sys,t,Utot,dUtot,~)
 

%% parameters of the system
mass = sys.parameters.mass;
stiffness = sys.parameters.stiffness;
damping = sys.parameters.damping;
force=sys.parameters.force;

GD=sys.parameters.GD;
GV=sys.parameters.GV;
GA=sys.parameters.GA;
GVd=sys.parameters.GVd;
HD=sys.parameters.HD;
HV=sys.parameters.HV;
HA=sys.parameters.HA;
HVd=sys.parameters.HVd;

n_modes=sys.parameters.n_modes;


%% variables of the system
% displacement x, lenght L
x = Utot(1:sys.nz/3,:);dx=dUtot(1:sys.nz/3,:);
% velocity y, lenght L
y = Utot(sys.nz/3+1:2*sys.nz/3,:);dy=dUtot(sys.nz/3+1:2*sys.nz/3,:);
% acceleration z, lenght L
z = Utot(2*sys.nz/3+1:sys.nz,:);%dz=dUtot(2*sys.nz/3+1:sys.nz,:);
% auxiliry variables, lenght 3*L^2
r = Utot(sys.nz+1:sys.nz+sys.nz_aux,:);
% 
lambda = Utot(end); % unused here because it is the forcing pulsation, used for the NNM



%% Auxiliary equations
Raux = zeros(sys.nz_aux,1);     dRaux = zeros(sys.nz_aux,1);
% The auxiliary variable r is a vector of length  nz_aux=3*L^2  
pick=@(j,k) (j-1)*n_modes+k;
% first L^2 rows of r are x^2
% second L^2 rows of r are y^2
% third L^2 rows of r are x*z
for j=1:n_modes
    for k=1:n_modes
        Raux(pick(j,k)) = r(pick(j,k)) - x(j)*x(k);
        Raux(sys.nz_aux/3+pick(j,k)) = r(sys.nz_aux/3+pick(j,k)) - y(j)*y(k);
        Raux(2*sys.nz_aux/3+pick(j,k)) = r(2*sys.nz_aux/3+pick(j,k)) - x(j)*z(k);
    end
end


%% Physical equations
%  NB the damping is present twice:
% in line 55 is the linear damping (to be multiplied by lambda when computing the NNM)
% in line 69 is the nonlinear damping (just put C to zero when computing the NNM)
R = zeros(sys.nz,1);
R(1:sys.nz/3)= y - dx;
R(sys.nz/3+1:2*sys.nz/3)= z - dy;
R(2*sys.nz/3+1:sys.nz)= - stiffness*x - damping*y -  mass*z;
for p=1:n_modes
for i=1:n_modes
    for j=1:n_modes
    	R(2*sys.nz/3+p) = R(2*sys.nz/3+p)...
    		- GD(p,i,j)*r(pick(i,j))...
    		- GV(p,i,j)*r(sys.nz_aux/3+pick(i,j))...
    		-GVd(p,i,j)*y(i)*x(j)...
    		- GA(p,i,j)*z(i)*x(j);
        for k=1:n_modes
            R(2*sys.nz/3+p) = R(2*sys.nz/3+p)...
                - HD(p,i,j,k)*x(i)*r(pick(j,k))...		%	HD, cubic term in displacement
                - HV(p,i,j,k)*r(sys.nz_aux/3+pick(i,j))*x(k)...	%	HV, cubic term in velocity
                - HVd(p,i,j,k)*y(i)*r(pick(j,k))...		%	HV, cubic term in velocity (damp)
                - HA(p,i,j,k)*z(i)*r(pick(j,k));		%	HA, cubic term in acceleration
    
        end
    end
end
end
%% All the equations of the quadratic system
Rtot = [R;Raux];

%% Forcing terms
% Should be written as if the forcing pulsation value is 1
% i.e. the forcing period is 2*pi
Forcing = zeros(2*sys.H+1,sys.nz_tot);
for kk=2*sys.nz/3+1:sys.nz
    Forcing(:,kk) = force(kk-2*sys.nz/3)*cos(t);
end

end
