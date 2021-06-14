function [Rtot,dRaux,Forcing] = DNF_SO_equations(sys,t,Utot,dUtot,~)


%% parameters of the system
mass = sys.parameters.mass;
stiffness = sys.parameters.stiffness;
damping = sys.parameters.damping;
force=sys.parameters.force;

B = sys.parameters.B;
AH = sys.parameters.AH;
C = sys.parameters.C;

n_modes=sys.parameters.n_modes;

%% variables of the system
x = Utot(1:sys.nz/2,:);dx=dUtot(1:sys.nz/2,:);
y = Utot(sys.nz/2+1:sys.nz,:);dy=dUtot(sys.nz/2+1:sys.nz,:);
r = Utot(sys.nz+1:sys.nz+sys.nz_aux,:);
lambda = Utot(end); % unused here because it is the forcing pulsation.



%% Auxiliary equations
Raux = zeros(sys.nz_aux,1);     dRaux = zeros(sys.nz_aux,1);
% The auxiliary variable r is a vector of length  nz_aux=2*N^2
% the first N^2 entries are in the form:
% r(1:N^2)={x1*x1,x1*x2,...,x1*xN, x2*x1,x2*x2,...,x2*xN, ..., xN*x1,xN*x2,...,xN*xN}
% function for picking the p-th element of r that corresponds to x(j)*x(k)
pick_xx=@(j,k) (j-1)*n_modes+k;
% the second N^2 entries are in the form:
% r(N^2+1:2*N^2)={y1*y1,y1*y2,...,y1*yN, y2*y1,y2*y2,...,y2*yN, ..., yN*y1,yN*y2,...,yN*yN}
% function for picking the p-th element of r that corresponds to y(j)*y(k)
pick_yy=@(j,k) n_modes^2+(j-1)*n_modes+k;
for j=1:n_modes
    for k=1:n_modes
        Raux(pick_xx(j,k)) = r(pick_xx(j,k)) - x(j)*x(k);
        Raux(pick_yy(j,k)) = r(pick_yy(j,k)) - y(j)*y(k);
    end
end


%% Physical equations
R = zeros(sys.nz,1);
R(1:sys.nz/2)= y - dx;
R(sys.nz/2+1:sys.nz)= - stiffness*x - lambda*damping*y -  mass*dy;

for p=1:n_modes
    for i=1:n_modes
        for j=1:n_modes
            for k=1:n_modes
                R(sys.nz/2+p) =R(sys.nz/2+p)...
                    -AH(p,i,j,k)*x(i)*r(pick_xx(j,k))...
                    -B(p,i,j,k)*x(i)*r(pick_yy(j,k))...
                    -0*C(p,i,j,k)*r(pick_xx(i,j))*y(k);
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
for kk=sys.nz/2+1:sys.nz
    Forcing(:,kk) = 0*force(kk-sys.nz/2)*cos(t);
end

end
