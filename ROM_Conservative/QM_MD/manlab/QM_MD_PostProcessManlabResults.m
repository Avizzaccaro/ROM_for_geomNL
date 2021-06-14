%% Postprocessing
%% Plot FRF or NNM at tip node where force is applied
clear all
clc
addpath('/home/av3116/work/MANLAB-4.0/SRC')  % Manlab SRC path


load('diagram.mat');
load('sys.mat');
sys.parameters.Phi=sys.parameters.Phi';

pts=length(diagram.ChckPoint);
Upp = diagram.ChckPoint{1,1}.Upp;
H=sys.H;
T=2*pi;
temps=0:T/(20*H+1):T;
nt=length(temps);
VectCos=cos(temps'*(1:H));
VectSin=sin(temps'*(1:H));
MAT=[ ones(nt,1) , VectCos , VectSin ];  % MAT(nt,2*H+1)
L=sys.nz/3;
nz_tot=sys.nz_tot;
nb_pt = size(Upp,2);


ampl=zeros(nb_pt*pts,1);
dof=sys.parameters.force_dof;
omega_nl=zeros(nb_pt*pts,1);


R=zeros(nt,L);
S=zeros(nt,L);

X_dof=zeros(1,nt);


for iChunk=1:pts
    Upp = diagram.ChckPoint{1,iChunk}.Upp;
    omega=zeros(nb_pt);
    lambda=omega;
    % Loop on the point representation of the series
    for iPoint = 1:nb_pt
        iCont=nb_pt*(iChunk-1)+iPoint;        
        
        Uk = Upp(:,iPoint);
        [Zk,omega(iPoint),lambda(iPoint)] = sys.get_Ztot(Uk);
        Ut=MAT*Zk;  % Ut(nt,nz)
        for p=1:L
            R(:,p)=Ut(:,p);
            S(:,p)=Ut(:,L+p);
        end
        
        X_dof=X_dof*0;
        for i=1:L
            X_dof=X_dof+sys.parameters.Phi(dof,i)*R(:,i)';
            %
            for j=1:L
                X_dof=X_dof+sys.parameters.theta(dof,i,j)*(R(:,i).*R(:,j))';
            end
        end
        
        omega_nl(iCont)=omega(iPoint);
        ampl(iCont)=(max(X_dof(:))-min(X_dof(:)))/2;

        
    end
end


%


figure(1)
plot(omega_nl,ampl,'Color',[0.6,0.3,.0],'LineWidth',2);
xlabel('\omega [rad/s]')
ylabel('A [m]')
%
legend('QM MD')

