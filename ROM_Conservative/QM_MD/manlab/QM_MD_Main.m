clear all
clc

global U Section Diagram    % Global variables to export point from the diagram. 

addpath('/home/av3116/work/MANLAB-4.0/SRC')  % Manlab SRC path


% Parameters of the system: 
parameters=QM_MD_import_results_from_code_aster();


%% Manlab variables
L=parameters.n_modes;
nz= 3*L;     % number of main equations of the DAE system
nz_aux =3*L^2;
H =15;     % number of harmonics


%% initialization of the system
sys=SystHBQ(nz,nz_aux,H,@QM_MD_equations,@point_display,@global_display,parameters,'autonomous');

%sys.zi_phase=1;
% starting point
omega=sqrt(parameters.stiffness(1,1));%  before resonance
lambda=0;
%
Z0  =rand(2*H+1,sys.nz_tot)*1e-11;
% Initial guess for backbone
% modal amplitude of mode 1 set to small:
q0=1e-4;
% cosine of Z0 set to equal to x0
Z0(2,1)=q0;
% velocities:
Z0(2+H,sys.nz/3+1) = -omega*q0;
% accelerations:
Z0(2,2*sys.nz/3+1) = -omega^2*q0;

U0tot = sys.init_U0(Z0,omega,lambda);

% The bifurcation diagram represents
% first cosine of x with respect to lambda and
% first sine of x with respect to lambda :
dispvars = [sys.neq 2; ... %blue curve is cos
    sys.neq 2+H];%red curve is sin
    %sys.neq sys.neq+1];% green is lambda 


%% Launch of Manlab with options
Manlab('sys'       ,sys , ...
    'U0value'         ,U0tot, ...
    'ANMthreshold'    ,1e-7 , ...
    'NRthreshold'     ,1e-7, ...
    'StabilityCheck'  ,0, ...
    'displayvariables',dispvars);     % MANLAB run


%% Remember to save both the diagram (from Manlab GUI) 
%% and sys.mat (from your workspace) 
%% or you cannot postprocess
