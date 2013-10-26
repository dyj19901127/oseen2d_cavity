savedata = 0;

% If  oseen_2d_active_matrices.m is not in your path, add it here
%path('/Path/To/oseen_2d_active_matrices',path);


% Set data location if not already set
if (~exist('dataLoc','var'))
  dataLoc = 'data';
  n = 10;
end

dataSavePath = '/Users/alattime/Documents/Projects/MR_Research/ModRed/irka_dae/';
if ~exist(dataSavePath,'dir')
  dataSaveFile = dataLoc;
else
  dataSaveFile = fullfile(dataSavePath, dataLoc);
end

%% Set required constants needed to linearize the flow.

Tref = 300;
dir_u = [0,0];
dir_T = [-0.5 0.5];
epsilon = 0;
th = 1e-16;

Re   = 100;
Pr   = 0.7;
Ra   = 10;

material.gx   = 0;
material.gy   =-1;
material.epsilon = epsilon;
material.f_function = @f_zero;

material.mu               = 1/Re;      % fluid viscosity
material.kappa            = 1/Re/Pr;   % fluid conductivity
material.beta             = Ra;        % bouyancy term


%% Import the data from the raw OpenFOAM files
lf.pmsg(lf.ERR,'Loading the OpenFOAM data');
[u_adv,v_adv,T_eq,idx_u,idx_T,idx_p,x,e_conn,~] = importFoamData(n,dataLoc,Tref,lf);

int_nodes = find(idx_u(:,1)>0);
u_nodes = idx_u(int_nodes,1);
v_nodes = idx_u(int_nodes,2);
T_int_nodes = find(idx_T>0);
T_nodes = idx_T(T_int_nodes);





%% Create the linear model
lf.pmsg(lf.ERR,'Building Oseen matrices...')
[A, B, M, C] = oseen_2d_active_matrices(x,e_conn,                               ...
                                        idx_u, idx_T, idx_p,                    ...
                                        dir_u, dir_T,                           ...
                                        material,                               ...
                                        u_adv, v_adv, T_eq                        );

lf.pmsg(lf.ERR,'completed.');

%% Remove singularity introduced because of the constant pressure
lf.pmsg(lf.ERR,'Break up A in order to remove singularity in A12 and assemble');
lf.pmsg(lf.ERR,'needed matrices.')
lf.pmsg(lf.PED,'== Getting matrix sizes')
a = size(A,1)-1; m = size(M,1);
p = a-m;
lf.pmsg(lf.PED,'== Splitting A');
A11 = A(1:m,1:m);
A12 = A(1:m,m+1:end);

lf.pmsg(lf.PED,'== Removing noise from A')
A11(abs(A11)<th)=0;
A12(abs(A12)<th)=0;
lf.pmsg(lf.PED,'== Resizing A12 to remove singularity.')
A12=A12(:,1:p);
lf.pmsg(lf.PED,'== Reassembling A.');
A=[A11,A12;A12',zeros(p)];

lf.pmsg(lf.PED,'== Building E.');
E = sparse(a,a);
E(1:m,1:m) = M;
lf.pmsg(lf.PED,'== Building B, C, D.');
D = 0;
B = B(1:a);
C = [C',zeros(1,p)];
nv = m;
lf.pmsg(lf.ERR,'Matrices constructed and singularity removed.')


%% Save data files
if savedata 
  lf.pmsg(lf.WARN,'Saving data to %s ... ',dataSaveFile);
  save(dataSaveFile,'A','B','C','D','E', 'nv');
  lf.pmsg(lf.WARN,'completed.')
else
  lf.pmsg(lf.WARN, 'Not saving the data')
end;

%break;



