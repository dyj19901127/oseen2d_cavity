clear;
clc;
%close all;


dataLoc = 'data';
Tref = 300;
n = 10;
dir_u = [0,0];
Re = 10;
epsilon = 0;
th = 1e-16;



[u_adv,v_adv,~,idx_u,~,idx_p,x,e_conn,~] = importFoamData(n,dataLoc,Tref);

int_nodes = find(idx_u(:,1)>0);
u_nodes = idx_u(int_nodes,1);
v_nodes = idx_u(int_nodes,2);





[A, B, M] = oseen_2d_matrices(x,e_conn,                           ...
                              idx_u, idx_p, dir_u,                ...
                              Re, epsilon,                        ...
                              @f_fcn,                             ...
                              u_adv, v_adv                          );

 

%break;

% Break up A in order to remove singularity in A12
a = size(A,1)-1; m = size(M,1);
p = a-m;
A11 = A(1:m,1:m);
A12 = A(1:m,m+1:end);

A11(A11<th)=0;
A12(A12<th)=0;
A12=A12(:,1:p);

A=[A11,A12;A12',zeros(p)];

E = sparse(a,a);
E(1:m,1:m) = M;
D = 0;
B = B(1:a);
C = [(1/m)*ones(1,m),zeros(1,p)];
nv = m;



save(dataLoc,'A','B','C','D','E', 'nv');
break

M_epsilon = 5;
M_tilde = [M zeros(m,a-m); zeros(a-m,m),M_epsilon*eye(a-m)];

C_tilde = [(1/m)*ones(1,m),zeros(1,a-m)];

opt = odeset('Mass', M_tilde, 'RelTol', 1.e-04, 'AbsTol', 1.e-06');
% odefcn = @(t,z) A*z + B*sin(pi*t);
odefcn = @(t,z) A*z + B*-1;

z0 = zeros(a,1);
z0(u_nodes) = u_adv(int_nodes);  z0(v_nodes) = v_adv(int_nodes);

t = 0:0.1:10;

[T, Z] = ode23(odefcn, t, z0, opt);
figure
u_res = zeros(121,1);v_res = zeros(121,1);
for k=1:size(Z,1)
  u_res(int_nodes) = Z(k,u_nodes);
  v_res(int_nodes) = Z(k,v_nodes);
  quiver(x(:,1),x(:,2),u_res,v_res);
  F(k) = getframe;
end;

movie(F,1);
