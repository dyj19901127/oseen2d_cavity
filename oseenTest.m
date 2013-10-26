clear;
clc;
close all;


logName = [datestr(now,'mmddyyyy') '.cavtst'];
lf = Msgcl(99,logName);
lf.pmsg(lf.ALL,'************************************************************************');
lf.pmsg(lf.ALL,'* Testing IRKA model reduction for the linearized Boussinesq System');
lf.pmsg(lf.ALL,'* over a 1m x 1m 2D cavity with temperature induced flow.');


dataLoc = 'data50_2';
n=50;

%% Build the matrices for the linear model
lf.pmsg(lf.ERR,'Creating the linearized model for u, T, and p');
create_uTp_matrices;


%% Solve the full linear model
Z0 = zeros(a,1);
%Z0(u_nodes) = u_adv(int_nodes);  Z0(v_nodes) = v_adv(int_nodes);
%Z0(T_nodes) = T_eq(T_int_nodes);
tspan = [0,10];

lf.pmsg(lf.ERR,'Solving the full system.')
[ t, Z, Y ] = daesol2( A, B, C, D, E, tspan, @u_fcn,Z0);

lf.pmsg(lf.ERR,'Plotting control.')
figure('name','Control Function');
plot(t,u_fcn(t),'r');

lf.pmsg(lf.ERR,'Plotting velocity vector field movie.');
figure('name','Velocity Vector Field');
u_res = zeros((n+1)^2,1);v_res = zeros((n+1)^2,1);
for k=1:size(Z,2)
  u_res(int_nodes) = Z(u_nodes,k);
  v_res(int_nodes) = Z(v_nodes,k);
  quiver(x(:,1),x(:,2),u_res,v_res);
  F(k) = getframe;
end;

%movie(F,1);

%% Create reduced order model
lf.pmsg(lf.ERR,'Creating reduced order model.')
% Set the test run parameters
loglevel = 10;
tol = [1e-4, 1e-8];
distr_type = 1;
distr = [-1,2,0,0];
testnum = 3;
r = 15;

[Ar,Br,Cr,Dr,Er,Vr,Wr,siter] = irka_dae2(A,B,C,D,E,nv,r,loglevel,tol,distr_type,distr);

lf.pmsg(lf.WARN,'Printing Bode plot of full system versus reduced model')
figure('name', 'Bode Plots', 'NumberTitle', 'off');
bode_plt(A,B,C,D,E,-4,10,100);
hold on
bode_plt(Ar,Br,Cr,Dr,Er,-4,10,100, '-ro');
hold off
legend('G','Reduced G')

lf.pmsg(lf.ALL,'* END of test. Exiting.');
lf.pmsg(lf.ALL,'************************************************************************');

%% Solve the reduced model
xr0 = zeros(r,1);
[ t, xr, yr ] = daesol2( Ar, Br, Cr, Dr, Er, tspan, @u_fcn,xr0);
xr2f = Vr*xr;

avh = figure('name', 'Average Vorticity');
plot(t,Y)
plot(t,yr,'ko');
legend('Full Model','Reduced Model','Location','SE');



figure('name','Velocity Vector Field (Reduced)');
u_res = zeros((n+1)^2,1);v_res = zeros((n+1)^2,1);
for k=1:size(Z,2)
  u_res(int_nodes) = xr2f(u_nodes,k);
  v_res(int_nodes) = xr2f(v_nodes,k);
  quiver(x(:,1),x(:,2),u_res,v_res);
  Fr(k) = getframe;
end;

figure('name','State Error Plot');
errM = (xr - Z(1:nv,:)).^2;
err = sqrt(sum(errM));
plot(t,err)

