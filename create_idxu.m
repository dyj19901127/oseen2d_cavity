function [u_adv,v_adv,idx_u] = create_idxu(u, v, bnd_incl)
%CREATE_IDXU creates U and idx_u for oseen_2d_matrices
%
% Inputs:
%   u = data in the x direction
%   v = data in the y derection
%   bnd_incl =  1 -> boundary nodes are included in u,v
%               2 -> boundary nodes not included in u,v
%
% Ouputs:
%   u_adv,v_adv = averaged data over the time steps.
%   idx_u = indices for velocity vector
%

u_adv = mean(u,2);
v_adv = mean(v,2);

if nargin == 2
  bnd_incl = 0;
end

n = round(sqrt(length(u(:,1))));


if bnd_incl==1 
  n=n-2;
end;

bnd_vec = -1*ones(n+2,1);
% Build matrices
ui = [bnd_vec [bnd_vec(1:n)'; ...
               flipud(reshape(1:2:2*(n)^2,n,n)'); ...
               bnd_vec(1:n)'] bnd_vec];
             
vi = [2*bnd_vec [2*bnd_vec(1:n)'; ...
                 flipud(reshape(2:2:2*(n^2),n,n)'); ...
                 2*bnd_vec(1:n)'] 2*bnd_vec];

%                % Build matrices
% ui = [bnd_vec [bnd_vec(1:n)'; ...
%                flipud(reshape(1:n^2,n,n)'); ...
%                bnd_vec(1:n)'] bnd_vec];
%              
% vi = [2*bnd_vec [2*bnd_vec(1:n)'; ...
%                  flipud(reshape((n^2)+1:2*(n^2),n,n)'); ...
%                  2*bnd_vec(1:n)'] 2*bnd_vec];
% 

% Convert the matrices to vectors
ui = reshape( flipud(ui)',[(n+2)^2,1]);
vi = reshape( flipud(vi)',[(n+2)^2,1]);

idx_u = [ui,vi];


end

