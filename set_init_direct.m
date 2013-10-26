function [ b, c ] = set_init_direct( A, B, C, D, E, sigma, lf )
%SET_INIT_DIRECT Creates the initial set of direction vectors for IRKA
%
%   We start with the following system:
%
%     E(x_dot) = Ax + Bu
%            y = Cx + Du
%
%   Inputs:
%     A,B,C,D,E - Matrices from the system above
%         sigma - sample points
%            lf - log file object
%
%   Outputs:
%     b - right tangential vectors
%     c - left tangential vectors
%
%--------------------------------------------------------------------------
% Author: Alan Lattimer
% email: alattime at vt dot edu
% Date: April 2013
%

uselog = 0;
if nargin < 7 
  fprintf('No log file\n');
else
  uselog=1;
  lf.pmsg(lf.WARN, '=> set_init_direct.m');
end


% Initialize 
if uselog 
  lf.pmsg(lf.PED,'  Initializing the parameters used.');
end
r = length(sigma);
m = size(B,2);  p = size(C,1);
b = zeros(m,r); 
c = zeros(p,r);

if uselog
  lf.pmsg(lf.PED,'  Parameter data as follows:');
  lf.pmsg(lf.PED,'    => r = %d',r);
  lf.pmsg(lf.PED,'    => m = %d',m);
  lf.pmsg(lf.PED,'    => p = %d',p);
end
  

% Assigning direction vectors
if uselog
  lf.pmsg(lf.WARN,'  Setting initial directions.');
end
if m==1 && p==1
  lf.pmsg(lf.PED,'  SISO system, set directions to 1.');
  b = ones(1,r);
  c = ones(1,r);
else
  lf.pmsg(lf.PED,'  MIMO system, directions are the right and left vectors of SVD')
  for k=1:r
    [U,~,V] = svd(C*(((sigma(k)*E)-A)\B));

    if p < m
      b(:,k) = V(:,1);
      c(:,k) = U(:,1);
%       b(:,k) = mean(V,2);
%       c(:,k) = mean(U,2);
    else
      b(:,k) = U(:,1);
      c(:,k) = V(:,1);
%       b(:,k) = mean(U,2);
%       c(:,k) = mean(V,2);
    end
    
    if uselog
      lf.pmsg(lf.PED,'  b_%d and c_%d have been set.',k,k);
    end
  end
end

if uselog
  lf.pmsg(lf.WARN,'<= set_init_direct.m')
end

end

