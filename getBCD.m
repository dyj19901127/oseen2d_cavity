function [B,C,Dk,Dkp] = getBCD(A11, A12, B1, B2, C1, C2, D, E11)
%GETBCD Summary of this function goes here
%   Detailed explanation goes here

% **************************************************************************
% = Function: getBCD
% = Inputs:
% = Outputs:
% Code from Serkan Geugercin
%


  % Compute B, C, D matrices
  nv = size(A11,1);
  np = size(A12,2);
  nb = size(B1,2);
  nc = size(C1,1); 
  K  = [ E11 , A12;
         A12', zeros(np,np)];
  X  = K \ [ zeros(nv,nb); B2 ];
  B  = B1 - A11*X(1:nv,:);
  Dk = D  - C1 *X(1:nv,:);

  if( norm(C2,1) > 0 )
      % C2 is not zero
      nc  = size(C2,1);
      X   = K \ [ zeros(nv,nc); C2' ];
      C   = C1 - X(1:nv,:)'*A11;
      Dk  = Dk - X(1:nv,:)'*B;
      Dkp = X(nv+1:nv+np,:)'*B2;
  else
      C   = C1;
      Dkp = sparse(size(D,1),size(B1,2));
  end
  
end



